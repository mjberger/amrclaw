c
c ---------------------------------------------------------------------
c
       subroutine SRD_cellMerge(q,nvar,irr,mitot,mjtot,qx,qy,lstgrd,
     .                      dx,dy,lwidth,xlow,ylow,istage,
     .                      ncount,numHoods,mptr,ffluxlen,gfluxlen)

       use amr_module
       implicit double precision (a-h, o-z)

       include "cuserdt.i"
       dimension q(nvar,mitot,mjtot),  irr(mitot,mjtot)
       dimension qx(nvar,mitot,mjtot), qy(nvar,mitot,mjtot)
       dimension gradmx(nvar,mitot,mjtot), gradmy(nvar,mitot,mjtot)
       dimension ffluxlen(mitot+1,mjtot+1),gfluxlen(mitot+1,mjtot+1)

       dimension valnew(nvar,mitot,mjtot), volMerge(mitot,mjtot)
       dimension qMerge(nvar,mitot,mjtot), numHoods(mitot,mjtot)
       dimension nCount(mitot,mjtot)
       dimension xcentMerge(mitot,mjtot),ycentMerge(mitot,mjtot)
       dimension fakeState(nvar), qm(nvar), rhs(5,nvar)
       dimension a(5,5),b(5),db(nvar)
       dimension nborList(25,2)
       character ch

       logical IS_OUTSIDE, NOT_OK_GHOST
       logical quad, nolimiter,verbose
       common /order2/ ssw, quad, nolimiter
       integer omp_get_max_threads
       integer maxthreads/1/


c  xlow,ylow   refers to the corner of grid (w/o ghost cells)
c  xlower,ylower refers to the domain

c  OLD WAY
c  istage will determine how many ghost cells can be trusted:
c  all in 1st stage, 2 less in 2nd stage

c NEW WAY - do one stage at a time then copy in ghost cells at intermediate
c stages for neighboring grids.  With 4 ghost cells, this means that
c interior cells can always trust 2 cells to the side after 1 stage
c never use cells from exterior to the domain for SRD

      IS_OUTSIDE(x,y) = (x .lt. xlower .or. x .gt. xupper .or.
     .                   y .lt. ylower .or. y .gt. yupper)
      
      NOT_OK_GHOST(i,j) = (i .lt. 3 .or. 
     .                     i .gt. mitot-2 .or.
     .                     j .lt. 3 .or. 
     .                     j .gt. mjtot-2)

c :::::::::::::;:::
c
c   try cell merging and reconstruction to stabilize updates in small cut cells
c   calculate all updates using provisional values, then actually do the update
c   (as in, jacobi rather than gauss seidel)
c ::::::::::::::::
      !!areaMin = 2.d0*ar(lstgrd)  
      areaMin = 0.5d0*ar(lstgrd)  
      !!areaMin = ar(lstgrd)  

c
       verbose = .false.
       !verbose = .true.
       !nterms = 5
       nterms = 2

c      some initializations
       ar(lstgrd) = dx*dy   ! area of regular grid cell 
       qMerge   = 0.d0
       dx2 = 2.d0*dx
       dy2 = 2.d0*dy

       ! put in 'normal' values to prevent errors e.g. in converting to prim vars
       ! these fake vals are in conserved variables
       fakeState(1) = 1.4d0
       fakeState(2) = 0.d0
       fakeState(3) = 0.d0
       fakeState(4) = 2.5d0 !(p/gm1, corresponding to fakestate in method)

c     first make neighborhoods - need count for each cells, and width (nhood above)
c     nCount is size of neighborhood, numHoods is number of merged nhoods each cells is in
      call makeNHood(volMerge,xcentMerge,ycentMerge,ncount,irr,numHoods,
     .               mitot,mjtot,lwidth,lstgrd,xlow,ylow,dx,dy,istage,
     .               mptr,ffluxlen,gfluxlen,areaMin)

       if (verbose) then
          totmass =  bigconck(q,irr,mitot,mjtot,lwidth,nvar)
          write(*,911) totmass
 911      format(/,"         mass before redistribution is ",e25.15)
       endif

c       form qMerge vals 
        do 10 j = lwidth-1, mjtot-lwidth+1
        do 10 i = lwidth-1, mitot-lwidth+1
            k = irr(i,j)
            if (k .eq. -1) go to 10 ! no solid cells
            call getCellCentroid(lstgrd,i,j,xc,yc,xlow,ylow,dx,dy,k)
            if (IS_OUTSIDE(xc,yc) .or. NOT_OK_GHOST(i,j)) then 
              qMerge(:,i,j) = rinfinity ! to make sure we dont use it
              go to 10 
            endif
            if (k.eq.lstgrd) then ! full cell is its own nhood
              qMerge(:,i,j) = q(:,i,j)
              go to 10
            endif
            if (ar(k) .gt. areaMin) then
              qMerge(:,i,j) = q(:,i,j)  ! cut cell is large enough
              go to 10
            endif
c
c           this routine currently assumes only one cell needed for
c           merging to make sufficiently large cell, saved in svi,svj
c
            ioff = svi(k)
            joff = svj(k)           
            koff = irr(i+ioff,j+joff)
            call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
     &                           xlow,ylow,dx,dy,koff)

            ! add the one other nbor cell and yourself to
            ! get merged val. volmerge should have same weighting
            qMerge(:,i,j) = (ar(k)*q(:,i,j)/numHoods(i,j) +
     .            ar(koff)*q(:,i+ioff,j+joff)/numHoods(i+ioff,j+joff))/
     .            volMerge(i,j)
 10     continue

        ! gradient of merging neighborhoods, initialized to 0. set using neighboring tiles
        if (ssw .eq. 0.d0) then
           gradmx = 0.d0
           gradmy = 0.d0
           go to 35   ! ssw = 0 means no slopes
        else
           gradmx = rinfinity   ! 0.d0 for debugging
           gradmy = rinfinity   ! 0.d0
        endif

        !do 20 j = 3, mjtot-2
        !do 20 i = 3, mitot-2
        jst  = max(2,lwidth/2)  ! cant go below 2
        jend = min(mjtot-1,mjtot-lwidth/2)
        ist  = max(2,lwidth/2)
        iend = min(mitot-1,mitot-lwidth/2)
        do 20 j = jst, jend
        do 20 i = ist, iend
            k = irr(i,j)
            if (k .eq. -1) go to 20 ! solid cells have no gradient
            if (k .eq. lstgrd) go to 20 ! wont need gradient for full cells
            if (ar(k) .gt. areaMin) go to 20 ! wont need gradient
            call getCellCentroid(lstgrd,i,j,xc,yc,
     &                           xlow,ylow,dx,dy,k)
            if (IS_OUTSIDE(xc,yc)) go to 20  ! exterior cells dont contribute

!--            if (ncount(i,j) .eq. 0) then 
!--            ! these should be all regular cells that are left
!--               nco = 1 ! need to use some kind of nhood to make gradient
!--            else 
!--               nco = ncount(i,j)
!--            endif

            !!!nco = 1 ! this means use 3 by 3 nhood to compute gradients of merging tiles
            nco = 2 ! this means use 5 by 5 nhood to compute gradients of merging tiles
            nborCount = 0 ! if not enough cant use quad, drop to linear

            ! can reconstruct high order poly here, but for now do linear reconstruction
            ! solving least squares problem with all neighboring merged vals. 
            ! Should preprocess the matrix and save
            rhs = 0.d0 ! initialize for accumulation
            a = 0.d0
            x0 = xcentMerge(i,j)
            y0 = ycentMerge(i,j)
            do 22 joff = -nco, nco
            do 22 ioff = -nco, nco
               if (ioff .eq. 0 .and. joff .eq. 0) go to 22 ! no eqn to solve
               koff = irr(i+ioff,j+joff)
               if (koff .eq. -1) go to 22

               call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
     &                              xlow,ylow,dx,dy,koff)
               if (IS_OUTSIDE(xcn,ycn)) go to 22

               nborCount = nborCount + 1
               nborList(nborCount,1) = i+ioff
               nborList(nborCount,2) = j+joff
               deltax = xcentMerge(i+ioff,j+joff) - x0
               deltay = ycentMerge(i+ioff,j+joff) - y0
               a(1,1) = a(1,1) + deltax*deltax
               a(1,2) = a(1,2) + deltax*deltay
               a(2,2) = a(2,2) + deltay*deltay
               rhs(1,:) = rhs(1,:) + deltax * 
     .                    (qMerge(:,i+ioff,j+joff) - qMerge(:,i,j))
               rhs(2,:) = rhs(2,:) + deltay * 
     .                    (qMerge(:,i+ioff,j+joff) - qMerge(:,i,j))
               if (nterms .eq. 5) then
                  a(1,3) = a(1,3) + deltax*deltax**2
                  a(1,4) = a(1,4) + deltax*deltax*deltay
                  a(1,5) = a(1,5) + deltax*deltay**2
                  a(2,3) = a(2,3) + deltay*deltax**2
                  a(2,4) = a(2,4) + deltay*deltax*deltay
                  a(2,5) = a(2,5) + deltay*deltay**2
                  a(3,3) = a(3,3) + deltax**2*deltax**2
                  a(3,4) = a(3,4) + deltax**2*deltax*deltay
                  a(3,5) = a(3,5) + deltax**2*deltay**2
                  a(4,4) = a(4,4) + deltax*deltay * deltax*deltay
                  a(4,5) = a(4,5) + deltax*deltay * deltay**2
                  a(5,5) = a(5,5) + deltay**4
                  rhs(3,:) = rhs(3,:) + deltax**2*
     .                      (qMerge(:,i+ioff,j+joff) - qMerge(:,i,j))
                  rhs(4,:) = rhs(4,:) + deltax*deltay*
     .                      (qMerge(:,i+ioff,j+joff) - qMerge(:,i,j))
                  rhs(5,:) = rhs(5,:) + deltay**2*
     .                      (qMerge(:,i+ioff,j+joff) - qMerge(:,i,j))
               endif
 22          continue

             if (ghost_ccg) then
                call getBndryInfo(alf,beta,k,bxpt,bypt)
                distx = x0 - bxpt
                disty = y0 - bypt
                d_dot_n = distx*alf + disty*beta
                dNormx = d_dot_n * alf
                dNormy = d_dot_n * beta
                dnx = 2.d0*dNormx
                dny = 2.d0*dNormy
                xghost = x0 - dnx
                yghost = y0 - dny

                !irow = irow+1
                !a(irow,1) = xghost - x0
                !a(irow,2) = yghost - y0
                a(1,1) = a(1,1) + dnx*dnx
                a(1,2) = a(1,2) + dnx*dny
                a(2,2) = a(2,2) + dny*dny
! check that signs cancel, since my normal is inward pointing
                uvel = qMerge(2,i,j)
                vvel = qMerge(3,i,j)
                dot = uvel*alf + vvel*beta
                ughost = uvel - 2.d0*dot*alf
                vghost = vvel - 2.d0*dot*beta

                db(1) = 0.d0 ! pw const extrap of density
                db(2) = ughost - uvel
                db(3) = vghost - vvel
                db(4) = 0.d0 ! and pressure

                rhs(1,:) = rhs(1,:) + dx * db(:)
                rhs(2,:) = rhs(2,:) + dy * db(:)
              endif

             ntermsToUse = nterms
             if (nterms .eq. 5 .and. nborCount .lt.7) then
                write(*,222) i,j,nborCount
 222            format("Cell i,j = ",2i5," has only ",i6," neighbors",
     &                 /,"Drop to first order gradient")
                ntermsToUse = 2

             endif

             if (ntermsToUse .eq. 2) then
               c11 = sqrt(a(1,1))
               c12 = a(1,2) / c11
               c22 = sqrt(a(2,2) - c12**2)

               if (c22 .ne. 0.d0) then ! might have to test c11 too?
                  do m = 1, nvar
                    b(1) = rhs(1,m) / c11
                    b(2) = (rhs(2,m) - c12*b(1))/c22
                    gradmy(m,i,j) = b(2) / c22
                    gradmx(m,i,j) = (b(1) - c12*gradmy(m,i,j))/c11
                  end do
               else
                 write(*,*)"found c22=0 for cell ",i,j
                 ! set gradient to zero?
               endif

             ! solve a^t*a  * grad = a^t rhs.  First form cholesky factors.
             ! will have to add robustness checks
             else ! larger matrix, larger cholesky, 
               c11 = sqrt(a(1,1))
               c12 = a(1,2) / c11   
               c13 = a(1,3) / c11   
               c14 = a(1,4) / c11   
               c15 = a(1,5) / c11   

               c22 = sqrt(a(2,2)-c12**2)
               c23 = (a(2,3)-c12*c13)/c22
               c24 = (a(2,4)-c12*c14)/c22
               c25 = (a(2,5)-c12*c15)/c22

               c33 = sqrt(a(3,3)-c13**2-c23**2)
               c34 = (a(3,4)- c13*c14-c23*c24)/c33
               c35 = (a(3,5)- c13*c15-c23*c25)/c33

               c44 = sqrt(a(4,4)-c14**2-c24**2-c34**2)
               c45 = (a(4,5)-c14*c15-c24*c25-c34*c35)/c44
               c55 = sqrt(a(5,5)-c15**2-c25**2-c35**2-c45**2)
                  
             ! solving c^t b = rhs, then c b = gradients
             do m = 1, nvar
                b(1) = rhs(1,m)/c11
                b(2) = (rhs(2,m)-c12*b(1))/c22
                b(3) = (rhs(3,m)-c13*b(1)-c23*b(2))/c33
                b(4) = (rhs(4,m)-c14*b(1)-c24*b(2)-c34*b(3))/c44
                b(5) = (rhs(5,m)-c15*b(1)-c25*b(2)-c35*b(3) -
     &                            -c45*b(4))/c55

                w5 = b(5)/c55
                w4 = (b(4)-c45*w5)/c44
                w3 = (b(3)-c35*w5-c34*w4)/c33
                w2 = (b(2)-c25*w5-c24*w4-c23*w3)/c22
                w1 = (b(1)-c15*w5-c14*w4-c13*w3-c12*w2)/c11

                gradmy(m,i,j) = w2
                gradmx(m,i,j) = w1
             end do

             endif
c
c      apply limiter if requested. Go over all neighbors, do BJ
        if (nolimiter) go to 20

       if (limitTile .eq. 1) then
          call limitTileCellBJ(qmerge,gradmx,gradmy,xcentMerge,
     &                         ycentMerge,irr,nvar,mitot,mjtot,
     &                         nborList,nborCount,i,j) 

        else
           call limitTileGradientLP(qmerge,gradmx,gradmy,xcentMerge,
     &                           ycentMerge,xlow,ylow,dx,dy,irr,lwidth,
     &                           nvar,mitot,mjtot,lstgrd,areaMin,mptr,
     &                           lpChoice,nborList,nborCount)
       endif
        
 20     continue

c
      ! redo neighborhood calc as in makeNHoods but putting merged vals INTO cells
      ! instead of getting FROm cells. have to do it this way because a cell doesn't 
      ! know which cells contribute TO it, only the giving cell knows.

 35   continue

      valnew = 0.d0  !  all cells initialized to 0

c     next loop indices are assuming maxnco <= 2, so
c     dont bother looking at ghost cells further away
c     these ghost cells will be set next stage.
c     need to look at some ghost cells because they may distribute
c     to the last real cell
      !do 50 j = lwidth-1, mjtot-lwidth+2
      !do 50 i = lwidth-1, mitot-lwidth+2
      do 50 j = lwidth-1, mjtot-lwidth+1
      do 50 i = lwidth-1, mitot-lwidth+1
          k = irr(i,j)
          if (k .eq. -1) then
             ! set valnew to 'robust' fake state
             valnew(:,i,j) = fakeState
             go to 50  ! does not contribute
          endif
          if (k .eq. lstgrd .or. ar(k) .gt. areaMin) then
             valnew(:,i,j) = valnew(:,i,j)+qMerge(:,i,j)/numHoods(i,j)
             go to 50
          endif
          call getCellCentroid(lstgrd,i,j,xc,yc,
     &                         xlow,ylow,dx,dy,k)
          if (IS_OUTSIDE(xc,yc) .or. NOT_OK_GHOST(i,j)) then
             valnew(:,i,j) = qMerge(:,i,j)  ! copy what came in  
             go to 50 
          endif

             ioff = svi(k)
             joff = svj(k)
             ! cell i,j gives its tile to the offset cell as well as itself
             qm(:) = qMerge(:,i,j) + 
     .               (xc-xcentMerge(i,j))*gradmx(:,i,j) +
     .               (yc-ycentMerge(i,j))*gradmy(:,i,j)
             pr = .4d0*(qm(4)- 0.5d0*(qm(2)**2+qm(3)**2)/qm(1))
             if (qm(1) .le. 0.d0 .or. pr .le. 0.d0) then
                 qm(:) = qMerge(:,i,j)  ! is this conservative
             endif

             valnew(:,i,j) = valnew(:,i,j) + qm(:)/numHoods(i,j)

c            now contribute to neighboring cell
             koff = irr(i+ioff,j+joff)
             call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
     &                           xlow,ylow,dx,dy,koff)
             qm(:) = qMerge(:,i,j) + 
     .               (xcn-xcentMerge(i,j))*gradmx(:,i,j) +
     .               (ycn-ycentMerge(i,j))*gradmy(:,i,j)
             pr = .4d0*(qm(4)- 0.5d0*(qm(2)**2+qm(3)**2)/qm(1))
             if (qm(1) .le. 0.d0 .or. pr .le. 0.d0) then
                qm(:) = qMerge(:,i,j) ! is this conservative
             endif
             valnew(:,i+ioff,j+joff) = valnew(:,i+ioff,j+joff) +
     .                                 qm(:)/numHoods(i+ioff,j+joff)


!--             pr = .4d0*(qm(4)- 0.5d0*(qm(2)**2+qm(3)**2)/qm(1))
!--             ! only 2 cases - i/joff 0,0 or the at most one nbor
!--             if (ioff .eq. svi(k) .and. joff .eq. svj(k)) then
!--                 qm(:) = qMerge(:,i,j)
!--             else if (ioff .eq. 0 .and. joff .eq. 0) then
!--                 qm(:) = qMerge(:,i,j) + 
!--     .                    (xcn-xcentMerge(i,j))*gradmx(:,i,j) +
!--     .                    (ycn-ycentMerge(i,j))*gradmy(:,i,j)
!--                pr = .4d0*(qm(4)- 0.5d0*(qm(2)**2+qm(3)**2)/qm(1))
!--                if ((qm(1) .le. 0.d0) .or. (pr .le. 0.d0)) then
!--                   write(*,900)qm(1),pr,mptr,i,j,k,istage
!-- 900               format("shouldnt happen, rho,pr,mptr,i,j,k,istage",
!--     .                     2e15.7,5i4)
!--                   write(*,901) ioff,joff,numHoods(i+ioff,j+joff),
!--     .             ncount(i,j)
!-- 901               format("         reconstructing to offsets ",2i5,
!--     .                   " with nhoods=",i5," & ncount(i,j)",i3)
!--                   write(*,902) qmerge(:,i,j),q(:,i,j)
!-- 902               format("qMerge ",4e15.7,/,"q      ",4e15.7)
!--                endif
!--             endif
!--c            each cell takes its own merged val
!--c            shouldnt valnew be 0 for full cells?
!--             valnew(:,i+ioff,j+joff) = valnew(:,i+ioff,j+joff) +
!--     .                                 qm(:)/numHoods(i+ioff,j+joff)
!-- 40       continue

 50    continue
c
       q = valnew
c
c      check positivity
       call checkPhys(q,irr,mitot,mjtot,mptr,istage,lstgrd,
     .                   'from SRD_merge',lwidth,mitot-lwidth+1,
     .                   lwidth,mjtot-lwidth+1)

       if (verbose) then
          totmass2 =  bigconck(q,irr,mitot,mjtot,lwidth,nvar)
          dif = totmass2 - totmass
          write(*,912) totmass2,dif
 912      format("         mass after  redistribution is ",e25.15,
     .           "  dif is ",e15.7)
       endif

       return
       end
c
c ----------------------------------------------------------------------------
c
      subroutine makeNHood(volMerge,xcentMerge,ycentMerge,ncount,irr,
     .                     numHoods,mitot,mjtot,lwidth,lstgrd,xlow,ylow,
     .                     dx,dy,istage,mptr,
     .                     ffluxlen,gfluxlen,areaMin)

      use amr_module
      implicit double precision (a-h, o-z)

      dimension numHoods(mitot,mjtot), volMerge(mitot,mjtot)
      dimension xcentMerge(mitot,mjtot), ycentMerge(mitot,mjtot)
      dimension ffluxlen(mitot+1,mjtot+1),gfluxlen(mitot+1,mjtot+1)
      dimension ncount(mitot,mjtot), irr(mitot,mjtot)
      logical  IS_OUTSIDE, firstTimeThru
      logical debug/.false./
      logical works/.true./
      character ch


      IS_OUTSIDE(x,y) = (x .lt. xlower .or. x .gt. xupper .or.
     .                   y .lt. ylower .or. y .gt. yupper)

      ! merge until vqmerge at least this big (analogous to 1d where 
      ! left and right nhoods each dx

      numHoods = 1  ! initialize, everyone at least a member of its own hood
      ncount = 0    ! and its nborhood only includes itself to start
      maxnco = 0    ! will keep track of max

c     this code below counts on fact that don't need more than
c     2 cells to a side for a merging neighborhood
      !do 10 j = 3, mjtot-2
      !do 10 i = 3, mitot-2
      do 10 j = lwidth/2, mjtot-lwidth/2
      do 10 i = lwidth/2, mitot-lwidth/2
         k = irr(i,j)  
         if (k .eq. -1) go to 10
         call getCellCentroid(lstgrd,i,j,xc,yc,
     &                         xlow,ylow,dx,dy,k)
         if (IS_OUTSIDE(xc,yc)) go to 10  
         svi(k) = 0
         svj(k) = 0
         if (k .eq. lstgrd) then
            go to 10 ! a full  flow cell doesnt have to merge with a nbor 
         endif
         if (ar(k) .gt. areaMin) then
           ! nothing needs to be merged 
           go to 10
         endif
         vqmerge = ar(k)
         firstTimeThru = .true.
         nco = 1   ! initial size of neighborhood, from -1 to 1 square centered on cell
         ! next lines are for unstable corner that havent figured out
         ! yet in channel problem upper right corner
         !if (i.ge.47 .and. j.ge. 44) then
         !   areaMin = 2.d0*ar(lstgrd)
         !else
         !   areaMin = 0.5d0*ar(lstgrd)  
         !endif

         ! set ioff,joff to neighbor cell in most normal direction
         ! to cut cell boundary
         call getAdjCell(i,j,k,ioff,joff,ffluxlen,gfluxlen,
     &                   mitot,mjtot,dx,dy)
                  
         koff = irr(i+ioff,j+joff)
         call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
     &                        xlow,ylow,dx,dy,koff)
                  
         vqmerge = vqmerge + ar(koff)
         if (firstTimeThru) then ! count everybody
              numHoods(i+ioff,j+joff)=numHoods(i+ioff,j+joff)+1
              svi(k) = ioff   ! save indices of merging nhood cell
              svj(k) = joff   ! assuming only one here 
         endif
 
        if (vqmerge > areaMin) then
           ncount(i,j) = 1  ! this subroutine wont work if need > 1
           maxnco = max(maxnco,nco)
           go to 10
        else   ! redo with larger neighborhood
           write(*,909)i,j,ioff,joff
 909       format("cell ",2i4," not large enough w/ nhbor",2i4)
           works = .false.
        endif
            
 10   continue
      if (.not. works) then
         write(*,*)"Stopping calculation"
         stop
      endif

c

!   initialize array with  most common case, overwritten below as needed
      volMerge = ar(lstgrd) 
      xcentMerge = 0.d0      
      ycentMerge = 0.d0

      !do 20 j = 3, mjtot-2
      !do 20 i = 3, mitot-2
      do 20 j = lwidth/2, mjtot-lwidth/2
      do 20 i = lwidth/2, mitot-lwidth/2
         k = irr(i,j)  
         if (k .eq. lstgrd) then
             call getCellCentroid(lstgrd,i,j,xcentMerge(i,j),
     .                            ycentMerge(i,j),xlow,ylow,dx,dy,k)
             go to 20 ! flow or large cut cell is its own merged nhood
         endif
         if (k .eq. -1) then
            volMerge(i,j) = 0.d0
            go to 20
         endif
         call getCellCentroid(lstgrd,i,j,xcentMerge(i,j),
     .                        ycentMerge(i,j),xlow,ylow,dx,dy,k)
         if (IS_OUTSIDE(xcentMerge(i,j),ycentMerge(i,j))) then
            volMerge(i,j) = 0.d0
            go to 20
         endif

         ! initialize, note it is diff variable than above, weighted by numhoods
         call getCellCentroid(lstgrd,i,j,xc,yc,xlow,ylow,dx,dy,k)
         vmerge = ar(k)/numHoods(i,j)
         xcent = ar(k)*xc/numHoods(i,j)   ! initial centroid
         ycent = ar(k)*yc/numHoods(i,j)

         if (ar(k).lt. areaMin) then  ! add second cell vol and centroid info
            ioff = svi(k)  ! add one more centroid, assuming max of one more only
            joff = svj(k)
            koff = irr(i+ioff,j+joff)
         
            call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,xlow,ylow,
     .                           dx,dy,koff)
            vmerge = vmerge + ar(koff)/numHoods(i+ioff,j+joff)
            xcent = xcent + xcn*ar(koff)/numHoods(i+ioff,j+joff)
            ycent = ycent + ycn*ar(koff)/numHoods(i+ioff,j+joff)
         endif

         ! save info 
         volMerge(i,j) = vmerge
         xcentMerge(i,j) = xcent/vmerge
         ycentMerge(i,j) = ycent/vmerge
             
 20   continue

      if (debug) then
        write(*,*)"makeNhood "
        do j = 1, mjtot
        do i = 1, mitot
            if (irr(i,j) .eq. -1) then
               ch = "*"
            else if (irr(i,j) .eq. lstgrd) then
               ch = " "
            else
               ch = "+"
            endif
            write(*,888)ch,i,j,ncount(i,j),numHoods(i,j),
     &                xcentMerge(i,j),ycentMerge(i,j),
     &                volMerge(i,j)
 888        format(A1,4i4,3e15.7)
        end do
        end do
      endif

      return
      end
