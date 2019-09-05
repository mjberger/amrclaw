c
c ---------------------------------------------------------------------
c
       subroutine SRD_cellMerge(q,nvar,irr,mitot,mjtot,qx,qy,lstgrd,
     .                      dx,dy,lwidth,xlow,ylow,istage,
     .                      ncount,numHoods,mptr,ffluxlen,gfluxlen)

       use amr_module
       implicit double precision (a-h, o-z)

       dimension q(nvar,mitot,mjtot),  irr(mitot,mjtot)
       dimension qx(nvar,mitot,mjtot), qy(nvar,mitot,mjtot)
       dimension gradmx(nvar,mitot,mjtot), gradmy(nvar,mitot,mjtot)
       dimension ffluxlen(mitot+1,mjtot+1),gfluxlen(mitot+1,mjtot+1)

       dimension valnew(nvar,mitot,mjtot), volMerge(mitot,mjtot)
       dimension qMerge(nvar,mitot,mjtot), numHoods(mitot,mjtot)
       dimension nCount(mitot,mjtot)
       dimension xcentMerge(mitot,mjtot),ycentMerge(mitot,mjtot)
       dimension fakeState(nvar), qm(nvar), rhs(2,nvar)
       dimension dumax(nvar),dumin(nvar),phimin(nvar)
       dimension graddot(nvar),alpha(nvar),recon(nvar)
       dimension a(2,2),b(2)
       character ch

       logical IS_OUTSIDE, REG_NBORS,NOT_VALID_VAL,NOT_OK_GHOST
       logical quad, nolimiter,verbose
       common /order2/ ssw, quad, nolimiter
       integer omp_get_max_threads
       integer maxthreads/1/
       logical noDiagonals


c  xlow,ylow   refers to the corner of grid (w/o ghost cells)
c  xlower,ylower refers to the domain

c  OLD WAY
c  istage will determine how many ghost cells can be trusted:
c  a1ll in 1st stage, 2 less in 2nd stage

c NEW WAY - do one stage at a time then copy in ghost cells at intermediate
c stages for neighboring grids.  With 4 ghost cells, this means that
c interior cells can always trust 2 cells to the side after 1 stage
c never use cells from exterior to the domain for SRD

       IS_OUTSIDE(x,y) = (x .lt. xlower .or. x .gt. xupper .or.
     .                  y .lt. ylower .or. y .gt. yupper)
       REG_NBORS(i,j,lstgrd) = (irr(i+1,j).eq.lstgrd .and. 
     .                          irr(i-1,j).eq.lstgrd .and. 
     .                          irr(i,j+1).eq.lstgrd .and. 
     .                          irr(i,j-1).eq.lstgrd)
      NOT_VALID_VAL(i,j) = (i>mitot-2*(istage-1) .or. i<2*(istage-1)+1
     .                .or.  j>mjtot-2*(istage-1) .or. j<2*(istage-1)+1)
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
c      set verbosity for debugging, but only if serial and one grid,
c      otherwise meaningless
c      maxthreads initialized to 1 above in case no openmp
!$     maxthreads = omp_get_max_threads()
       if (maxthreads .gt. 1 .or. numgrids(1) .gt. 1) then
          verbose = .false.
       else
          verbose = .true.
       endif
       verbose = .false.
       !verbose = .true.
       noDiagonals = .true.
       !noDiagonals = .false.
       eps = 1d-4

c      some initializations
       ar(lstgrd) = dx*dy   ! area of regular grid cell 
       qMerge   = 0.d0
       dx2 = 2.d0*dx
       dy2 = 2.d0*dy

       ! put in 'normal' values to prevent errors e.g. in converting to prim vars
       ! these fake vals are in conserved variables
       fakeState(1) = 1.d0
       fakeState(2) = 0.d0
       fakeState(3) = 0.d0
       fakeState(4) = 2.5d0 !(p/gm1, corresponding to fakestate in method)

c     first make neighborhoods - need count for each cells, and width (nhood above)
c     nCount is size of neighborhood, numHoods is number of merged nhoods each cells is in
      call makeNHood(volMerge,xcentMerge,ycentMerge,ncount,irr,numHoods,
     .               mitot,mjtot,lwidth,lstgrd,xlow,ylow,dx,dy,istage,
     .               mptr,noDiagonals,ffluxlen,gfluxlen,areaMin)

       if (verbose) then
          totmass =  bigconck(q,irr,mitot,mjtot,lwidth,nvar)
          write(*,911) totmass
 911      format(/,"         mass before redistribution is ",e25.15)
       endif

c       form qMerge vals 
        do 10 j = 1, mjtot
        do 10 i = 1, mitot
            k = irr(i,j)
            if (k .eq. -1) go to 10 ! no solid cells
            call getCellCentroid(lstgrd,i,j,xc,yc,xlow,ylow,dx,dy,k)
            if (IS_OUTSIDE(xc,yc) .or. NOT_OK_GHOST(i,j)) then 
              qMerge(:,i,j) = rinfinity ! to make sure we dont use it
              go to 10 
            endif
            if (k.eq.lstgrd) then
              qMerge(:,i,j) = q(:,i,j)
              go to 10
            endif
            if (ar(k) .gt. areaMin) then
              qMerge(:,i,j) = q(:,i,j)
              go to 10
            endif
c
c           sum blocks to left & rt using only valid cells 
c
            ioff = svi(k)
            joff = svj(k)           
            koff = irr(i+ioff,j+joff)
            call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
     &                           xlow,ylow,dx,dy,koff)

            ! add the max of one other nbor cells and yourself
            ! to get merged val
            qMerge(:,i,j) = ar(k)*q(:,i,j)/numHoods(i,j) +
     .            ar(koff)*q(:,i+ioff,j+joff)/numHoods(i+ioff,j+joff)

            qMerge(:,i,j) = qMerge(:,i,j) / volMerge(i,j)

 10     continue

        ! gradient of merge neighborhoods, initialized to 0. set using neighboring merged tiels
        gradmx = 0.d0
        gradmy = 0.d0

        do 20 j = 3, mjtot-2
        do 20 i = 3, mitot-2
            k = irr(i,j)
            if (k .eq. -1) go to 20 ! solid cells have no gradient
            call getCellCentroid(lstgrd,i,j,xc,yc,
     &                           xlow,ylow,dx,dy,k)
            if (IS_OUTSIDE(xc,yc)) go to 20  ! exterior cells dont contribute

            ! if completely regular, use regular qmerge gradients
            ! even reg cells need gradients to put their val in cut cells
            if (k .eq. lstgrd .and. REG_NBORS(i,j,lstgrd)) then
               gradmx(:,i,j) = (qMerge(:,i+1,j)- qMerge(:,i-1,j))/dx2
               gradmy(:,i,j) = (qMerge(:,i,j+1)- qMerge(:,i,j-1))/dy2
               go to 20
            endif

            if (ncount(i,j) .eq. 0) then 
            ! these should be all regular cells that are left
               nco = 1 ! need to use some kind of nhood to make gradient
            else 
               nco = ncount(i,j)
            endif

            ! can reconstruct high order poly here, but for now do linear reconstruction
            ! solving least squares problem with all neighboring merged vals. Should preprocess the matrix and save
            rhs = 0.d0 ! initialize for accumulation
            a = 0.d0
            x0 = xcentMerge(i,j)
            y0 = ycentMerge(i,j)
            do 22 joff = -nco, nco
            do 22 ioff = -nco, nco
                if (ioff .eq. 0 .and. joff .eq. 0) go to 22 ! no eqn to solve
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 22
c               if (noDiagonals) then
c                 if (abs(ioff).eq.1 .and. abs(joff).eq.1) go to 22
c               endif
                call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
     &                               xlow,ylow,dx,dy,koff)
                if (IS_OUTSIDE(xcn,ycn)) go to 22
               deltax = xcentMerge(i+ioff,j+joff) - x0
               deltay = ycentMerge(i+ioff,j+joff) - y0
               a(1,1) = a(1,1) + deltax*deltax
               a(1,2) = a(1,2) + deltax*deltay
               a(2,2) = a(2,2) + deltay*deltay
               rhs(1,:) = rhs(1,:) + deltax * 
     .                    (qMerge(:,i+ioff,j+joff) - qMerge(:,i,j))
               rhs(2,:) = rhs(2,:) + deltay * 
     .                    (qMerge(:,i+ioff,j+joff) - qMerge(:,i,j))
 22          continue

             ! solve a^t*a  * grad = a^t rhs.  First form cholesky factors.
             ! will have to add robustness checks
             c11 = sqrt(a(1,1))
             c12 = a(1,2) / c11
             c22 = sqrt(a(2,2) - c12**2)

             ! now back solve (C^t C = rhs of A^tdu ) to get x and y gradient for all variables
             if (c22 .ne. 0.d0) then ! might have to test c11 too?
                do m = 1, nvar
                  b(1) = rhs(1,m) / c11
                  b(2) = (rhs(2,m) - c12*b(1))/c22
                  gradmy(m,i,j) = b(2) / c22
                  gradmx(m,i,j) = (b(1) - c12*gradmy(m,i,j))/c11
                end do
             endif
 20     continue
c
c      apply limiter if requested. Go over all neighbors, do BJ
        if (nolimiter) go to 35
        do 30 j = 3, mjtot-2
        do 30 i = 3, mitot-2
            k = irr(i,j)
            if (k .eq. -1) go to 30 ! solid cells have no gradient
            if (numHoods(i,j).eq.1 .and. ncount(i,j) .eq.0) go to 30  ! CHECK THAT NOTHING TO DO AND VAL NOT CHANGED
            nco = ncount(i,j)

            ! find max and min needed for BJ limiting
            dumax = 0.d0 
            dumin = 0.d0
            do 31 joff = -nco, nco
            do 31 ioff = -nco, nco
              if (ioff .eq. 0 .and. joff .eq. 0) go to 31
              koff = irr(i+ioff,j+joff)
              if (koff .eq. -1) go to 31
c             if (noDiagonals) then
c                if (abs(ioff).eq.1 .and. abs(joff).eq.1) go to 31
c             endif
              call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
     &                             xlow,ylow,dx,dy,koff)
              if (IS_OUTSIDE(xcn,ycn)) go to 31
              dumax = max(dumax,qmerge(:,i+ioff,j+joff)-qmerge(:,i,j))
              dumin = min(dumin,qmerge(:,i+ioff,j+joff)-qmerge(:,i,j))
 31         continue

            phimin = 1.d0
            do 32 joff = -nco, nco
            do 32 ioff = -nco, nco
                !! there may be an equation here since have to reconstruct
                !! to cell center, not merge center
                !! but presumably one would only check positivity since nothing
                !! to limit using same qMerge nhood
                if (ioff .eq. 0 .and. joff .eq. 0) go to 32 
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 32
c               if (noDiagonals) then
c                  if (abs(ioff).eq.1 .and. abs(joff).eq.1) go to 32
c               endif
                call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
     &                               xlow,ylow,dx,dy,koff)
                if (IS_OUTSIDE(xcn,ycn)) go to 32
                ! use merged val at ioff,joff to limit
                !diffx = xcentMerge(i+ioff,j+joff)-xcentMerge(i,j)
                !diffy = ycentMerge(i+ioff,j+joff)-ycentMerge(i,j)
                ! these next lines limit at cell where will evaluate
                ! previous lines limit at cells used in making gradient
                diffx = xcn-xcentMerge(i,j)
                diffy = ycn-ycentMerge(i,j)
                graddot  = gradmx(:,i,j)*diffx + gradmy(:,i,j)*diffy
                recon = qMerge(:,i,j) + graddot  
                  do m = 1,4
                     if (graddot(m) > 0.d0) then
                        alpha(m) = min(1.d0, dumax(m)/graddot(m))
                     else if (graddot(m) < 0.d0) then
                        alpha(m) = min(1.d0, dumin(m)/graddot(m))
                     else
                        alpha(m) = 1.d0
                     endif
                  end do
                  ! one last check for positivity
                  xymomsq = recon(2)**2+recon(3)**2
                  press = .4d0*(recon(4)-0.5d0*xymomsq/recon(1))
                  if (recon(1).le.eps .or. press.le.eps) alpha = 0.d0
                  phimin = min(phimin, alpha)
 32         continue
            gradmx(:,i,j) = gradmx(:,i,j)*phimin(:)
            gradmy(:,i,j) = gradmy(:,i,j)*phimin(:)
 30     continue

 
c      gradmx = 0.d0  
c      gradmy = 0.d0 
c
      ! redo neighborhood calc as in makeNHoods but putting merged vals INTO cells
      ! instead of getting FROm cells. have to do it this way because a cell doesn't 
      ! know which cells contribute TO it, only the giving cell knows.

 35   continue

      valnew = 0.d0  !  all cells initialized to 0

c     do 50 j = 1, mjtot
c     do 50 i = 1, mitot
c     next loop indices are assuming maxnco <= 2, so
c     dont bother looking at ghost cells further away
c     these ghost cells will be set next stage.
c     dont to prevent verbose error that happen in ghost cells
      do 50 j = lwidth-1, mjtot-lwidth+2
      do 50 i = lwidth-1, mitot-lwidth+2
          k = irr(i,j)
          if (k .eq. -1) then
             ! set valnew to 'robust' fake state
             valnew(:,i,j) = fakeState
             go to 50  ! does not contribute
          endif
          if (k .eq. lstgrd) then
             valnew(:,i,j) = valnew(:,i,j)+qMerge(:,i,j)/numHoods(i,j)
             go to 50
          endif
          call getCellCentroid(lstgrd,i,j,xc,yc,
     &                         xlow,ylow,dx,dy,k)
          if (IS_OUTSIDE(xc,yc) .or. NOT_OK_GHOST(i,j)) then
             valnew(:,i,j) = qMerge(:,i,j)  ! copy what came in  
             go to 50 
          endif

          do 40 joff = -ncount(i,j), ncount(i,j)
          do 40 ioff = -ncount(i,j), ncount(i,j)
             koff = irr(i+ioff,j+joff)
             if (koff .eq. -1) go to 40
             if (abs(ioff).eq.1 .and. abs(joff).eq.1) go to 40
             call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
     &                            xlow,ylow,dx,dy,koff)
             if (IS_OUTSIDE(xcn,ycn)) go to 40  

             ! only 2 cases - i/joff 0,0 or the at most one nbor
             if ((ioff .eq. svi(k) .and. joff .eq. svj(k)) .or.
     .           (ioff .eq. 0      .and. joff .eq. 0))     then
c            form qm, check if state is physical 
                qm(:) = qMerge(:,i,j) + 
     .                    (xcn-xcentMerge(i,j))*gradmx(:,i,j) +
     .                    (ycn-ycentMerge(i,j))*gradmy(:,i,j)
                pr = .4d0*(qm(4)- 0.5d0*(qm(2)**2+qm(3)**2)/qm(1))
                if ((qm(1) .le. 0.d0) .or. (pr .le. 0.d0)) then
                   write(*,900)qm(1),pr,mptr,i,j,istage
 900               format("should not happen, rho,pr,mptr,i,j,istage",
     .                     2e15.7,4i4)
                   write(*,901) ioff,joff,numHoods(i+ioff,j+joff)
 901               format("         reconstructing to offsets ",2i5,
     .                   " with ",i5," nhoods")
                endif
c               each cell takes its own merged val
c               shouldnt valnew be 0 for full cells?
                valnew(:,i+ioff,j+joff) = valnew(:,i+ioff,j+joff) +
     .                                    qm(:)/numHoods(i+ioff,j+joff)
             endif
 40       continue

 50    continue
c
       q = valnew
c
c      check positivity
       call checkPhysInt(q,mitot,mjtot,mptr,istage,
     .                   lwidth,'from SRD_merge')

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
      double precision function bigconck(q,irr,mitot,mjtot,lwidth,nvar)
c
      use amr_module
      implicit double precision (a-h,o-z)

      dimension q(nvar,mitot,mjtot), irr(mitot,mjtot)

      totmass = 0.d0

      do 10 j = lwidth+1, mjtot-lwidth
      do 10 i = lwidth+1, mitot-lwidth

         k = irr(i,j)
         if (k .eq. -1) cycle 
         totmass = totmass + q(1,i,j)*ar(k)

 10   continue

      bigconck = totmass

      return
      end

c
c ----------------------------------------------------------------------------
c
      subroutine makeNHood(volMerge,xcentMerge,ycentMerge,ncount,irr,
     .                     numHoods,mitot,mjtot,lwidth,lstgrd,xlow,ylow,
     .                     dx,dy,istage,mptr,noDiagonals,
     .                     ffluxlen,gfluxlen,areaMin)

      use amr_module
      implicit double precision (a-h, o-z)

      dimension numHoods(mitot,mjtot), volMerge(mitot,mjtot)
      dimension xcentMerge(mitot,mjtot), ycentMerge(mitot,mjtot)
      dimension ffluxlen(mitot+1,mjtot+1),gfluxlen(mitot+1,mjtot+1)
      dimension ncount(mitot,mjtot), irr(mitot,mjtot)
      logical noDiagonals
      logical NOT_OK_GHOST, IS_OUTSIDE, firstTimeThru
      logical debug/.false./, IS_OUT_OF_RANGE
      logical works/.true./
      character ch

      NOT_OK_GHOST(i,j) = (i .lt. 3 .or. 
     .                     i .gt. mitot-2 .or.
     .                     j .lt. 3 .or. 
     .                     j .gt. mjtot-2)
      IS_OUTSIDE(x,y) = (x .lt. xlower .or. x .gt. xupper .or.
     .                   y .lt. ylower .or. y .gt. yupper)
      IS_OUT_OF_RANGE(i,j)  = (i<1 .or. i>mitot .or. j<1 .or. j>mjtot)

      ! merge until vqmerge at least this big (analogous to 1d where 
      ! left and right nhoods each dx

      numHoods = 1  ! initialize, everyone at least a member of its own hood
      ncount = 0    ! and its nborhood only includes itself to start
      maxnco = 0    ! will keep track of max
      svi = 0       
      svj = 0

c     this code below counts on fact that don't need more than
c     2 cells to a side for a merging neighborhood
      do 10 j = 3, mjtot-2
      do 10 i = 3, mitot-2
         k = irr(i,j)  
         if (k .eq. -1) go to 10
         call getCellCentroid(lstgrd,i,j,xc,yc,
     &                         xlow,ylow,dx,dy,k)
         if (IS_OUTSIDE(xc,yc)) go to 10  
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

      do 20 j = 3, mjtot-2
      do 20 i = 3, mitot-2
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
c
c ------------------------------------------------------------------------
c
      subroutine checkPhysInt(q,mitot,mjtot,mptr,istage,nghost,str)

      implicit real*8 (a-h,o-z)
      dimension q(4,mitot,mjtot)
      character*14 str

c     this only checks the interior of the grid and
c     excluse ghost cells
c
      gamma1 = .4d0
      do j = nghost+1, mjtot-nghost
      do i = nghost+1, mitot-nghost
         rho = q(1,i,j)
         u = q(2,i,j)/rho
         v = q(3,i,j)/rho
         pr = gamma1*(q(4,i,j)-.5d0*rho*(u*u+v*v))
         if (rho <= 0 .or. pr <= 0) then
            write(*,901)rho,pr,i,j,mptr,istage,str
 901        format("non-physical den/pr",2e15.7," at i,j grid stage ",
     .           4i5,2x,a14)
         endif
      end do
      end do

      return
      end
 

c
c ------------------------------------------------------------------------
c 
      subroutine getAdjCell(i,j,k,ioff,joff,ffluxlen,gfluxlen,
     .                      mitot,mjtot,dx,dy)

      use amr_module
      implicit double precision (a-h, o-z)
      dimension ffluxlen(mitot+1,mjtot+1),gfluxlen(mitot+1,mjtot+1)
      dimension rl(4),iset(4),jset(4)
      data iset/-1, 1, 0, 0/
      data jset/0, 0, -1, 1/


c  for cell i,j return offsets for most normal neighbor
c  either look at normal direction, or full edge lengths
      
      rl(1) = ffluxlen(i,j)/dy
      rl(2) = ffluxlen(i+1,j)/dy
      rl(3) = gfluxlen(i,j)/dx
      rl(4) = gfluxlen(i,j+1)/dx

c take max over all and thats the neighbor
      ispot = 1
      do index = 2, 4
         if (rl(index) .gt. rl(ispot)) ispot = index  ! always positive
      end do
      ioff = iset(ispot)
      joff = jset(ispot)

!!    will need to chceck that i+ioff,j+joff is inside
!!    domain and valid, but for cylinder case for now
c     skip it

      return
      end
      
