c
c ---------------------------------------------------------------------
c
       subroutine SRD_cellMerge(q,nvar,irr,mitot,mjtot,lstgrd,
     .                      dx,dy,lwidth,xlow,ylow,istage,
     .                      ncount,numHoods,mptr)

       use amr_module
       implicit double precision (a-h, o-z)

       include "cuserdt.i"
       dimension q(nvar,mitot,mjtot),  irr(mitot,mjtot)
       dimension gradmx(nvar,mitot,mjtot), gradmy(nvar,mitot,mjtot)

       dimension valnew(nvar,mitot,mjtot)
       dimension qMerge(nvar,mitot,mjtot), numHoods(mitot,mjtot)
       dimension nCount(mitot,mjtot)
       dimension mioff(mitot,mjtot),mjoff(mitot,mjtot)
       dimension fakeState(nvar), qm(nvar), rhs(5,nvar)
       dimension a(5,5),b(5),db(nvar)
       dimension nborList(35,2)
       character c2

       logical IS_OUTSIDE, NOT_OK_GHOST
       logical OUT_OF_RANGE
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

      OUT_OF_RANGE(i,j) = (i .lt. 1 .or. i .gt. mitot .or.
     &                     j .lt. 1 .or. j .gt. mjtot)

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
      call makeNHood(ncount,irr,numHoods,
     .               mitot,mjtot,lwidth,lstgrd,xlow,ylow,dx,dy,istage,
     .               mptr,areaMin)

      call merge_shifts(irr,mitot,mjtot,lwidth,dx,dy,xlow,ylow,lstgrd,
     .                  numHoods,ncount)


c       form qMerge vals 
        do 10 j = lwidth/2, mjtot-lwidth/2
        do 10 i = lwidth/2, mitot-lwidth/2
            k = irr(i,j)
            if (k .eq. -1) go to 10 ! no solid cells
            call getCellCentroid(lstgrd,i,j,xc,yc,xlow,ylow,dx,dy,k)
            if (IS_OUTSIDE(xc,yc) .or. OUT_OF_RANGE(i,j)) then 
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
c           call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
c    &                           xlow,ylow,dx,dy,koff)

            ! add the one other nbor cell and yourself to
            ! get merged val. volmerge should have same weighting
            qMerge(:,i,j) = (ar(k)*q(:,i,j)/numHoods(i,j) +
     .            ar(koff)*q(:,i+ioff,j+joff)/numHoods(i+ioff,j+joff))/
     .            volMerge(k)
 10     continue

        ! gradient of merging neighborhoods, initialized to 0. set using neighboring tiles
        if (ssw .eq. 0.d0) then
           gradmx = 0.d0  ! so can use same code below with gradients
           gradmy = 0.d0
           go to 35   ! ssw = 0 means no slopes
        else
           gradmx = rinfinity   ! 0.d0 for debugging to make sure set
           gradmy = rinfinity   ! 0.d0
        endif

        ! compute stable neighborhood for gradient on merging tiles
        call makeMergeNHood(irr,lwidth,
     &                      mitot,mjtot,lstgrd,dx,dy,xlow,ylow,mptr,
     &                      mioff,mjoff)

        call qmslopes(irr,mitot,mjtot,lwidth,dx,dy,xlow,ylow,lstgrd,
     &                numHoods,ncount,areaMin,
     &                mioff,mjoff,qMerge,nvar,gradmx,gradmy)
c
c      apply limiter if requested. Go over all neighbors, do BJ
        if (nolimiter) go to 20

       if (limitTile .eq. 1) then
c         CELL based version
c         call limitTileBJ2(qmerge,gradmx,gradmy,
c    &                         xlow,ylow,dx,dy,
c    &                         irr,nvar,mitot,mjtot,lstgrd,
c    &                         i,j,k,mioff,mjoff,lwidth) 
c         limit all of grid all at once
          call limitTileGradientBJ(qmerge,gradmx,gradmy,
     &                         xlow,ylow,dx,dy,
     &                         mioff,mjoff, 
     &                         irr,nvar,mitot,mjtot,lstgrd,lwidth)

        else
           call limitTileGradientLP(qmerge,gradmx,gradmy,
     &                           xlow,ylow,dx,dy,irr,lwidth,
     &                           nvar,mitot,mjtot,lstgrd,areaMin,mptr,
     &                           lpChoice,mioff,mjoff,
     &                           nborList,nborCount)
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
      do 50 j = lwidth-1, mjtot-lwidth+2
      do 50 i = lwidth-1, mitot-lwidth+2
      !do 50 j = lwidth, mjtot-lwidth
      !do 50 i = lwidth, mitot-lwidth
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
          if (IS_OUTSIDE(xc,yc) .or. OUT_OF_RANGE(i,j)) then
             valnew(:,i,j) = qMerge(:,i,j)  ! copy what came in  
             go to 50 
          endif

             ioff = svi(k)
             joff = svj(k)
             ! cell i,j gives its tile to the offset cell as well as itself
             qm(:) = qMerge(:,i,j) + 
     .               (xc-xcentMerge(k))*gradmx(:,i,j) +
     .               (yc-ycentMerge(k))*gradmy(:,i,j)
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
     .               (xcn-xcentMerge(k))*gradmx(:,i,j) +
     .               (ycn-ycentMerge(k))*gradmy(:,i,j)
c            pr = .4d0*(qm(4)- 0.5d0*(qm(2)**2+qm(3)**2)/qm(1))
c            if (qm(1) .le. 0.d0 .or. pr .le. 0.d0) then
c               qm(:) = qMerge(:,i,j) ! is this conservative
c            endif
             valnew(:,i+ioff,j+joff) = valnew(:,i+ioff,j+joff) +
     .                                 qm(:)/numHoods(i+ioff,j+joff)

 50    continue
c
       q = valnew
c
c      check positivity
       call checkPhys(q,irr,mitot,mjtot,mptr,istage,lstgrd,
     .                   'from SRD_merge',lwidth+1,mitot-lwidth,
     .                   lwidth+1,mjtot-lwidth)

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
      subroutine makeNHood(ncount,irr,
     .                     numHoods,mitot,mjtot,lwidth,lstgrd,xlow,ylow,
     .                     dx,dy,istage,mptr,
     .                     areaMin)

      use amr_module
      implicit double precision (a-h, o-z)

      dimension numHoods(mitot,mjtot)
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
      do 10 j = lwidth/2, mjtot-lwidth/2
      do 10 i = lwidth/2, mitot-lwidth/2
      !do 10 j = lwidth+1, mjtot-lwidth
      !do 10 i = lwidth+1, mitot-lwidth
         k = irr(i,j)  
         if (k .eq. -1) go to 10
         call getCellCentroid(lstgrd,i,j,xc,yc,
     &                         xlow,ylow,dx,dy,k)
         if (IS_OUTSIDE(xc,yc)) then
            if (k .ne. lstgrd) then ! set centroid for later testing but shouldnt be using
              xcentMerge(k) = xc
              ycentMerge(k) = yc
              go to 10  
            endif
         endif
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
         call getAdjCell(i,j,k,ioff,joff,
     &                   mitot,mjtot,dx,dy)
                  
         koff = irr(i+ioff,j+joff)
         call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
     &                        xlow,ylow,dx,dy,koff)
                  
         vqmerge = vqmerge + ar(koff)
         if (firstTimeThru) then ! count everybody
              numHoods(i+ioff,j+joff)=numHoods(i+ioff,j+joff)+1
              svi(k) = ioff   ! save indices of merging nhood cell
              svj(k) = joff   ! assuming only one here 
              if (IS_OUTSIDE(xcn,ycn)) then
                 write(*,*)"do need to fix this"
                 stop
              endif
         endif
 
        if (vqmerge >= areaMin) then
           ncount(i,j) = 1  ! this subroutine wont work if need > 1
           maxnco = max(maxnco,nco)
           go to 10
        else   ! redo with larger neighborhood
           write(*,909)i,j,ioff,joff,vqmerge,areaMin
 909       format("cell ",2i4," not large enough w/ nhbor",2i4,
     .            " vqmerge ",e15.7," areamin r", e15.7)
           works = .false.
        endif
            
 10   continue
      if (.not. works) then
         write(*,*)"Stopping calculation"
         stop
      endif

c

!   initialize array with  most common case, overwritten below as needed

      do 20 j = lwidth/2, mjtot-lwidth/2
      do 20 i = lwidth/2, mitot-lwidth/2
      !do 20 j = lwidth+1, mjtot-lwidth
      !do 20 i = lwidth+1, mitot-lwidth
         k = irr(i,j)  
         if (k .eq. lstgrd) go to 20 ! flow cell is its own merged nhood
         if (k .eq. -1) go to 20 ! solid cell is its own merged nhood
             ! flow or large cut cell is its own merged nhood
         call getCellCentroid(lstgrd,i,j,xc,
     .                        yc,xlow,ylow,dx,dy,k)
         if (IS_OUTSIDE(xc,yc)) go to 20

         ! initialize, note it is diff variable than above, weighted by numhoods
         ! small cut cell - compute merging vol and centroids
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
         volMerge(k) = vmerge
         xcentMerge(k) = xcent/vmerge
         ycentMerge(k) = ycent/vmerge
             
 20   continue

      if (debug) then
        volMerge(lstgrd) = dx*dy
        xcentMerge(lstgrd) = 0.d0
        ycentMerge(lstgrd) = 0.d0
        write(*,*)"makeNhood "
        do j = 1, mjtot
        do i = 1, mitot
            k = irr(i,j)
            if (k .eq. -1) then
               ch = "*"
               write(*,888)ch,i,j
            else if (k .eq. lstgrd) then
               ch = " "
               write(*,888)ch,i,j
            else
               ch = "+"
               write(*,888)ch,i,j,ncount(i,j),numHoods(i,j),
     &                   xcentMerge(k),ycentMerge(k),volMerge(k)
 888        format(A1,4i4,3e15.7)
            endif
        end do
        end do
      endif

      return
      end
c
c ----------------------------------------------------------------------
c
        subroutine  makeMergeNHood(irr,lwidth,mitot,mjtot,lstgrd,dx,dy,
     &                             xlow,ylow,mptr,mioff,mjoff)

        use amr_module
        implicit double precision (a-h,o-z)
        include "cuserdt.i"
        dimension irr(mitot,mjtot)
        dimension mioff(mitot,mjtot),mjoff(mitot,mjtot)
        logical NOT_OK_GHOST, OKstencil, IS_OUTSIDE, OUT_OF_RANGE
        logical debug/.false./

        IS_OUTSIDE(x,y) = (x .lt. xlower .or. x .gt. xupper .or.
     .                     y .lt. ylower .or. y .gt. yupper)

        NOT_OK_GHOST(i,j) = (i .lt. 3 .or. 
     .                       i .gt. mitot-2 .or.
     .                       j .lt. 3 .or. 
     .                       j .gt. mjtot-2)


        OUT_OF_RANGE(i,j) =  (i .lt. 1 .or. i .gt. mitot .or.
     .                        j .lt. 1 .or. j .gt. mjtot)


        ! to get same as previous behavior uncomment next line
        ! if (numMergeTerms .eq. 2) return  

        if (numMergeTerms .eq. 2) then
          initVal = 1   ! initial neighborhood is 1 to each side
          tolx = 0.5d0 * dx
          toly = 0.5d0 * dy
        else ! numMergeTerms is 5, want 2nd gradient
          tolx = 1.5d0*dx 
          toly = 1.5d0*dy 
          initVal = 2   ! to fit quadratic need more cells 
        endif
        mioff = initVal
        mjoff = initVal

        !do j = lwidth+1, mjtot-lwidth
        !do i = lwidth+1, mitot-lwidth
        do j = lwidth-1, mjtot-lwidth+2
        do i = lwidth-1, mitot-lwidth+2
           k = irr(i,j)
           if (k .eq. -1 .or. k .eq. lstgrd) cycle !default mi/mjoff works

           OKstencil = .false.
           do while (.not. OKstencil) 
             diffx = -1.d0
             diffy = -1.d0
             ! form nhood, see if stable or need to expand
             OKstencil = .true.
             do joff = -mjoff(i,j), mjoff(i,j)
             do ioff = -mioff(i,j), mioff(i,j)
                if (OUT_OF_RANGE(i+ioff,j+joff)) cycle
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) cycle
                if (koff .eq. lstgrd) then
                   xm = xlow + (i+ioff-0.5d0)*dx
                   ym = ylow + (j+joff-0.5d0)*dy
                else
                   xm = xcentMerge(koff)
                   ym = ycentMerge(koff)
                endif
                if (IS_OUTSIDE(xm,ym) .or.
     &              NOT_OK_GHOST(i+ioff,j+joff)) cycle
                diffx = max(diffx, dabs(xm-xcentMerge(k)))
                diffy = max(diffy, dabs(ym-ycentMerge(k)))
             end do
             end do

             if (diffx .lt. tolx) then
                mioff(i,j) = mioff(i,j) + 1
                if (debug) write(*,900) i, j
 900            format("increasing mioff for cell ",2i5)
                OKstencil = .false.
             endif
             if (diffy .lt. toly) then
                mjoff(i,j) = mjoff(i,j) + 1
                if (debug) write(*,901) i, j
 901            format("increasing mjoff for cell ",2i5)
                OKstencil = .false.
             endif
             end do
        end do
        end do

        return
        end
