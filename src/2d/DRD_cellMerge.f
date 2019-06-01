c
c ---------------------------------------------------------------------
c
       subroutine DRD_cellMerge(q,nvar,irr,mitot,mjtot,qx,qy,lstgrd,
     .                      dx,dy,lwidth,xlow,ylow,istage,
     .                      ncount,numHoods)

       use amr_module
       implicit double precision (a-h, o-z)

       dimension q(mitot,mjtot,nvar),  irr(mitot,mjtot)
       dimension qx(mitot,mjtot,nvar), qy(mitot,mjtot,nvar)
       dimension gradmx(mitot,mjtot,nvar), gradmy(mitot,mjtot,nvar)

       dimension delta(mitot,mjtot,nvar), volMerge(mitot,mjtot)
       dimension denvolMerge(mitot,mjtot)
       dimension qMerge(mitot,mjtot,nvar), numHoods(mitot,mjtot)
       dimension nCount(mitot,mjtot)
       dimension xcentMerge(mitot,mjtot),ycentMerge(mitot,mjtot)
       dimension fakeState(nvar), qm(nvar), rhs(2,nvar)
       dimension dumax(nvar),dumin(nvar),phimin(nvar)
       dimension graddot(nvar),alpha(nvar),recon(nvar)
       dimension a(2,2),b(2)
       real*8 minmod

       logical IS_GHOST, IS_FAR_GHOST, verbose/.true./
       logical quad, nolimiter
       common /order2/ ssw, quad, nolimiter

c  this next statement considers a ghost cell to be anything beyond the
c  loop indices, since they change according to the stage
       IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                  j .le. lwidth .or. j .gt. mjtot-lwidth)

c :::::::::::::;:::
c
c   try cell merging and reconstruction to stabilize updates in small cut cells
c   calculate all updates using provisional values, then actually do the update
c   (as in, jacobi rather than gauss seidel)
c ::::::::::::::::

c
       write(*,*)"Havent fixed indexing yet"
       stop

       ar(-1) = 0.d0        ! zero out area of solid cells for this loop
       ar(lstgrd) = dx*dy   ! area of regular grid cell 
       qMerge   = 0.d0

       ! put in 'normal' values to prevent errors e.g. in converting to prim vars
       fakeState(1) = 1.d0
       fakeState(2) = 0.d0
       fakeState(3) = 0.d0
       fakeState(4) = 2.5d0

c     first make neighborhoods - need count for each cells, and width (nhood above)
c     nCount is size of neighborhood, numHoods is number of merged nhoods each cells is in
      call make_drdHood(volMerge,xcentMerge,ycentMerge,ncount,irr,
     .                  numHoods,mitot,mjtot,lwidth,lstgrd,xlow,ylow,
     .                  dx,dy)   

       if (verbose) then
          totmass =  bigconck(q,irr,mitot,mjtot,lwidth,nvar)
          write(*,911) totmass
 911      format(/,"         mass before redistribution is ",e30.20)
       endif

c       form qMerge vals for this step 
c       also form denvolMerge since needs density
        denvolMerge = 0.d0
        do 10 j = 1, mjtot
        do 10 i = 1, mitot
            if (irr(i,j) .eq. -1) go to 10 ! no solid cells
            if (irr(i,j) .eq. lstgrd .or. IS_GHOST(i,j)) then 
              qMerge(i,j,:) = q(i,j,:)
              go to 10 
            endif
c
c           sum in blocks of size 2*ncount on a side using only valid cells 
c
             do 27 joff = -ncount(i,j), ncount(i,j)
             do 27 ioff = -ncount(i,j), ncount(i,j)
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 27  ! solid cells dont contribute
                if (IS_GHOST(i+ioff,j+joff)) go to 27  ! nor ghost cells
c               count this cell 
                aoff = ar(koff)
                qMerge(i,j,:) = qMerge(i,j,:) + aoff*q(i+ioff,j+joff,:)
                denvolMerge(i,j) = denvolMerge(i,j) + 
     .                             aoff*q(i+ioff,j+joff,1)
 27          continue
             qMerge(i,j,:) = qMerge(i,j,:) / volMerge(i,j)
 10     continue

        ! gradient of merged neighborhoods, initialized to 0. 
        ! set using neighboring merged tiles
        gradmx = 0.d0
        gradmy = 0.d0

!!       if (ssw .eq. 0.d0) go to 35  ! nogradients needed

        do 20 j = lwidth+1, mjtot-lwidth
        do 20 i = lwidth+1, mitot-lwidth
            k = irr(i,j)
            if (k .eq. -1) go to 20 ! solid cells have no gradient
            if (ncount(i,j) .eq. 0) then 
            ! these should be all regular cells that are left
               nco = 1 ! need to use some kind of nhood to make gradient
            else 
               nco = ncount(i,j)
            endif

            ! can reconstruct high order poly here, but for now do linear reconstruction
            ! solving least squares problem with all neighboring merged vals. 
            ! Should preprocess the matrix and save
            rhs = 0.d0 ! initialize for accumulation
            a = 0.d0
            x0 = xcentMerge(i,j)
            y0 = ycentMerge(i,j)
            do 22 joff = -nco, nco
            do 22 ioff = -nco, nco
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 22
                if (IS_GHOST(i+ioff,j+joff)) go to 22
                if (ioff .eq. 0 .and. joff .eq. 0) go to 22 ! no eqn to solve
               deltax = xcentMerge(i+ioff,j+joff) - x0
               deltay = ycentMerge(i+ioff,j+joff) - y0
               a(1,1) = a(1,1) + deltax*deltax
               a(1,2) = a(1,2) + deltax*deltay
               a(2,2) = a(2,2) + deltay*deltay
               rhs(1,:) = rhs(1,:) + deltax * 
     .                    (qMerge(i+ioff,j+joff,:) - qMerge(i,j,:))
               rhs(2,:) = rhs(2,:) + deltay * 
     .                    (qMerge(i+ioff,j+joff,:) - qMerge(i,j,:))
 22          continue

             ! solve a^t*a  * grad = a^t rhs.  First form cholesky factors.
             ! will have to add robustness checks
             c11 = sqrt(a(1,1))
             c12 = a(1,2) / c11
             c22 = sqrt(a(2,2) - c12**2)

             ! now back solve (C^t C = rhs of A^tdu ) to get x and y gradient for all variables
             do m = 1, nvar
               b(1) = rhs(1,m) / c11
               b(2) = (rhs(2,m) - c12*b(1))/c22
               gradmy(i,j,m) = b(2) / c22
               gradmx(i,j,m) = (b(1) - c12*gradmy(i,j,m))/c11
             end do
 20     continue

c
c      apply limiter if requested. Go over all neighbors, do BJ
        if (nolimiter) go to 35
        do 30 j = lwidth+1, mjtot-lwidth
        do 30 i = lwidth+1, mitot-lwidth
            k = irr(i,j)
            if (k .eq. -1) go to 30 ! solid cells have no gradient
            if (numHoods(i,j) .eq. 1) go to 30  ! CHECK THAT NOTHING TO DO AND VAL NOT CHANGED
            nco = ncount(i,j)

            ! find max and min needed for BJ limiting
            dumax = 0.d0 
            dumin = 0.d0
            do 31 joff = -nco, nco
            do 31 ioff = -nco, nco
              if (ioff .eq. 0 .and. joff .eq. 0) go to 31
              if (IS_GHOST(i+ioff,j+joff)) go to 31
              koff = irr(i+ioff,j+joff)
              if (koff .eq. -1) go to 31
              dumax = max(dumax,qmerge(i+ioff,j+joff,:)-qmerge(i,j,:))
              dumin = min(dumin,qmerge(i+ioff,j+joff,:)-qmerge(i,j,:))
 31         continue

            phimin = 1.d0
            do 32 joff = -nco, nco
            do 32 ioff = -nco, nco
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 32
                if (IS_GHOST(i+ioff,j+joff)) go to 32
                if (ioff .eq. 0 .and. joff .eq. 0) go to 32 
                ! use merged val at ioff,joff to limit
                !diffx = xcentMerge(i+ioff,j+joff)-xcentMerge(i,j)
                !diffy = ycentMerge(i+ioff,j+joff)-ycentMerge(i,j)
                ! these next lines limit at cell where will evaluate
                ! previous lines limit at cells used in making gradient
                call getCellCentroid(lstgrd,i+ioff,j+joff,xc,yc,xlow,
     .                          ylow,dx,dy,koff)
                diffx = xc-xcentMerge(i,j)
                diffy = yc-ycentMerge(i,j)
                graddot  = gradmx(i,j,:)*diffx + gradmy(i,j,:)*diffy
                recon = qMerge(i,j,:) + graddot  
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
                  if (recon(1) .le. 0.d0) alpha = 0.d0
                  velsq = recon(2)**2+recon(3)**2
                  press = .4d0*(recon(4)-0.5d0*velsq/recon(1))
                  if (press .le. 0.d0) alpha = 0.d0
                  phimin = min(phimin, alpha)
 32         continue
            gradmx(i,j,:) = gradmx(i,j,:)*phimin(:)
            gradmy(i,j,:) = gradmy(i,j,:)*phimin(:)
 30     continue

 
c      gradmx = 0.d0  
c      gradmy = 0.d0 
c
      ! redo neighborhood calc as in makeNHoods but putting merged vals INTO cells
      ! instead of getting FROM cells. have to do it this way because a cell doesn't 
      ! know which cells contribute TO it, only the giving cell knows.

 35   continue

      do 50 j = 1, mjtot
      do 50 i = 1, mitot
          k = irr(i,j)
          if (k .eq. -1 .or. IS_GHOST(i,j)) then
             delta(i,j,:) = 0.d0
             go to 50  ! does not contribute
          endif
          call getCellCentroid(lstgrd,i,j,xc,yc,xlow,
     .                         ylow,dx,dy,k)
          qm(:) = qMerge(i,j,:) + 
     .                   (xc-xcentMerge(i,j))*gradmx(i,j,:)  +
     .                   (yc-ycentMerge(i,j))*gradmy(i,j,:)

          pr = .4d0 * (qm(4)- 0.5d0*(qm(2)**2+qm(3)**2)/qm(1))
          if ((qm(1) .le. 0.d0) .or. (pr .le. 0.d0)) then
             write(*,*)" should not happen"
          endif
         !! compute conservative diff between new and old vals to redistribute
          delta(i,j,:) = ar(k) * (qm(:) - q(i,j,:))
 50    continue
c
       call distribute(q,delta,volMerge,numHoods,ncount,
     .                 irr,mitot,mjtot,nvar,lstgrd,
     .                 denvolMerge)
c
c      q comes back updated
c
       if (verbose) then
          totmass2 =  bigconck(q,irr,mitot,mjtot,lwidth,nvar)
          dif = totmass2 - totmass

          write(*,912) totmass2,dif
 912      format("         mass after  redistribution is ",e30.20,
     .           "  dif is ",e15.7)
       endif

       return
       end
c
c ----------------------------------------------------------------------------
c
      subroutine make_drdHood(volMerge,xcentMerge,ycentMerge,ncount,irr,
     .                     numHoods,mitot,mjtot,lwidth,lstgrd,xlow,ylow,
     ,                     dx,dy)

       use amr_module
      implicit double precision (a-h, o-z)

      dimension numHoods(mitot,mjtot), volMerge(mitot,mjtot)
      dimension xcentMerge(mitot,mjtot), ycentMerge(mitot,mjtot)
      dimension ncount(mitot,mjtot), irr(mitot,mjtot)

      logical IS_GHOST, firstTimeThru
      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                 j .le. lwidth .or. j .gt. mjtot-lwidth)

      ! merge until vqmerge at least this big (analogous to 1d where left and right nhoods each dx
      !!areaMin = 2.d0*ar(lstgrd)  
      !!areaMin = 0.5d0*ar(lstgrd)  
      areaMin = ar(lstgrd)  
      numHoods = 0  ! initialize, loop below will add each cell to its own nhood
      ncount = 0
      ar(-1) = 0.d0  ! reset here to remind us
      volMerge = 0.d0
      eps = 1.d-12
      eps = -1.d-12

      do 10 j = lwidth+1, mjtot-lwidth
      do 10 i = lwidth+1, mitot-lwidth

         k = irr(i,j)  
         if (k .eq. -1) go to 10
         if (k .eq. lstgrd) then
            numHoods(i,j) =  numHoods(i,j) + 1
            go to 10 ! a full  flow cell is its own merged neighborhood
         endif
         vqmerge = 0.d0
         firstTimeThru = .true.
         nco = 0   ! initial size of neighborhood, from -1 to 1 square centered on cell

            do while (vqmerge < areaMin) 
               do 15 joff = -nco, nco
               do 15 ioff = -nco, nco
                   if (IS_GHOST(i+ioff,j+joff)) go to 15  
                   koff = irr(i+ioff,j+joff)
                   if (koff .eq. -1) go to 15  ! solid cells dont help
                   vqmerge = vqmerge + ar(koff)
                   if (firstTimeThru) then ! count everybody
                      numHoods(i+ioff,j+joff)=numHoods(i+ioff,j+joff)+1
                   else ! only add to new cells-on newly enlarged nhood border
                     if (abs(ioff).eq. nco .or. abs(joff).eq.nco)
     .                numHoods(i+ioff,j+joff)=numHoods(i+ioff,j+joff)+1
                   endif
 15             continue
                if (vqmerge >= areaMin*(1.d0-eps)) then
                   ncount(i,j) = nco
                   volMerge(i,j) = vqmerge ! in DRD it is reg. volume, not div by num nhoods
                   go to 10
                else   ! redo with larger neighborhood
                   nco = nco + 1
                   vqmerge = 0.d0
                   firstTimeThru = .false.
                endif
            end do
             
 10   continue      


c     needed number of neighbhoods to compute volMerge = which is not
c     the real volume of the merging neighborhood
c

!   initialize array with  most common case, overwritten below as needed
      xcentMerge = 0.d0      
      ycentMerge = 0.d0

      do 20 j = 1, mjtot
      do 20 i = 1, mitot
         k = irr(i,j)  
         if (k .eq. lstgrd) then
             call getCellCentroid(lstgrd,i,j,xcentMerge(i,j),
     .                            ycentMerge(i,j),
     .                            xlow,ylow,dx,dy,k)
             go to 20 ! a full  flow cell is its own merged neighborhood
         endif
         if (IS_GHOST(i,j) .or. k .eq. -1) then
            volMerge(i,j) = 0.d0
            go to 20
         endif

         xcent = 0.d0
         ycent = 0.d0
         nco = ncount(i,j)
         do 25 joff = -nco, nco
         do 25 ioff = -nco, nco
            if (IS_GHOST(i+ioff,j+joff)) go to 25  
            koff = irr(i+ioff,j+joff)
            if (koff .eq. -1) go to 25 ! solid cells dont help
            call getCellCentroid(lstgrd,i+ioff,j+joff,xc,yc,xlow,ylow,
     .                           dx,dy,koff)
            xcent = xcent + xc*ar(koff)
            ycent = ycent + yc*ar(koff)

 25      continue
         xcentMerge(i,j) = xcent/volMerge(i,j)
         ycentMerge(i,j) = ycent/volMerge(i,j)
             
 20   continue



      return
      end
c
c -------------------------------------------------------------------
c
      subroutine distribute(q,delta,volMerge,numHoods,ncount,
     .                      irr,mitot,mjtot,nvar,lstgrd,
     .                      denvolMerge)

      use amr_module
      implicit double precision (a-h, o-z)

      dimension numHoods(mitot,mjtot), volMerge(mitot,mjtot)
      dimension denvolMerge(mitot,mjtot)
      dimension ncount(mitot,mjtot), irr(mitot,mjtot)
      dimension q(mitot,mjtot,nvar), delta(mitot,mjtot,nvar)
      dimension qold(mitot,mjtot,1)
      logical IS_GHOST

      IS_GHOST(i,j) = (i .le. nghost .or. i .gt. mitot-nghost .or.
     .                 j .le. nghost .or. j .gt. mjtot-nghost)


      qold(:,:,1) = q(:,:,1)  ! save for density weighted distrib
      do j = 1, mjtot
      do i = 1, mitot
          k = irr(i,j)
          if (k .eq. -1) cycle
          if (IS_GHOST(i,j)) cycle
          ! only cut cells that merged have something to distribute
          if (k .eq. lstgrd) cycle  
          nco = ncount(i,j)
          if (nco .eq. 0) then
             if (delta(i,j,1) .ne. 0.d0) then
                write(*,*)"error"
             endif
             cycle  ! large enough cut cells dont play either
          endif

          dsum = 0.d0
          do 10 joff = -nco, nco
          do 10 ioff = -nco, nco
            koff = irr(i+ioff,j+joff)
            if (koff .eq. -1 .or. IS_GHOST(i+ioff,j+joff)) go to 10
            q(i+ioff,j+joff,:) = q(i+ioff,j+joff,:)-
     &                           delta(i,j,:)/volMerge(i,j)
c           dsum = dsum + ar(koff)*qold(i+ioff,j+joff,1)
c           q(i+ioff,j+joff,:) = q(i+ioff,j+joff,:)- delta(i,j,:)*
c    &                         qold(i+ioff,j+joff,1)/denvolMerge(i,j)
 10       continue
c         if (dabs(dsum-denvolMerge(i,j)) .gt. 1.d-14) then
c            write(*,*)"error here"
c         endif

      end do
      end do

      return
      end