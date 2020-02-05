c
c------------------------------------------------------------
c
       subroutine slopes(qp,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar,mptr,istage)

       use amr_module
       implicit double precision(a-h,o-z)
       dimension qp(nvar,mitot,mjtot),qx(nvar,mitot,mjtot),
     &           qy(nvar,mitot,mjtot),irr(mitot,mjtot)

       dimension dumax(nvar),dumin(nvar),phimin(nvar), graddot(nvar)
       dimension dalpha(nvar), recon(nvar), b(2)
       logical  regular, quad, nolimiter,pwconst
       logical MC/.true./
       common /userdt/ cflcart,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold,pwconst
       common /order2/ ssw, quad, nolimiter

       logical IS_GHOST
       dimension a(2,2), rhs(2,nvar)

c
c      # ssw = slope switch (1. for slopes, 0 for donor cell 0 slopes)
c      # now set in amrcart
c
c   ssw = 0 for ctu (no slopes), 1 for muscl (second order).
c   compute slopes using muscl limiter
c   qp contains the primitive variables
c
      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                  j .le. lwidth .or. j .gt. mjtot-lwidth)
c
      isloperecon = 4
      !isloperecon = 2
      if (isloperecon .eq. 4 .and. nolimiter) then
         isloperecon = 5
      endif
      ! 0 is no slopes
      ! 1 is central
      ! 2 is MC
      ! 3 is LSQ + BJ (8 neighbors)
      ! 4 is LSQ + BJ (4 neighbors)
      ! 5 is LSQ + no limiter (4 neighbors)

      if(isloperecon .eq. 0) then ! no slopes
          qx = 0.d0
          qy = 0.d0

      else if(isloperecon .eq. 1) then ! central
       do 10 j = 1, mjtot
       do 10 i = 2, mitot-1
       do 10 m = 1, nvar
          if (irr(i,j) .ne. lstgrd) goto 10
          ducc = qp(m,i+1,j) - qp(m,i-1,j)
          qx(m,i,j) = .5d0*ducc/hx
 10    continue

       do 20 j = 2, mjtot-1
       do 20 i = 1, mitot
       do 20 m = 1, nvar
       if (irr(i,j) .ne. lstgrd) goto 20
          ducc = qp(m,i,j+1) - qp(m,i,j-1)
          qy(m,i,j) = .5d0*ducc/hy
 20    continue

      else if(isloperecon .eq. 2) then ! MC

       do 30 j = 1, mjtot
       do 30 i = 2, mitot-1
       do 30 m = 1, nvar
          if (irr(i,j) .ne. lstgrd) goto 30
          ducc = qp(m,i+1,j) - qp(m,i-1,j)
          dupc = qp(m,i+1,j) - qp(m,i,j)
          dumc = qp(m,i,j)   - qp(m,i-1,j)
          du   = dmin1(dabs(dupc),dabs(dumc))
          du   = dmin1(2.d0*du, .5d0*dabs(ducc))
          fl = dmax1(0.d0,dsign(1.d0, dupc*dumc))*ssw
          qx(m,i,j) = du*dsign(1.d0,ducc)*fl/hx
 30    continue
       do 40 j = 2, mjtot-1
       do 40 i = 1, mitot
       do 40 m = 1, nvar
          if (irr(i,j) .ne. lstgrd) goto 40
          ducc = qp(m,i,j+1) - qp(m,i,j-1)
          dupc = qp(m,i,j+1) - qp(m,i,j)
          dumc = qp(m,i,j)   - qp(m,i,j-1)
          du   = dmin1(dabs(dupc),dabs(dumc))
          du   = dmin1(2.d0*du, .5d0*dabs(ducc))
          fl = dmax1(0.d0,dsign(1.d0, dupc*dumc))*ssw
          qy(m,i,j) = du*dsign(1.d0,ducc)*fl/hy
 40    continue

      else if (isloperecon .eq. 3 .or. isloperecon .eq. 4
     .    .or. isloperecon .eq. 5) then
      ! LSQ (4 or 8 neighbors) with/without limiting
       do 34 j = 2, mjtot-1
       do 34 i = 2, mitot-1
            if (irr(i,j) .ne. lstgrd) goto 34
            rhs = 0.d0 ! initialize for accumulation
            a = 0.d0
            k = irr(i,j)
            if(k .eq. lstgrd) then
                x0 = xlow + (i-.5d0)*hx
                y0 = ylow + (j-.5d0)*hy
            else
                x0 = xcirr(k)
                y0 = ycirr(k)
            endif

            do 22 joff = -1, 1
            do 22 ioff = -1, 1
                if (ioff .eq. 0 .and. joff .eq. 0) goto 22 ! no eqn to solve

                if (isloperecon .eq. 4 .or. isloperecon .eq. 5) then ! only 4 neighbors
                    if (abs(joff).eq.1 .and. abs(ioff).eq.1) goto 22
                endif

                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) goto 22

                if(koff .eq. lstgrd) then
                    xoff = xlow + (i+ioff-.5d0)*hx
                    yoff = ylow + (j+joff-.5d0)*hy
                else
                    xoff = xcirr(koff)
                    yoff = ycirr(koff)
                endif

               deltax =  xoff - x0
               deltay =  yoff - y0
               a(1,1) = a(1,1) + deltax*deltax
               a(1,2) = a(1,2) + deltax*deltay
               a(2,2) = a(2,2) + deltay*deltay
               rhs(1,:) = rhs(1,:) + deltax *
     .                    (qp(:,i+ioff,j+joff) - qp(:,i,j))
               rhs(2,:) = rhs(2,:) + deltay *
     .                    (qp(:,i+ioff,j+joff) - qp(:,i,j))
 22          continue

             ! solve a^t*a  * grad = a^t rhs.  First form cholesky factors.
             ! will have to add robustness checks
             c11 = sqrt(a(1,1))
             c12 = a(1,2) / c11
             c22 = sqrt(a(2,2) - c12**2)

             ! now back solve (C^t C = rhs of A^tdu ) to get x and y gradient for all variables
             do mm = 1, nvar
               b(1) = rhs(1,mm) / c11
               b(2) = (rhs(2,mm) - c12*b(1))/c22
               qy(mm,i,j) = b(2) / c22
               qx(mm,i,j) = (b(1) - c12*qy(mm,i,j))/c11
             end do

 34    continue
      endif

      if(isloperecon .eq. 3 .or. isloperecon .eq. 4) then
        do 50 i = 2, mitot-1
        do 50 j = 2, mjtot-1
            k = irr(i,j)
            if (k .eq. -1) goto 50 ! solid cells have no gradient

            ! find max and min needed for BJ limiting
            dumax = 0.d0
            dumin = 0.d0
            do 31 joff = -1, 1
            do 31 ioff = -1, 1
              if (ioff .eq. 0 .and. joff .eq. 0) goto 31
c              if (IS_GHOST(i+ioff,j+joff)) goto 31

                if(isloperecon .eq. 4) then
                    if(abs(joff) .eq. 1 .and. abs(ioff) .eq. 1) goto 31
                endif

              koff = irr(i+ioff,j+joff)
              if (koff .eq. -1) goto 31 ! ignore solid cells
              dumax = max(dumax,qp(:,i+ioff,j+joff)-qp(:,i,j))
              dumin = min(dumin,qp(:,i+ioff,j+joff)-qp(:,i,j))
 31         continue

            phimin = 1.d0

            do 93 iface = 1,4
                if(iface .eq. 1) then
                    diffx = 0.5d0*hx
                    diffy = 0.d0
                else if(iface .eq. 2) then
                    diffx = 0.d0
                    diffy = 0.5d0*hy
                else if(iface .eq. 3) then
                    diffx = -0.5d0*hx
                    diffy = 0.d0
                else if(iface .eq. 4) then
                    diffx = 0.d0
                    diffy = -0.5d0*hy
                endif

                graddot  = qx(:,i,j)*diffx + qy(:,i,j)*diffy

                do m = 1,4
                   if (graddot(m) > 1.e-16) then
                      dalpha(m) = min(1.d0, dumax(m)/graddot(m))
                   else if (graddot(m) < -1.e-16) then
                      dalpha(m) = min(1.d0, dumin(m)/graddot(m))
                   else
                      dalpha(m) = 1.d0
                   endif
                end do

                phimin = min(phimin, dalpha)

 93         continue

            do 934 iface = 1,4
                if(iface .eq. 1) then
                    diffx = 0.5d0*hx
                    diffy = 0.d0
                else if(iface .eq. 2) then
                    diffx = 0.d0
                    diffy = 0.5d0*hy
                else if(iface .eq. 3) then
                    diffx = -0.5d0*hx
                    diffy = 0.d0
                else if(iface .eq. 4) then
                    diffx = 0.d0
                    diffy = -0.5d0*hy
                endif

                graddot  = qx(:,i,j)*diffx + qy(:,i,j)*diffy
                recon = qp(:,i,j) + phimin * graddot
                ! one last check for positivity
                if (recon(1) .le. 0.d0 .or. recon(4) .le. 0.d0) then
                  dalpha = 0.d0
                  write(*,900) i,j,recon(1),recon(4),mptr,istage
 900              format("domain nonphysical: cell ",2i5,
     .            " rho, pr",2e15.7," grid ",i4," istage ",i2)
                endif
                phimin = min(phimin, dalpha)
 934         continue


            qx(:,i,j) = qx(:,i,j)*phimin(:)
            qy(:,i,j) = qy(:,i,j)*phimin(:)
 50     continue
      endif


 99    continue
       return
       end
