c
c ---------------------------------------------------------------------
c
       subroutine qslopes(qp,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,
     &                    hx,hy,xlow,ylow,mptr,nvar,istage)

      use amr_module
      implicit double precision(a-h,o-z)

      dimension qp(nvar,mitot,mjtot),qx(nvar,mitot,mjtot),
     &          qy(nvar,mitot,mjtot),irr(mitot,mjtot)

      dimension a(2,2), rhs(2,nvar), b(2)
      dimension recon(nvar), graddot(nvar), dalpha(nvar)
      dimension dumin(nvar),dumax(nvar), phimin(nvar)
      common /order2/ ssw, quad, nolimiter
      common /userdt/cflcart,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold,pwconst

      logical IS_GHOST
      integer isloperecon/4/
      ! 3 is LSQ + BJ (8 neighbors)
      ! 4 is LSQ + BJ (4 neighbors)
      ! 5 is LSQ + no limiter (4 neighbors)
c
      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                 j .le. lwidth .or. j .gt. mjtot-lwidth)

        nco = 1
        do 20 j = 2, mjtot-1
        do 20 i = 2, mitot-1
            k = irr(i,j)
            if(k .eq. lstgrd .or. k .eq. -1) goto 20

            ! can reconstruct high order poly here, but for now do linear reconstruction
            ! solving least squares problem with all neighboring merged vals. Should preprocess the matrix and save
            rhs = 0.d0 ! initialize for accumulation
            a = 0.d0
            x0 = xcirr(k)
            y0 = ycirr(k)

            do 22 joff = -nco, nco
            do 22 ioff = -nco, nco
                if (ioff .eq. 0 .and. joff .eq. 0) go to 22 ! no eqn to solve
                if(isloperecon .eq. 4 .or. isloperecon .eq. 5) then
                    if(abs(joff) .eq. 1 .and. abs(ioff) .eq. 1) goto 22
                endif

c                if (IS_GHOST(i+ioff,j+joff)) go to 22
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 22

                if (koff .eq. lstgrd) then
                    xoff = xlow + (i+ioff-.5d0)*hx
                    yoff = ylow + (j+joff-.5d0)*hy
                else
                    xoff = xcirr(koff)
                    yoff = ycirr(koff)
                endif

                deltax = xoff - x0
                deltay = yoff - y0
                a(1,1) = a(1,1) + deltax*deltax
                a(1,2) = a(1,2) + deltax*deltay
                a(2,2) = a(2,2) + deltay*deltay
                rhs(1,:) = rhs(1,:) + deltax *
     .                     (qp(:,i+ioff,j+joff) - qp(:,i,j))
                rhs(2,:) = rhs(2,:) + deltay *
     .                     (qp(:,i+ioff,j+joff) - qp(:,i,j))
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
               qy(m,i,j) = b(2) / c22
               qx(m,i,j) = (b(1) - c12*qy(m,i,j))/c11
             end do


            if(isloperecon .eq. 3 .or. isloperecon .eq. 4) then

            ! find max and min needed for BJ limiting
            dumax = 0.d0
            dumin = 0.d0
            do 31 joff = -1, 1
            do 31 ioff = -1, 1
              if (ioff .eq. 0 .and. joff .eq. 0) go to 31
c              if (IS_GHOST(i+ioff,j+joff)) go to 31

              if(isloperecon .eq. 4) then
                  if(abs(joff) .eq. 1 .and. abs(ioff) .eq. 1) goto 31
              endif

              koff = irr(i+ioff,j+joff)
              if (koff .eq. -1) go to 31
              dumax = max(dumax,qp(:,i+ioff,j+joff)-qp(:,i,j))
              dumin = min(dumin,qp(:,i+ioff,j+joff)-qp(:,i,j))
 31         continue

            phimin = 1.d0

            do 93 kside = 1,6
               if (poly(kside+1,1,k) .ne. -11)then
                  sidex2 = poly(kside,1,k)
                  sidey2 = poly(kside,2,k)
                  sidex1 = poly(kside+1,1,k)
                  sidey1 = poly(kside+1,2,k)

                  dmidx = 0.5d0*(sidex1 + sidex2)
                  dmidy = 0.5d0*(sidey1 + sidey2)

                  diffx = dmidx - xcirr(k)
                  diffy = dmidy - ycirr(k)
               else
                    exit
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

            do 934 kside = 1,6
               if (poly(kside+1,1,k) .ne. -11)then
                  sidex2 = poly(kside,1,k)
                  sidey2 = poly(kside,2,k)
                  sidex1 = poly(kside+1,1,k)
                  sidey1 = poly(kside+1,2,k)

                  dmidx = 0.5d0*(sidex1 + sidex2)
                  dmidy = 0.5d0*(sidey1 + sidey2)

                  diffx = dmidx - xcirr(k)
                  diffy = dmidy - ycirr(k)
               else
                    exit
               endif

                graddot  = qx(:,i,j)*diffx + qy(:,i,j)*diffy
                recon = qp(:,i,j) + phimin * graddot
                ! one last check for positivity
                if (recon(1) .le. 0.d0 .or. recon(4) .le. 0.d0) then
                  dalpha = 0.d0
                  write(*,900) i,j,mptr,istage
 900              format("boundary nonphysical: cell ",2i5,
     .                   " grid ",i4," istage ",i2)
                endif
                phimin = min(phimin, dalpha)
 934         continue


            qx(:,i,j) = qx(:,i,j)*phimin(:)
            qy(:,i,j) = qy(:,i,j)*phimin(:)

       endif  ! slope recon is 3 or 4


 20    continue

       end
