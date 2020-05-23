c
c -----------------------------------------------------------
c
c barth jespersen on cut cells and on their neighbors
c
      subroutine limitCellBJ(qp,qx,qy,mitot,mjtot,irr,nvar,hx,hy,
     .                    lwidth,lstgrd,xlow,ylow)
       
       use amr_module
       implicit double precision(a-h,o-z)
       dimension qp(nvar,mitot,mjtot),qx(nvar,mitot,mjtot),
     &           qy(nvar,mitot,mjtot),irr(mitot,mjtot)

      logical quad, nolimiter
      common /order2/ ssw, quad, nolimiter
      include "cuserdt.i"

      logical   prflag, OUT_OF_RANGE

      dimension dumax(nvar),dumin(nvar),phimin(nvar)
      dimension graddot(nvar),alpha(nvar),recon(nvar)
      logical ALL_NBORS_EXIST, IS_GHOST 

      IS_GHOST(i,j) = (i .le. lwidth .or. i .gt. mitot-lwidth .or.
     .                 j .le. lwidth .or. j .gt. mjtot-lwidth)


      OUT_OF_RANGE(i,j) = ((i .lt. 1) .or. (i.gt. mitot) .or.
     .                     (j .lt. 1) .or. (j.gt. mjtot))


      ALL_NBORS_EXIST(i,j) = (i .gt. 1 .and. i .lt. mitot .and.
     .                        j .gt. 1 .and. j .lt. mjtot)

c
c  **************************************
c
c   LIMIT CUT CELLS AND NEXT DOORS USING LINEARITY PRESERVING BARTH-JESPERSEN OVER RECONSTRUCTION NEIGHBORHOOD
c

      do 30 j = lwidth-1, mjtot-lwidth+2
      do 30 i = lwidth-1, mitot-lwidth+2
         k = irr(i,j)
         if (k .eq. -1) cycle ! skips solid/regular cells
         if (k .eq. lstgrd .and. ALL_NBORS_EXIST(i,j)) then
           if (irr(i+1,j) .eq. lstgrd .and.
     .         irr(i,j+1) .eq. lstgrd .and.
     .         irr(i,j-1) .eq. lstgrd .and.
     .         irr(i-1,j) .eq. lstgrd) go to 30
         endif
         if (ar(k)/ar(lstgrd) .lt. gradThreshold) then
            ! leave 0 gradient:  more stable for teeny cells w/o slopes
            qx(:,i,j) = 0.d0
            qy(:,i,j) = 0.d0
            go to 30          
         endif
c        write(*,900) i,j,qx(1,i,j),qy(1,i,j)
 900     format(2i4,2e15.7)

         call getCellCentroid(lstgrd,i,j,xcent,ycent,xlow,
     .                         ylow,hx,hy,k)


          ! find max and min needed for BJ limiting
          dumax = 0.d0
          dumin = 0.d0
          if (k .eq. lstgrd) then 
             nco = 1
          else
             nco = 2
          endif

          do 31 joff =  -nco, nco
          do 31 ioff =  -nco, nco
            if (ioff .eq. 0 .and. joff .eq. 0) go to 31
            if (OUT_OF_RANGE(i+ioff,j+joff)) go to 31
c           call getCellCentroid(lstgrd,i,j,xcent,ycent,xlow,
c    .                         ylow,hx,hy,k)
c           if (OUT_OF_RANGE(xcent,ycent)) go to 31
c           if (IS_GHOST(i,j)) go to 31
            koff = irr(i+ioff,j+joff)
            if (koff .eq. -1) go to 31

            dumax = max(dumax,qp(:,i+ioff,j+joff)-qp(:,i,j))
            dumin = min(dumin,qp(:,i+ioff,j+joff)-qp(:,i,j))
  31      continue

          phimin = 1.d0
          do 32 joff = -nco, nco
          do 32 ioff = -nco, nco
c            if (ioff .eq. 0 .and. joff .eq. 0) go to 32 ! no eqn to solve
             if (OUT_OF_RANGE(i+ioff,j+joff)) go to 32
c            if (IS_GHOST(i,j)) go to 32
             koff = irr(i+ioff,j+joff)
             if (koff .eq. -1) go to 32
             call getCellCentroid(lstgrd,i+ioff,j+joff,xc,yc,xlow,
     .                            ylow,hx,hy,koff)
              diffx = xc-xcent
              diffy = yc-ycent
              graddot  = qx(:,i,j)*diffx + qy(:,i,j)*diffy
              recon = qp(:,i,j) + graddot
                do m = 1,4
                   if (graddot(m) > 1.d-10) then
                      alpha(m) = min(1.d0, dumax(m)/graddot(m))
                   else if (graddot(m) < -1.d-10) then
                      alpha(m) = min(1.d0, dumin(m)/graddot(m))
                   else
                      alpha(m) = 1.d0
                   endif
                end do
                if (recon(1) .le. 0.d0) alpha(1) = 0.d0
                if (recon(4) .le. 0.d0) alpha(4) = 0.d0
                phimin = min(phimin, alpha)
  32      continue
          qx(:,i,j) = qx(:,i,j)*phimin(:)
          qy(:,i,j) = qy(:,i,j)*phimin(:)
c        write(*,901) qx(1,i,j),qy(1,i,j)
 901     format(8x,2e15.7)

         go to 30 ! skip testing for positivity at quadrature points for now
         ! only testing at neighboring cells above

      ! check for positivity
      if (irr(i,j) .eq. lstgrd) then
         alpha = 1.d0

         diffx = hx/2.d0
         diffy = 0.
         graddot  = qx(:,i,j)*diffx + qy(:,i,j)*diffy
         recon = qp(:,i,j) + graddot



         if (recon(1) .le. 0.d0) alpha(1) = 0.d0
         if (recon(4) .le. 0.d0) alpha(4) = 0.d0

         diffx = -hx/2.d0
         diffy = 0.
         graddot  = qx(:,i,j)*diffx + qy(:,i,j)*diffy
         recon = qp(:,i,j) + graddot


         if (recon(1) .le. 0.d0) alpha(1) = 0.d0
         if (recon(4) .le. 0.d0) alpha(4) = 0.d0

         diffx = 0.d0
         diffy = hy/2.d0
         graddot  = qx(:,i,j)*diffx + qy(:,i,j)*diffy
         recon = qp(:,i,j) + graddot


         if (recon(1) .le. 0.d0) alpha(1) = 0.d0
         if (recon(4) .le. 0.d0) alpha(4) = 0.d0

         diffx = 0.d0
         diffy = -hy/2.d0
         graddot  = qx(:,i,j)*diffx + qy(:,i,j)*diffy
         recon = qp(:,i,j) + graddot


         if (recon(1) .le. 0.d0) alpha(1) = 0.d0
         if (recon(4) .le. 0.d0) alpha(4) = 0.d0



         qx(1,i,j) = qx(1,i,j)*alpha(1)
         qy(4,i,j) = qy(4,i,j)*alpha(4)
      else
      ! positivity at quadrature points is required, midpoints because p is always = 1.
         do 510 kside=1,6
               x1 = poly(kside,1,k)
               y1 = poly(kside,2,k)
               x2 = poly(kside+1,1,k)
               y2 = poly(kside+1,2,k)

               d_bxp = (x1+x2)/2.d0
               d_byp = (y1+y2)/2.d0

               diffx = d_bxp-xcent
               diffy = d_byp-ycent
               graddot  = qx(:,i,j)*diffx + qy(:,i,j)*diffy
               recon = qp(:,i,j) + graddot


              alpha = 1.d0
              if (recon(1) .le. 0.d0) alpha(1) = 0.d0
              if (recon(4) .le. 0.d0) alpha(4) = 0.d0


              qx(1,i,j) = qx(1,i,j)*alpha(1)
              qy(4,i,j) = qy(4,i,j)*alpha(4)
              if (poly(kside+2,1,k).eq.-11.) exit
 510      continue
      endif

  30  continue

      return
      end
