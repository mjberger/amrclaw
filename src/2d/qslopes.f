c
c ---------------------------------------------------------------------
c
       subroutine qslopes(qp,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,
     &                    hx,hy,xlow,ylow,mptr,nvar,istage)

      use amr_module
      implicit double precision(a-h,o-z)

      dimension qp(nvar,mitot,mjtot),qx(nvar,mitot,mjtot),
     &          qy(nvar,mitot,mjtot),irr(mitot,mjtot)

      logical  pwconst
      common /order2/ ssw, quad, nolimiter
      include cuserdt.i 

      dimension rhsmax(nvar), rhsmin(nvar)
      dimension dumax(nvar),dumin(nvar),phimin(nvar)
      dimension a(25,5),at(5,25),c(5,5),b(25,nvar),d(25,nvar)
      dimension rhs(25,nvar), w(5,nvar)
      dimension nlist(25,2)
      logical   prflag, quad, enufNbor, outside, almost_outside
      data      prflag/.false./
      logical   nolimiter, ghost_ccg/.false./
      logical ALL_NBORS_EXIST
c
      outside(i,j) = ((i.lt.1) .or. (j.lt.1.) .or. 
     .                (i.gt.mitot) .or. (j.gt.mjtot))

      ALL_NBORS_EXIST(i,j) = (i .gt. 1 .and. i .lt.mitot .and.
     .                        j .gt. 1 .and. j .lt. mjtot)


c
c   ##########
c   #  compute slopes for cut cells using least squares approach
c   ##########
c   
c   ssw = 0 for ctu (no slopes), 1 for muscl (second order).
c   qp contains primitive variables, regular cell slopes already set 
c   all slopes initialized to zero in regular slopes subr.
c
c
      quad = .false.
      nterms = 2 


c      do 110 iy0 = lwidth-2, mjtot-lwidth+3
c      do 110 ix0 = lwidth-2, mitot-lwidth+3
      do 110 iy0 = 1, mjtot
      do 110 ix0 = 1, mitot
         k = irr(ix0,iy0)
         if (k .eq. -1) go to 110
         if (k .eq. lstgrd .and. ALL_NBORS_EXIST(ix0,iy0)) then
           if (irr(ix0+1,iy0) .eq. lstgrd .and. 
     .         irr(ix0,iy0+1) .eq. lstgrd .and. 
     .         irr(ix0,iy0-1) .eq. lstgrd .and. 
     .         irr(ix0-1,iy0) .eq. lstgrd) go to 110
         endif
c     
         if (ar(k)/ar(lstgrd) .lt. gradThreshold) then
            qx(:,ix0,iy0) = 0.d0
            qy(:,ix0,iy0) = 0.d0
            go to 110           ! leave 0 gradient:  more stable for teeny cells w/o slopes
         endif
c     
c      # this cell needs a slope
         if (k .ne. lstgrd) then
            x0 = xcirr(k)
            y0 = ycirr(k)
            do 37 kside = 1, 6
               if (poly(kside+2,1,k).eq.-11)then
                  sidex2 = poly(kside,1,k)
                  sidey2 = poly(kside,2,k)
                  sidex1 = poly(kside+1,1,k)
                  sidey1 = poly(kside+1,2,k)
                  go to 39
               endif
 37         continue
 39         rlenb = dsqrt((sidey1-sidey2)**2 + (sidex1-sidex2)**2)
            alf  = (sidey1-sidey2)/rlenb
            beta = (sidex2-sidex1)/rlenb
            bxpt = .5d0*(sidex1+sidex2)
            bypt = .5d0*(sidey1+sidey2)
         else
            x0 = xlow + (ix0-.5d0)*hx
            y0 = ylow + (iy0-.5d0)*hy
         endif
         
c     if (outside(x0,y0)) go to 110   ! to match cart3d, no gradients in ghost cells
         nlist(1,1) = ix0
         nlist(1,2) = iy0
         nst        = 1  
         nend       = 1  
         call addneighbors(irr,nlist,nst,nend,newend,mitot,mjtot,
     .        lstgrd,quad,xlow,ylow,hx,hy)
         if (newend .gt. 3) go to 16 ! since have at least 2 nbors
c     in the linear case 2 nbors could be ill conditioned. add diagonal nbor
c     first have to find diag nbor. this code doesnt test for thin bodies,etc
c        ixn = nlist(2,1)
c        iyn = nlist(2,2)
c        if (nlist(2,1) .eq. ix0) ixn = nlist(3,1)
c        if (nlist(2,2) .eq. iy0) iyn = nlist(3,2)
c        newend = newend + 1
c        write(*,*)"qslopes ix0,iy0,ixn,iyn ",ix0,iy0,ixn,iyn
c        if (.not. outside(ixn,iyn)) then
c          nlist(newend,1) = ixn 
c          nlist(newend,2) = iyn
c          go to 16
c        else
           !not enough neighbors for gradient
            enufNbor = .false.
            go to 110
c        endif


c        ::: not enough neighbors for quadratic fit
c        ::: add neighbors of neighbors    
c
c         if ((irow .lt. 2)) then
c          write(*,*) "no slope for cut cell ",ix0,iy0
c          go to 110    ! not enough neighbors for gradient recon
c         endif

 16      irow = 0
         do 22 n = 2, newend
            irow = irow + 1
            ixn = nlist(n,1)
            iyn = nlist(n,2)
            kn =  irr(ixn,iyn)
            if (kn .ne. lstgrd) then
               xn = xcirr(kn)
               yn = ycirr(kn)
            else
               xn = xlow + (ixn-.5d0)*hx
               yn = ylow + (iyn-.5d0)*hy
            endif

c     # treating q as pt.wise values, wont work for quadratic fit
               a(irow,1) = (xn - x0)
               a(irow,2) = (yn - y0)

               do m = 1, nvar
                  b(irow,m) = qp(m,ixn,iyn) - qp(m,ix0,iy0)
               end do
 22         continue
c
c 
          do 30 it = 1, irow
          do 30 jt = 1, nterms
             at(jt,it) = a(it,jt)
 30       continue
c     
          do 40  m = 1, nvar
             rhsmax(m) = b(1,m)
             rhsmin(m) = b(1,m)
             do 40 it = 1, irow
                rhs(it,m) = b(it,m)
                rhsmax(m) = dmax1(rhsmax(m),b(it,m))
                rhsmin(m) = dmin1(rhsmin(m),b(it,m))
 40          continue
c     
             do 50 it = 1, nterms
                do 50 jt = 1, nterms
                   c(it,jt) = 0.d0
                   do m = 1, nvar
                      d(it,m)  = 0.d0
                   end do
                   do 45 kt = 1, irow
                      c(it,jt) = c(it,jt) + at(it,kt)*a(kt,jt)
                      do m = 1, nvar
                      d(it,m) = d(it,m) + at(it,kt)*b(kt,m)
                      end do
 45                continue
 50             continue

c
c
c      #  do linear fit

c     # solve C*w = d for least squares slopes. use cholesky
c     # put factors back in a
           a(1,1) = dsqrt(c(1,1))
           a(1,2) = c(1,2)/a(1,1)
           a(2,2) = dsqrt(c(2,2)-a(1,2)**2)
c     
           do 61 m = 1, nvar
c     
c     # at*a = c. solve at*b = d, aw = b.  reuse b.
              b(1,m) = d(1,m) / a(1,1)
              b(2,m) = (d(2,m) - a(1,2)*b(1,m)) / a(2,2)
              w2 =   b(2,m)/a(2,2)
              qy(m,ix0,iy0) =  w2
              qx(m,ix0,iy0) = (b(1,m)-a(1,2)*w2)/a(1,1)
 61        continue


c       should we limit?
        if (nolimiter)  go to 110

c

c      ::: now some kind of limiting of more accurate slope computed above
c      ::: fix up matrix first - transpose remained
c     
       do 65 i = 1, irow
         a(i,1) = at(1,i)
         a(i,2) = at(2,i)
 65    continue
c
c      ::::: check if solution feasible (doesn't need limiting)
c            see if exceed neighboring values at corners of cell
c
       if (k .eq. lstgrd) then  ! put in poly for this regular cell
          call makep(poly(1,1,lstgrd),ix0,iy0,xlow,ylow,hx,hy)
       endif

       if (.false.) then !limit at midpoints of cell edges
       do 95 m = 1, nvar
          frac = 1.d0
          do 66 iside = 1, 10
             if (poly(iside+1,1,k) .eq. -11) go to 67
             x = poly(iside,1,k)
             y = poly(iside,2,k)
             xdif = (x-x0)
             ydif = (y-y0)
             dd = xdif*qx(m,ix0,iy0)+ydif*qy(m,ix0,iy0)
             if (dd .gt. rhsmax(m)) frac=dmin1(frac,rhsmax(m)/dd)
             if (dd .lt. rhsmin(m)) frac=dmin1(frac,rhsmin(m)/dd)
             if ((dd .lt. 0.d0 .and. rhsmin(m) .gt. 0.d0) .or.
     .            (dd .gt. 0.d0 .and. rhsmax(m) .lt. 0.d0)) frac = 0.d0
 66       continue
 67       continue
       qx(m,ix0,iy0) = frac*qx(m,ix0,iy0)
       qy(m,ix0,iy0) = frac*qy(m,ix0,iy0)
c     
 95   continue
      endif
c   
c do it again to limit against neighboring cell centers. (as in minmod vs previous mc limiter)

       ! BJ limiting
       dumax = 0.0d0
       dumin = 0.0d0
       phimin = 1.d0
       do it = 1, irow
          ixn = nlist(1+it,1)
          iyn = nlist(1+it,2)
          dumax = max(dumax(:),qp(:,ixn,iyn) - qp(:,ix0,iy0))
          dumin = min(dumin(:),qp(:,ixn,iyn) - qp(:,ix0,iy0))
       end do
       do 107 it = 1, irow
         ixn = nlist(1+it,1)
         iyn = nlist(1+it,2)
         xdif = a(it,1)  ! (x-x0)
         ydif = a(it,2)  ! (y-y0)

         do 97 m = 1, nvar
             dd = xdif*qx(m,ix0,iy0)+ydif*qy(m,ix0,iy0)

             if (dd .lt. 0.d0) then
                frac = min(1.d0, dumin(m)/dd)
             else if (dd .gt. 0.d0) then
                frac = min(1.d0, dumax(m)/dd)
             else
                frac = 1.d0
             endif
           phimin(m) = min(frac,phimin(m))
 97        continue
 107      continue
          qx(:,ix0,iy0) = qx(:,ix0,iy0)*phimin(:)
          qy(:,ix0,iy0) = qy(:,ix0,iy0)*phimin(:)
c  
 110  continue
c
      if (prflag) then
         write(21,*)' qx '
         do 180 i = 2, mitot-1
         do 180 j = 2, mjtot-1
            if (irr(i,j) .ne. -1 .and. irr(i,j).ne.lstgrd) then
               write(21,190)i,j,(qx(m,i,j),m=1,nvar)
 190           format('  i,j  ',2i4,4e14.6)
            endif
 180     continue
         write(21,*)' qy '
         do 181 i = 2, mitot-1
         do 181 j = 2, mjtot-1
            if (irr(i,j) .ne. -1 .and. irr(i,j).ne.lstgrd) then
               write(21,190)i,j,(qy(m,i,j),m=1,nvar)
            endif
 181     continue
      endif
c
 99    return
       end
