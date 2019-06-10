c
c---------------------------------------------------------------------
c
c
      subroutine irregflux(q,firreg,irr,mitot,mjtot,dx,dy,
     &                 lstgrd,xlow,ylow,mptr,qx,qy,lwidth,nvar)

      use amr_module

      implicit double precision (a-h,o-z)

      dimension q(nvar,mitot,mjtot), irr(mitot,mjtot)
      dimension qx(nvar,mitot,mjtot),qy(nvar,mitot,mjtot)
      dimension firreg(nvar,-1:irrsize)
c
c     # Compute irregular flux in cut cells. Solve Riemann problems
c     # assume there are less than mtsize irregular cell sides on 
c     # any particular grid.
c
c ###
     
      logical  debug,quad, nolimiter
      integer  icell,jcell

      common  /order2/  ssw, quad, nolimiter

      data  debug/.false./
c
c
c       ========================================================
c
c     # loop around boundary: compute cut cell flux
c
c       =======================================================
c
      firreg(:,lstgrd) = 0.d0
      k = lstgrd
      k = iabs(nxtirr(k))

      do while (k .ne. 0) 
c
c     # march around the sides of this cell to find the irregular side
c
         do 20 kside=1,6
            if (poly(kside+2,1,k).eq.-11.) then
               x1 = poly(kside,1,k)
               y1 = poly(kside,2,k)
               x2 = poly(kside+1,1,k)
               y2 = poly(kside+1,2,k)
               go to 25 
            endif
 20      continue
 25      continue
c
c        # boundary segment face:
c        ------------------------
c
         ix0 = ixg(k)
         iy0 = iyg(k)
c     # compute (alf,beta) = unit normal to boundary pointing in.
         hsx1 = x2
         hsx2 = x1
         hsy1 = y2
         hsy2 = y1
         rlen = dsqrt((hsy1-hsy2)**2 + (hsx1-hsx2)**2)
         alf = (hsy1-hsy2)/rlen
         beta = (hsx2-hsx1)/rlen
c     
c     # compute boundary data  - q already has primitive variables 
c     
         rho_c    = q(1,ix0,iy0)
         u_c      = q(2,ix0,iy0)
         v_c      = q(3,ix0,iy0)
         pr_c     = q(4,ix0,iy0)

         bxpt = .5d0*(hsx1+hsx2)
         bypt = .5d0*(hsy1+hsy2)
         x0   = xcirr(k)
         y0   = ycirr(k)
         xdif = (bxpt-x0)
         ydif = (bypt-y0)


         rho = rho_c + xdif*qx(1,ix0,iy0) + ydif*qy(1,ix0,iy0)
         u   = u_c +   xdif*qx(2,ix0,iy0) + ydif*qy(2,ix0,iy0)
         v   = v_c +   xdif*qx(3,ix0,iy0) + ydif*qy(3,ix0,iy0)
         pr  = pr_c +  xdif*qx(4,ix0,iy0) + ydif*qy(4,ix0,iy0)

c     
c     NEW WAY: NO RP, JUST PRESSURE AT BNDRY, NO RP
c     
         firreg(1,k) = 0.d0
         firreg(2,k) = rlen*alf*pr
         firreg(3,k) = rlen*beta*pr
         firreg(4,k) = 0.d0
c
         k = iabs(nxtirr(k))
      end do
c     
      return
      end      
