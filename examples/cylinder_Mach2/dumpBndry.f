c
c ----------------------------------------------------------
c
      subroutine dumpBndry(qp,qx,qy,irr,mitot,mjtot,i,j,nvar,
     &                     ibunit,mptr)

      use amr_module
      implicit double precision (a-h,o-z)

      dimension qp(nvar,mitot,mjtot) 
      dimension qx(nvar,mitot,mjtot),qy(nvar,mitot,mjtot)
      integer irr(mitot,mjtot) 
      

      pi = 3.14159265358979d0
      pi2 = pi * 2.d0

      k = irr(i,j)
      ! get midpoint of bndry face
      do 20 kside=1,6
        if (poly(kside+2,1,k).eq.-11.) then
           x1 = poly(kside,1,k)
           y1 = poly(kside,2,k)
           x2 = poly(kside+1,1,k)
           y2 = poly(kside+1,2,k)
           go to 25
        endif
 20   continue

 25   continue

      bxpt = .5d0*(x2+x1)
      bypt = .5d0*(y2+y1)
      x0   = xcirr(k)
      y0   = ycirr(k)
      xdif = (bxpt-x0)
      ydif = (bypt-y0)

c     reconstruct
      rhob = qp(1,i,j) + xdif*qx(1,i,j) + ydif*qy(1,i,j)
      ub   = qp(2,i,j) + xdif*qx(2,i,j) + ydif*qy(2,i,j)
      vb   = qp(3,i,j) + xdif*qx(3,i,j) + ydif*qy(3,i,j)
      pb   = qp(4,i,j) + xdif*qx(4,i,j) + ydif*qy(4,i,j)

      ! center of circle is 2,2. 
      ! compute angle in radians for plotting
      yd = bypt - 0.5d0
      xd = bxpt - 0.5d0
      angle = atan2(yd,xd)
c
c     shift to match previous paper's ragne
c     if (angle  .lt. 0.d0) then
c          angle = angle + pi2
c     endif

      write(ibunit,100) angle,bxpt,bypt,rhob,ub,vb,pb,i,j,mptr
 100  format(7e12.4,3i4)
 
      return
      end
