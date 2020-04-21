
c
c
c
c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,xcorn,ycorn,
     &                   dx,dy,q,maux,aux,lstgrd,irr)
c     =====================================================
c
c     # Set initial conditions for q.
c     # rotated channel problem
c
       use amr_module
       implicit double precision (a-h,o-z)

       dimension q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
       dimension aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
       integer irr(1-mbc:mx+mbc, 1-mbc:my+mbc)

       common /moredata/ sloc

c
c      mach 10 flow 30 degree ramp
c
       pi = 3.14159265358979d0
c      alpha = 30.d0 * pi/180.d0
       alpha = pi/6.d0
       dxhalf = 0.5d0*dx
       dyhalf = 0.5d0*dy
       timeInit = 0.d0

       rhol = 8.d0
       ul   = cos(alpha)*8.25d0
       vl   = -sin(alpha)*8.25d0
       pl   = 116.5d0

       rhor = 1.4d0
       ur   = 0.d0
       vr   = 0.d0
       pr   = 1.0d0

c      do 20 i=1,mx
c      do 20 j=1,my
c      uncomment next 2 lines to initialize ghost cells too.
c      helps with plotting
       do 20 i=1-mbc,mx+mbc
       do 20 j=1-mbc,my+mbc
          kuse = irr(i,j)
          if (kuse .ne. -1 .and. kuse .ne. lstgrd) then
             xcen = xcirr(kuse)
             ycen = ycirr(kuse)
          else
             xcen = xcorn + (i-0.5d0)*dx
             ycen = ycorn + (j-0.5d0)*dy
          endif

          ! cellave approach, for regular cells only
          !call cellave(xcen-dxhalf,ycen-dyhalf,dx,dy,wl,timeInit)

          ! for this cells y, where is shock in x
          !rho = (1.d0-wl)*rhor + wl*rhol
          !u   = (1.d0-wl)*ur   + wl*ul
          !v   = (1.d0-wl)*vr   + wl*vl
          !p   = (1.d0-wl)*pr   + wl*pl
          !q(1,i,j) = rho
          !q(2,i,j) = rho * u
          !q(3,i,j) = rho * v
          !q(4,i,j) = p/.4d0 + 0.5d0*rho*(u*u+v*v)
          !go to 20

          !pointwise approach
          xshock =  sloc + ycen/sqrt(3.0d0)
          if (xcen .lt. xshock) then
            q(1,i,j) = rhol
            q(2,i,j) = rhol * ul
            q(3,i,j) = rhol * vl
            q(4,i,j) = pl/.4d0 + 0.5d0*rhol*(ul**2 + vl**2)
          else
            q(1,i,j) = rhor
            q(2,i,j) = rhor * ur
            q(3,i,j) = rhor * vr
            q(4,i,j) = pr/.4d0 + 0.5d0*rhor*(ur**2 + vr**2)
          endif
  20   continue

       return
       end
