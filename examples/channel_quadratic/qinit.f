
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

c
       pl = 1.0d0
       ul = 0.1d0
       vl = 0.01d0
c       ul = 0.0d0
c       vl = 0.0d0

c      fill ghost cells too so can more easily plot before
c      time stepping starts
       do 20 i = 1-mbc, mx+mbc
       do 20 j = 1-mbc, my+mbc
          kuse = irr(i,j)
          if (kuse .ne. -1 .and. kuse .ne. lstgrd) then
             xcen = xcirr(kuse)
             ycen = ycirr(kuse)
          else
             xcen = xcorn + (i-0.5d0)*dx
             ycen = ycorn + (j-0.5d0)*dy
          endif

          !rhol = ycen - 0.1d0*xcen + 0.5d0
          rhol = quad(xcen,ycen)
          !rhol = - 0.1d0*xcen + 0.5d0
          !rhol = 1.d0
          if (kuse .eq. -1) then  ! for easier debugging
            rhol = 1.4d0 
          endif

          q(1,i,j) = rhol
          q(2,i,j) = rhol * ul
          q(3,i,j) = rhol * vl
          q(4,i,j) = pl/.4d0 + 0.5d0*rhol*(ul**2 + vl**2)
  20   continue

       return
       end
