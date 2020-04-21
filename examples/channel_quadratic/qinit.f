
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
c      fill ghost cells too so can more easily plot before
c      time stepping starts
       do 20 i = 1-mbc, mx+mbc
       do 20 j = 1-mbc, my+mbc
          kuse = irr(i,j)
          if (kuse .ne. -1) then
            if (kuse .ne. lstgrd) then
               xcen = xcirr(kuse)
               ycen = ycirr(kuse)
            else
               xcen = xcorn + (i-0.5d0)*dx
               ycen = ycorn + (j-0.5d0)*dy
               call makep(poly(1,1,lstgrd),i,j,xcorn,ycorn,dx,dy)
            endif
c           call channelPtInit(xcen,ycen,rhol,ul,vl,pl)
            call channelAvgInit(poly(1,1,kuse),xcen,ycen,rhol,ul,vl,pl)
          else
            rhol = 1.4d0
            ul = 0.d0
            vl = 0.d0
            pl = 1.d0
          endif

          q(1,i,j) = rhol
          q(2,i,j) = rhol * ul
          q(3,i,j) = rhol * vl
          q(4,i,j) = pl/.4d0 + 0.5d0*rhol*(ul**2 + vl**2)
  20   continue

       return
       end
