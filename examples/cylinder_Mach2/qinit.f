
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
c      mach 2 shock 
c
       rhol = 1.4d0*8.d0/3.d0
       ul   = 14.d0/3.d0  / rhol
       vl   = 0.00d0
       pl   = 4.5d0  


       rhor = 1.4d0
       ur   = 0.d0
       vr   = 0.d0
       pr   = 1.0d0

       sloc = 0.2d0

       do 20 i=1,mx
       do 20 j=1,my
          kuse = irr(i,j)
          if (kuse .ne. -1 .and. kuse .ne. lstgrd) then
             xcen = xcirr(kuse)
             ycen = ycirr(kuse)
          else
             xcen = xcorn + (i-0.5d0)*dx
             ycen = ycorn + (j-0.5d0)*dy
          endif

          ! put in quadratic to look only at initial gradients
          ! q(1,i,j) = xcen**2 - 2*ycen**2 + 5.d0 + (xcen-1.d0)*ycen
          ! q(2,i,j) = 0.d0
          ! q(3,i,j) = 0.d0
          ! q(4,i,j) = 10.d0
          ! go to 20

          if (xcen .lt. sloc) then
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

 99    return
       end
