
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

c  put these vals in for solid cells
       rhoSolid = 1.d0
       uSolid   = 0.0d0
       vSolid   = 0.0d0
       pSolid   = 1.0d0

c      use pointwise evaluation of exact solution
c      for higher order methods would have to compute integral average

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
          call ssvInit(xcen,ycen,rhot,ut,vt,pt)             
          if (kuse .eq. -1) then
            q(1,i,j) = rhoSolid
            q(2,i,j) = 0.d0
            q(3,i,j) = 0.d0
            q(4,i,j) = pSolid/.4d0 
          else
            q(1,i,j) = rhot
            q(2,i,j) = rhot * ut
            q(3,i,j) = rhot * vt
            q(4,i,j) = pt/.4d0 + 0.5d0*rhot*(ut**2 + vt**2)
          endif
  20   continue

       return
       end
