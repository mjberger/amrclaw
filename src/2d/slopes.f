c
c------------------------------------------------------------
c
       subroutine slopes(qp,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,
     .                   hx,hy,xlow,ylow,nvar,mptr,istage)

       use amr_module
       implicit double precision(a-h,o-z)

       dimension qp(nvar,mitot,mjtot),qx(nvar,mitot,mjtot),
     &           qy(nvar,mitot,mjtot),irr(mitot,mjtot)
       logical  regular, quad, nolimiter, pwconst
      common   /userdt/ cflcart,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold,pwconst
       common /order2/ ssw, quad, nolimiter

c      regular(i,j) = ((i.gt. lwidth).and.(i.le.mitot-lwidth).and.
c    &                 (j.gt. lwidth).and.(j.le.mjtot-lwidth))

c
c      # ssw = slope switch (1. for slopes, 0 for donor cell 0 slopes)
c      # now set in amrcart
c
c   ssw = 0 for ctu (no slopes), 1 for muscl (second order).
c   compute slopes using muscl limiter
c   qp contains the primitive variables
c

c
       if (ssw .eq. 0.d0) go to 99
c
       do 10 j = 1, mjtot
       do 10 i = 2, mitot-1
       do 11 m = 1, nvar
          if (irr(i,j) .ne. lstgrd) go to 10
          ducc = qp(m,i+1,j) - qp(m,i-1,j)
          dupc = qp(m,i+1,j) - qp(m,i,j)
          dumc = qp(m,i,j)   - qp(m,i-1,j)

c         one sided at domain boundaries
c         x = xlow + (dfloat(i)-.5d0)*hx          ! cell center of this cell
c         if (x+hx .gt. xprob) ducc = 2.d0*dumc ! last pt has rt nbor outside domain
c         if (x-hx .lt. 0.d0)  ducc = 2.d0*dupc
c         turn off limiting if executing next 2 lines, with go to
          if (nolimiter) then
             qx(m,i,j) = .5d0*ducc/hx
c                # adjust for irregularity, if not done in qslopes
c$$$             if (irr(i+1,j) .ne. lstgrd) then
c$$$                 qx(m,i,j) = dumc/hx
c$$$             else if (irr(i-1,j) .ne. lstgrd) then
c$$$                 qx(m,i,j) = dupc/hx
c$$$             endif
             go to 11
          endif
c
         du   = dmin1(dabs(dupc),dabs(dumc))
         du   = dmin1(2.d0*du, .5d0*dabs(ducc))
c
         fl = dmax1(0.d0,dsign(1.d0, dupc*dumc))*ssw
         qx(m,i,j) = du*dsign(1.d0,ducc)*fl/hx
 11    continue
 10    continue
c
       do 20 j = 2, mjtot-1
       do 20 i = 1, mitot
       do 20 m = 1, nvar
       if (irr(i,j) .ne. lstgrd) go to 20
          ducc = qp(m,i,j+1) - qp(m,i,j-1)
          dupc = qp(m,i,j+1) - qp(m,i,j)
          dumc = qp(m,i,j)   - qp(m,i,j-1)
c         turn off limiting if execute next 2 lines, with go to
          if (nolimiter) then
             qy(m,i,j) = .5d0*ducc/hy
c                # adjust for irregularity, if not done in qslopes
c$$$             if (irr(i,j+1) .ne. lstgrd) then
c$$$                qy(m,i,j) = dumc/hy
c$$$             else if (irr(i,j-1) .ne. lstgrd) then
c$$$                qy(m,i,j) = dupc/hy
c$$$             endif
             go to 20
          endif
c 
          du   = dmin1(dabs(dupc),dabs(dumc))
          du   = dmin1(2.d0*du, .5d0*dabs(ducc))
c
          fl = dmax1(0.d0,dsign(1.d0, dupc*dumc))*ssw
          qy(m,i,j) = du*dsign(1.d0,ducc)*fl/hy
 20    continue
c
 99    continue
       return
       end
