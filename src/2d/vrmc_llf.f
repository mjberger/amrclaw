c
c -------------------------------------------------------------
c
      subroutine vrmc(qr,ql,rx,ixmin,ixmax,iflip,msize)

      implicit double precision (a-h,o-z)
      include "cuserdt.i"


      dimension qr(4,msize),ql(4,msize)
      dimension ur(4), ul(4)
      dimension rx(4,msize),fl(4),fr(4)
c
c     local lax friedrich with input in conserved vars
c
      if (ixmax .gt. msize) then
         write(6,*)" need to increase m1size in vrm to ", ixmax
         stop
      endif

c
c     # vector Riemann solver
c     # Given conserved variables, compute
c     # conserved variables at interface by solving riemann problem
c     # only solve neighboring vertical interfaces
c     # as in Van Leer MUSCL paper # 5 (1979), JCP
c
c      iflip = 1 indicates the normal velocity is in
c                qr(k,2) and the tangential in qr(k,3)
c            = 0 indicates the reverse.
c
      inorm = 3 - iflip
      itan  = 2 + iflip
      do 10 k = ixmin, ixmax
         rhol = ql(1,k)
         rhoul   = ql(inorm,k)
         rhoutl  = ql(itan,k)
         El   = ql(4,k)
         pl = gamma1*(El-.5d0*(rhoul*rhoul+rhoutl*rhoutl)/rhol)
         cl = dsqrt(gamma * pl / rhol)
         sl = max(dabs(rhoul/rhol + cl), dabs(rhoul/rhol - cl) )

         rhor = qr(1,k)
         rhour   = qr(inorm,k)
         rhoutr  = qr(itan,k)
         Er   = qr(4,k)
         pr = gamma1*(Er-.5d0*(rhour*rhour+rhoutr*rhoutr)/rhor)
         cr = dsqrt(gamma * pr / rhor)
         sr = max(dabs(rhour/rhor + cr), dabs(rhour/rhor - cr) )

         s = max(sl,sr)

         ul(1) = rhol
         ul(inorm) = rhoul
         ul(itan) = rhoutl
         ul(4) = El

         ur(1) = rhor
         ur(inorm) = rhour
         ur(itan) = rhoutr
         ur(4) = Er

         fl(1) = rhoul
         fl(inorm) = rhoul*rhoul/rhol + pl
         fl(itan) = rhoul*rhoutl/rhol
         fl(4) = (El + pl) * rhoul/rhol

         fr(1) = rhour
         fr(inorm) = rhour*rhour/rhor + pr
         fr(itan) = rhour*rhoutr/rhor
         fr(4) = (Er + pr) * rhour/rhor

         rx(:,k) = 0.5d0*(fl(:)+fr(:)) + 0.5d0*s*(ul(:)-ur(:))

  10     continue

      return
      end
