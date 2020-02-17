c
c -------------------------------------------------------------
c
      subroutine vrm(qr,ql,rx,ixmin,ixmax,iflip,msize)
      implicit double precision (a-h,o-z)
      include "cuserdt.i"
      dimension qr(4,msize),ql(4,msize)
      dimension rx(4,msize), rlx(4), rrx(4), stateR(4), stateL(4)
      dimension rhor(msize),ur(msize),pr(msize),utr(msize)
      dimension rhol(msize),ul(msize),pl(msize),utl(msize)

      data      itno/2/, small/1.d-6/
c
c
c     # HLL vector Riemann solver
c     # Given primitive variables, compute
c     # conserved variables at interface by solving riemann problem
c     # only solve neighboring vertical interfaces
c
c      iflip = 1 indicates the normal velocity is in
c                qr(k,2) and the tangential in qr(k,3)
c            = 0 indicates the reverse.
c   done so a single riemann solver can be used for both columns and rows.
c
      if (ixmax .gt. msize) then
         write(6,*)" need to increase msize in vrm to ", ixmax
         stop
      endif

      inorm = 3 - iflip
      itan  = 2 + iflip

      do 10 k = ixmin, ixmax
         rhor(k) = qr(1,k)
         ur(k)   = qr(inorm,k)
         utr(k)  = qr(itan,k)
         pr(k)   = qr(4,k)
         rhol(k) = ql(1,k)
         ul(k)   = ql(inorm,k)
         utl(k)  = ql(itan,k)
         pl(k)   = ql(4,k)
  10     continue
c
      do 20 k = ixmin, ixmax
        aL  = dsqrt(gamma*pl(k)/rhol(k))
        aR  = dsqrt(gamma*pr(k)/rhor(k))

c       sL = ul(k) - aL; 
c       sR = ur(k) + aR; 
        sL = min((ul(k) - aL), (ur(k) - aR)) ! fancier option  
        sR = max((ul(k) + aL), (ur(k) + aR))

        if (0.d0 .lt. sL) then   ! use left flux

         rx(1,k)     = rhol(k)*ul(k)
         rx(inorm,k) = rhol(k)*ul(k)*ul(k) + pl(k)
         rx(itan,k)  = rhol(k)*ul(k)*utl(k)
         etot        = .5d0*rhol(k)*(ul(k)*ul(k)+utl(k)*utl(k)) + 
     .                  pl(k)/gamma1
         rx(4,k)     = ul(k)*(etot + pl(k))

        else if (0. .gt. sR) then   ! flux is right flux

         rx(1,k)     = rhor(k)*ur(k)
         rx(inorm,k) = rhor(k)*ur(k)*ur(k) + pr(k)
         rx(itan,k)  = rhor(k)*ur(k)*utr(k)
         etot        = .5d0*rhor(k)*(ur(k)*ur(k)+utr(k)*utr(k)) + 
     .                  pr(k)/gamma1
         rx(4,k)     = ur(k)*(etot + pr(k))

        else
        ! middle flux
         rlx(1)     = rhol(k)*ul(k)
         rlx(inorm) = rhol(k)*ul(k)*ul(k) + pl(k)
         rlx(itan)  = rhol(k)*ul(k)*utl(k)
         etotl        = .5d0*rhol(k)*(ul(k)*ul(k)+utl(k)*utl(k)) +
     .                  pl(k)/gamma1
         rlx(4)     = ul(k)*(etotl + pl(k))

         rrx(1)     = rhor(k)*ur(k)
         rrx(inorm) = rhor(k)*ur(k)*ur(k) + pr(k)
         rrx(itan)  = rhor(k)*ur(k)*utr(k)
         etotr        = .5d0*rhor(k)*(ur(k)*ur(k)+utr(k)*utr(k)) +
     .                  pr(k)/gamma1
         rrx(4)     = ur(k)*(etotr + pr(k))

         stateL(1)     = rhol(k)
         stateL(inorm) = rhol(k)*ul(k)
         stateL(itan)  = rhol(k)*utl(k)
         stateL(4)     = etotl

         stateR(1)     = rhor(k)
         stateR(inorm) = rhor(k)*ur(k)
         stateR(itan)  = rhor(k)*utr(k)
         stateR(4)     = etotr

         rx(:,k)=(sR*rlx(:)-sL*rrx(:)
     .                       +sL*sR*( stateR(:)-stateL(:) ) )/ (sR-sL)
        endif
 20   continue

c
       return
       end
