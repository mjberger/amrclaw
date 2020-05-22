c
c -------------------------------------------------------------
c
      subroutine vrm(qr,ql,rx,ixmin,ixmax,iflip,msize)

      implicit double precision (a-h,o-z)
      include "cuserdt.i"


      parameter ( psmall = 1d-07 , rsmall = 1d-07)

      dimension qr(4,msize),ql(4,msize)
      dimension rx(4,msize)
      dimension rhor(msize),ur(msize),pr(msize),utr(msize)
      dimension rhol(msize),ul(msize),pl(msize),utl(msize)
      dimension clsqp(msize),clsqm(msize),pstar(msize)
      dimension zp(msize),zm(msize),ustarm(msize),ustarp(msize)
      dimension zsum(msize)
      data      itno/2/, small/1.d-6/
c
      if (ixmax .gt. msize) then
         write(6,*)" need to increase msize in vrm to ", ixmax
         stop
      endif

      gp1g2i = (gamma+1.d0)/(2.d0*gamma)
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
c     # iterate to find pstar, ustar:
c
      do 20 k = ixmin, ixmax
         clsqp(k) = gamma*rhol(k)*pl(k)
         clsqm(k) = gamma*rhor(k)*pr(k)
         wp = dsqrt(clsqp(k))
         wm = dsqrt(clsqm(k))
         pstar(k) = (wp*pr(k)+ wm*pl(k)- wm*wp
     &      *(ur(k)- ul(k)))/(wp+ wm)
         pstar(k) = dmax1(pstar(k),small)
  20     continue
c
      do 30 iter = 1, itno
         do 30 k = ixmin, ixmax
           wsqp = clsqp(k)*(1.d0+gp1g2i*(pstar(k)/pl(k)-1.d0))
           wp   = dsqrt(wsqp)
           wsqm = clsqm(k)*(1.d0+gp1g2i*(pstar(k)/pr(k)-1.d0))
           wm   = dsqrt(wsqm)
c 
           zp(k) = 2.d0*wp*wsqp/(wsqp+clsqp(k))
           zm(k) = 2.d0*wm*wsqm/(wsqm+clsqm(k))
c 
           ustarm(k) = ur(k)-(pr(k)-pstar(k))/wm
           ustarp(k) = ul(k)+(pl(k)-pstar(k))/wp
           zsum(k) = zp(k)+zm(k)
           epsp = zp(k)*zm(k)*(ustarm(k)-ustarp(k))/zsum(k)
           pstar(k) = pstar(k)-epsp
           pstar(k) = dmax1(pstar(k),small)
  30       continue
c
      do 40 k = ixmin, ixmax
         ustar = (zp(k)*ustarm(k)+zm(k)*ustarp(k))/zsum(k)
         sein  = -dsign(1.d0,ustar)
         if (sein .ge. 0.d0) then
            uo = ur(k)
            po = pr(k)
            ro = rhor(k)
            utgdnv = utr(k)
         else
            uo = ul(k)
            po = pl(k)
            ro = rhol(k)
            utgdnv = utl(k)
         endif
         co = dsqrt(po*gamma/ro)

         dummy = pstar(k)-po
         wo = gamma*ro*po*(1.d0+gp1g2i*dummy/po)
         rstar = ro/(1.d0-ro*dummy/wo)
         wo = dsqrt(wo)
         cstar = dsqrt(gamma*pstar(k)/rstar)
         wso = sein*uo+co
         wsi = sein*ustar+cstar
         ushok = sein*ustar+wo/rstar
         if (dummy .ge. 0.d0) then
            wsi = ushok
            wso = ushok
         else
            wsi = wsi
            wso = wso
         endif
         t4 = (wso+wsi)/dmax1(wso-wsi,wso+wsi,small)
         t4 = .5d0*(t4+1.d0)
         t5  = 1.d0-t4
         t3  = t4*cstar+t5*co
         rgdnv = t4*rstar+t5*ro
         ugdnv = t4*ustar+t5*uo
         pgdnv = t3*t3*rgdnv/gamma
         
         if (wso .ge. 0.d0) then
            pgdnv = pgdnv
            rgdnv = rgdnv
            ugdnv = ugdnv
         else
            pgdnv = po
            rgdnv = ro
            ugdnv = uo
         endif
         
         if (wsi .lt. 0.d0) then
            pgdnv = pgdnv
            rgdnv = rgdnv
            ugdnv = ugdnv
         else
            pgdnv = pstar(k)
            rgdnv = rstar
            ugdnv = ustar
         endif

         rx(1,k)     = rgdnv*ugdnv
         rx(inorm,k) = rx(1,k)*ugdnv + pgdnv
         rx(itan,k)  = rx(1,k)*utgdnv
         rx(4,k)     = rx(1,k)*(pgdnv*gamma/(rgdnv*gamma1)+
     &                 .5d0*(ugdnv**2 + utgdnv**2))
 40   continue
c     
      return
      end
