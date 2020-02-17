c
c -------------------------------------------------------------
c
      subroutine vrmc(qr,ql,rx,ixmin,ixmax,iflip,msize)

      implicit double precision (a-h,o-z)
      include "cuserdt.i"
      dimension qr(4,msize),ql(4,msize)
      dimension rx(4,msize)
      dimension rhor(msize),ur(msize),pr(msize),utr(msize)
      dimension rhol(msize),ul(msize),pl(msize),utl(msize)

c
c     # van Leer vector Riemann solver
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
         ur(k)   = qr(inorm,k)/qr(1,k)
         utr(k)  = qr(itan,k)/qr(1,k)
         pr(k)   = gamma1*(qr(4,k)-0.5d0*(qr(2,k)**2
     &             + qr(3,k)**2)/qr(1,k))
         rhol(k) = ql(1,k)
         ul(k)   = ql(inorm,k)/ql(1,k)
         utl(k)  = ql(itan,k)/ql(1,k)
         pl(k)   = gamma1*(ql(4,k)-0.5d0*(ql(2,k)**2
     &             + ql(3,k)**2)/ql(1,k))
  10     continue
c
      do 20 k = ixmin, ixmax

        orhor = 1.d0/rhor(k)
        unr   = ur(k)
        utanr = utr(k)

        orhol = 1.d0/rhol(k)
        unl   = ul(k)
        utanl = utl(k)

        cl2 = gamma*pl(k)*orhol
        cr2 = gamma*pr(k)*orhor
        cl = sqrt(cl2)
        cr = sqrt(cr2)
   
        odenom = .5d0/(gamma*gamma-1.d0)

        if (unl .ge. cl) then   ! completely supersonic in left state
          fl_rho     = rhol(k)*unl
          fl_thisMom = rhol(k)*unl*unl + pl(k)
          fl_tanMom  = rhol(k)*unl*utanl
          etot       = .5d0*rhol(k)*(unl*unl+utanl*utanl) + pl(k)/gamma1
          fl_rhoE    = unl*(etot + pl(k))
        else if (unl > -cl) then  ! subsonic left state (either in or out)
          f1l = gamma1*unl + 2.d0*cl
          fplus  = rhol(k) * .25d0 * (unl+cl)*(unl+cl)/cl
          fl_rho = fplus
          fl_thisMom = fplus * f1l /gamma
          fl_tanMom  = fplus * utanl
          fl_rhoE    = fplus * (utanl*utanl*.5d0 +  f1l*f1l*odenom)
        else 
          fl_rho     = 0.d0
          fl_thisMom = 0.d0
          fl_tanMom  = 0.d0
          fl_rhoE    = 0.d0
        endif

        if (unr .le. -cr) then   ! completely supersonic in right state
          fr_rho     = rhor(k)*unr
          fr_thisMom = rhor(k)*unr*unr + pr(k)
          fr_tanMom  = rhor(k)*unr*utanr
          etot       = .5d0*rhor(k)*(unr*unr+utanr*utanr) + pr(k)/gamma1
          fr_rhoE    = unr*(etot + pr(k))
        else if (unr < cr) then  ! subsonoic right state (either in or out)
          f1r = gamma1*unr - 2.d0*cr
          fminus = -rhor(k)*.25d0 * (unr-cr)*(unr-cr)/cr
          fr_rho  = fminus
          fr_thisMom = fminus*f1r /gamma
          fr_tanMom  = fminus * utanr
          fr_rhoE    = fminus * (utanr*utanr*.5d0 + f1r*f1r*odenom)
        else
          fr_rho     = 0.d0
          fr_thisMom = 0.d0
          fr_tanMom  = 0.d0
          fr_rhoE    = 0.d0
        endif


         rx(1,k)     = fl_rho + fr_rho
         rx(inorm,k) = fl_thisMom + fr_thisMom
         rx(itan,k)  = fl_tanMom  + fr_tanMom
         rx(4,k)     = fl_rhoE + fr_rhoE

 20   continue

c
       return
       end
