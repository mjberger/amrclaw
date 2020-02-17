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
      ggm1 = gamma/(gamma-1.d0)

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

c  get Roe average speeds
        dd = rhol(k)/rhor(k)
        odd = 1.d0/(dd+1.d0)
        roeAvg_normvel = ur(k) + dd*ul(k)*odd
        roeAvg_tanvel  = utr(k) + dd*utl(k)*odd
        hL = (pl(k)*ggm1)/rhol(k) + 0.5d0*(ul(k)**2+utl(k)**2)
        hR = (pR(k)*ggm1)/rhor(k) + 0.5d0*(ur(k)**2+utr(k)**2)
        avg_H = (hR+dd*hL)*odd
        q2 = roeAvg_normvel**2 + roeAvg_tanvel**2
        roeAvg_ssp = sqrt(max(gamma1*(avg_H - 0.5d0*q2),1d-6))

        sL = min((ul(k) - aL), (roeAvg_normvel - roeAvg_ssp)) ! fancier option  
        sR = max((ur(k) + aL), (roeAvg_normvel + roeAvg_ssp))
        sM = (rhor(k)*ur(k) * (sR-ur(k)) - rhol(k)*ul(k)*(sL-ul(k)) +
     .       pl(k)-pr(k)) / (rhor(k)*(sR-ur(k))-rhol(k)*(sL-ul(k)))

c       if (.not.((sL .le. sM) .and. (sM.le.sR))) then
c           write(*,*)"sM outside roe avg speeds"
c       endif

        sL = min((ul(k) - aL), (uR(k)-aR))
        sR = max((ul(k) + aL), (uR(k)+aR))
        sM = (rhor(k)*ur(k) * (sR-ur(k)) - rhol(k)*ul(k)*(sL-ul(k)) +
     .        pl(k)-pr(k)) / (rhor(k)*(sR-ur(k))-rhol(k)*(sL-ul(k)))

       if (.not.((sL .le. sM) .and. (sM.le.sR))) then
            write(*,*)"sM still outside more diffusive wave speeds"
        endif

        if (0.d0 .lt. sL) then   ! use left flux

         rx(1,k)     = rhol(k)*ul(k)
         rx(inorm,k) = rhol(k)*ul(k)*ul(k) + pl(k)
         rx(itan,k)  = rhol(k)*ul(k)*utl(k)
         etot        = .5d0*rhol(k)*(ul(k)*ul(k)+utl(k)*utl(k)) + 
     .                  pl(k)/gamma1
         rx(4,k)     = ul(k)*(etot + pl(k))

        else if (0. .gt. sR) then   

         rx(1,k)     = rhor(k)*ur(k)
         rx(inorm,k) = rhor(k)*ur(k)*ur(k) + pr(k)
         rx(itan,k)  = rhor(k)*ur(k)*utr(k)
         etot        = .5d0*rhor(k)*(ur(k)*ur(k)+utr(k)*utr(k)) + 
     .                  pr(k)/gamma1
         rx(4,k)     = ur(k)*(etot + pr(k))

        else 

         pStar  = rhol(k)*(ul(k)-sL)*(ul(k)-sM) + pl(k) 
         

         if (sM .ge. 0.d0) then
            rlx(1)     = rhol(k)*ul(k)
            rlx(inorm) = rhol(k)*ul(k)*ul(k) + pl(k)
            rlx(itan)  = rhol(k)*ul(k)*utl(k)
            etotl        = .5d0*rhol(k)*(ul(k)*ul(k)+utl(k)*utl(k)) +
     .           pl(k)/gamma1
            rlx(4)     = ul(k)*(etotl + pl(k))

            omegaL = 1.d0/(sL-sM)
            rhoStar = omegaL*rhol(k)*(sL-ul(k))
            rhoUnStar = omegaL*(rhol(k)*(sL-ul(k))*ul(k) + pStar-pl(k))
            rhotanStar = omegaL*rhol(k)*(sL-ul(k))*utl(k)
            eStar = omegaL*((sL-ul(k))*etotL - pl(k)*ul(k) + pStar*sM)

            rx(1,k)     = rlx(1) + sL*(rhoStar - rhol(k))
            rx(inorm,k) = rlx(inorm) + sL*(rhoUnstar - rhol(k)*ul(k))
            rx(itan,k)  = rlx(itan)  + sL*(rhotanStar - rhol(k)*utl(k))
            rx(4,k)       = rlx(4) + sL*(eStar - etotL)

         else

            rrx(1)     = rhor(k)*ur(k)
            rrx(inorm) = rhor(k)*ur(k)*ur(k) + pr(k)
            rrx(itan)  = rhor(k)*ur(k)*utr(k)
            etotR        = .5d0*rhor(k)*(ur(k)*ur(k)+utr(k)*utr(k)) +
     .           pr(k)/gamma1
            rrx(4)     = ur(k)*(etotr + pr(k))

            omegaR = 1.d0/(sR-sM)
            rhoStar = omegaR*rhor(k)*(sR-ur(k))
            rhoUnStar = omegaR*(rhor(k)*(sR-ur(k))*ur(k) + pStar-pr(k))
            rhotanStar = omegaR*rhor(k)*(sR-ur(k))*utr(k)
            eStar = omegaR*((sR-ur(k))*etotR - pr(k)*uR(k) + pStar*sM)

            rx(1,k)     = rrx(1)     + sR*(rhoStar - rhor(k))
            rx(inorm,k) = rrx(inorm) + sR*(rhoUnstar - rhor(k)*ur(k))
            rx(itan,k)  = rrx(itan)  + sR*(rhotanStar - rhor(k)*utr(k))
            rx(4,k)       = rrx(4)     + sR*(eStar - etotR)

          endif
         endif

 20   continue

c
      return
      end
