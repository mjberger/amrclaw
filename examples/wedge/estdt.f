c
c-----------------------------------------------------------------------
c
       subroutine estdt(val,irr,mitot,mjtot,nvar,dx,dy,dtgrid,
     1                  lwidth,aux,naux,cfl)
c
       implicit double precision (a-h, o-z)
       dimension val(nvar,mitot,mjtot), irr(mitot,mjtot)
       common /userdt/ cflcart,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                 ismp,gradThreshold
c
       dt = 1.d+20
       do 20 j = lwidth+1, mjtot-lwidth
       do 10 i = lwidth+1, mitot-lwidth
         if (irr(i,j) .ne. -1) then
             rho  = val(1,i,j)
             xmom   = dabs(val(2,i,j))
             ymom   = dabs(val(3,i,j))
             etot  = val(4,i,j)
             velsq  = (xmom**2 + ymom**2 ) /rho
             eint = etot - 0.5d0*velsq
             p    = eint*gamma1
             c    = dsqrt(gamma*p/rho)
             dt1  = dx/(c + xmom/rho)
             dt2  = dy/(c + ymom/rho)
             dt   = dmin1(dt,dt1,dt2)
         endif
 10    continue
 20    continue
c
       dtgrid = cflcart * dt
c
       return
       end
