c
c -------------------------------------------------------------
c
c
       subroutine vprmtoc(q,mitot,mjtot,nvar)

       implicit double precision(a-h,o-z)
       dimension q(nvar,mitot,mjtot)
       include "cuserdt.i"
c
c NB: new dimensioning so can use on scrath arrays and actual soln arrays
c
       do 20 j = 1, mjtot
       do 20 i = 1, mitot
         rho  = q(1,i,j)
         u    = q(2,i,j)
         v    = q(3,i,j)
         pr   = q(4,i,j)
         ei   = pr/(rho*gamma1)
         en   = ei + .5d0*(u**2+v**2)
         q(2,i,j) = u *rho
         q(3,i,j) = v *rho
         q(4,i,j) = en*rho
 20    continue
c     
       return
       end
