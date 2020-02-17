c
c ----------------------------------------------------------
c
       subroutine vctoprm(q,qp,mitot,mjtot,nvar)
c
c      convert to primitie variables. note for future
c      that may be overwriting conserved variables
c
       implicit double precision(a-h,o-z)
       include "cuserdt.i" 
       dimension q(nvar,mitot,mjtot),qp(nvar,mitot,mjtot)
c
       gamma1 = .4d0
       do 10 j = 1, mjtot
       do 10 i = 1, mitot
          rho       = q(1,i,j)
          u         = q(2,i,j)/rho
          v         = q(3,i,j)/rho
c         ei        = q(4,i,j)/rho - .5d0*(u**2+v**2)
c         pr        = gamma1*rho*ei 
          pr        = gamma1*(q(4,i,j)-.5d0*(q(2,i,j)**2+
     .                q(3,i,j)**2)/rho)
          qp(1,i,j) = rho
          qp(2,i,j) = u
          qp(3,i,j) = v
          qp(4,i,j) = pr
 10    continue
c
       return
       end
