c
c
c
c
      double precision function fbody(x,y)

      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
c
c  supersonic vortex problem
c  negative inside the body (exterior to the domain), positive otherwise.
c
      r1 = 1.d0
      r2 = 1.3840
      fbody = -1.d0

      rsq = (x**2 + y**2)
      if (rsq .ge. r1*r1 .and. rsq .le. r2*r2) fbody = 1.d0

      return
      end
c
