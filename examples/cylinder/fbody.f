c
c
c
c
      double precision function fbody(x,y)

      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
c
c  negative inside the body (exterior to the domain), positive otherwise.
c

c  cylinder geometry
      pi = 3.14159265358979d0
      ptDist2 = (x-2.d0)**2 + (y-2.d0)**2
      radius = 0.5d0
      radiusSquared = radius**2
      

      if (ptDist2 .le. radiusSquared) then
           fbody = -1.d0
      else
           fbody =  1.d0
      endif

      return
      end
c
