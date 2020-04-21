c
c
c
c
      double precision function fbody(x,y)

      implicit double precision (a-h,o-z)
c     common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
c    .                ismp,gradThreshold
c
c  negative inside the body (exterior to the domain), positive otherwise.
c

c  wedge geomtry - flat, then 30 degree ramp            
      pi = 3.14159265358979d0

c     xwall = 0.506756756756750d0
      xwall = 1.d0/6.d0
      ywall = (3.d0-xwall)*tan(30.d0*pi/180.d0)

      if (x .le. xwall) then
         if (y .gt.0.d0) then
             fbody = 1.d0
         else
             fbody = -1.d0
         endif
      else
        slope = (ywall - 0.d0)/(3.00d0-xwall)
        ybndry = (x-xwall)*slope
        if (y .gt. ybndry) then
           fbody = 1.d0
        else
           fbody = -1.d0
        endif
      endif

      return
      end
c
