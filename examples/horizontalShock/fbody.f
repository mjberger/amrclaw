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

c no geometry
c     fbody = 1.d0
c     return


c  flat geometry - for easy debugging

      if (x .le. 4.0d0 .and. y .ge. .00001d0) then
             fbody = 1.d0
         else
             fbody = -1.d0
         endif
      return
      end
c
