c
c
c
c
      double precision function fbody(x,y)

      use amr_module,only : xlower
      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
c
c  negative inside the body (exterior to the domain), positive otherwise.
c

      fbody = 1.d0
      return

c     not quite working with cut cells due to lower left boundary
      fbody = 1.d0
      if ((x .gt. xlower .and. y .lt. .00001d0) .or.
     .    (x .lt. xlower .and. y .lt. (x-xlower)/sqrt(3.d0))) then
	 fbody = -1.d0
      endif

      return
      end
c
