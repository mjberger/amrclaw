c
c
c
c
      double precision function fbody(x,y)

      use amr_module,only : xlower
      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      common /moredata/ sloc
c
c  negative inside the body (exterior to the domain), positive otherwise.
c

      sloc = 0.39d0
      !bottom = .00001d0
      bottom = .01001d0

c     no geometry - bottom of domain is (mostly) wall bndry condition
c     fbody = 1.d0
c     return

c     next lines use cut cells at bottom
      fbody = 1.d0
      if (x .gt. sloc .and. y .lt. bottom) then
         fbody = -1.d0
      endif
      if (x .le. sloc .and. y .lt. -(x-sloc)/sqrt(3.d0)+bottom) then
         fbody = -1.d0
      endif

      return
      end
c
