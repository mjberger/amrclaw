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

c   adjusts top bndry to eliminate unstable corner
c     if ((y .ge. .1d0*x+.11d0) .and. (y .le. .1d0*x+.809d0)) then
c   adjusts bottom bndry to try to force another  unstable corner
c     if ((y .ge. .1d0*x+.0892).and.(y.le. .1d0*x+.809d0)) then

c  original channel geom  with unstable top right corner
      if ((y .ge. .1d0*x+.11d0) .and. (y .le. .1d0*x+.81d0)) then
        fbody = 1.0d0
      else
        fbody = -1.d0
      endif

      return
      end
c
