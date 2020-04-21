c
c     =================================================
      function fdisc(x,y,time)
c     =================================================
c
c     # for computing cell averages for initial data that has a
c     # discontinuity along some curve.  fdisc should be negative to the 
c     # left of the curve and positive to the right
c     # idisc specifies the nature of the discontinuity for two
c     # particular cases (a straight line and circle) but this routine
c     # can be modified for any other curve.

c     thanks for Randy LeVeque for his cellave routine.
c
c
      use amr_module, only:  ylower
      implicit double precision (a-h,o-z)
      include "cdisc.i"
      common /moredata/ sloc

      idisc = 10  ! for straight line shock problem
      go to (10,20) idisc
c
   10 continue
c     # straight line through (x0,y0) with normal (alf,beta) pointing 
c     # into right state
c
      ! mach 10 30 degree ramp 
      pi = 4.d0*atan(1.0d0)
      alf = cos(pi/10.d0)
      beta = -sin(pi/10.d0)
      y0 = ylower
      x0 = sloc + 20.d0*time/sqrt(3.d0)

      fdisc = (x-x0)*alf + (y-y0)*beta
      return
c
   20 continue
c     # circle of radius r0:
      fdisc = (x-x0)**2 + (y-y0)**2 - r0**2
c
      return
      end
