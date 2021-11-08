c
c-----------------------------------------------------------------------
c
      subroutine ssvInit(xcen,ycen,rho,u,v,p)
      implicit double precision(a-h,o-z)
      logical real/.true./
c
c     # For the point (x,y), compute the true solution at this point 
c
       if (real) then
       refrho   = 1.d0
       refmach  = 2.25d0
       refc     = 1.d0

       radius = dsqrt(xcen*xcen + ycen*ycen)
       if (radius .lt. 0.725d0) radius = .725d0
       ratio = (1.d0/radius)**2
       theta = datan2(ycen,xcen)
       rho  =  refrho*
     &         (1.d0+.4d0*.5d0*(refmach**2)*(1.-ratio))**(1.d0/.4d0)
       vel = refc*refmach*(1.d0/radius)
       u = vel*dsin(theta)
       v = -vel*dcos(theta)
       p = (1.d0/1.4d0)*rho**1.4d0

       else
       !linear test case for debugging
       rho = xcen + 2.d0*ycen + 5.d0
       !rho = xcen  + 5.d0
       !rho =  5.d0
       u = 0.d0
       v = 0.d0
       p = 10.d0 
       endif

       return
       end
