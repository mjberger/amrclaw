c
c-----------------------------------------------------------------------
c
      subroutine p19fn(xcen,ycen,q)

      implicit double precision(a-h,o-z)
      dimension q(4)
c
c     # For the point (x,y), compute the true solution at this point 
c
       refrhol = 1.d0
       refmach = 2.25d0
       refc    = 1.d0

       radius = dsqrt(xcen*xcen + ycen*ycen)
       if (radius .lt. .725d0) radius = .725d0
       ratio = (1.d0/radius)**2
       theta = datan2(ycen,xcen)
       rho  =  refrhol*
     &         (1.d0+.4d0*.5d0*(refmach**2)*(1.d0-ratio))**(1.d0/.4d0)
       vel = refc*refmach*(1.d0/radius)
       u = vel*dsin(theta)
       v = -vel*dcos(theta)
       p = (1.d0/1.4d0)*rho**1.4d0

       q(1) = rho
       q(2) = u
       q(3) = v
       q(4) = p

       return
       end
