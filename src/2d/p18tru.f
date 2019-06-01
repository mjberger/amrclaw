c
c-----------------------------------------------------------------------
c
c
      subroutine p18tru(x,y,rho,u,v,p)
      implicit double precision(a-h,o-z)
      common/p18c/ delta,pmu0,theta1
c
c     # For the point (x,y), compute the true solution at this point for
c     # simple wave flow around a bend as in Whitham, p. 204
c
      pi2 = datan(1.d0)*2.d0
      delta = dsqrt(2.4d0/0.4d0)
      amach = 1.30255d0
      q0 = 1.d0
      a0 = q0/amach
      p0 = 9.04545d0
      rho0 = 1.4d0*p0/a0**2
      dmu0 = dasin(1.d0/amach)
      pmu0 = delta*datan(delta*dtan(dmu0)) - dmu0
c
c     # compute the wall location xw so that (xw,yw) and (x,y) are on
c     # the same characteristic.  Then the solution at (x,y) is equal
c     # to the solution at (xw,yw)
c
      xw = p18xw(x,y)
      theta = p18th(xw)
      dmu = p18mu(theta)
      aratio = dsqrt((1.d0+0.4d0/(2.d0*dsin(dmu0)**2)) /
     &        (1.d0+0.4d0/(2.d0*dsin(dmu)**2)))
      q = a0*aratio/dsin(dmu)
      aa = a0*aratio
      p = p0*((1.d0+0.4d0/(2.d0*dsin(dmu0)**2)) /
     &        (1.d0+0.4d0/(2.d0*dsin(dmu)**2)))**(1.4d0/0.4d0)
      rho = 1.4d0*p/aa**2
      u = q*dcos(theta)
      v = q*dsin(theta)
c     rho = 1.d0
c     u = 3.d0*x + 4.d0*y
c     v = 0.d0
c     p = 4.
      return
      end
