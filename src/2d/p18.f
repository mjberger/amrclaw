c
c-----------------------------------------------------------------------
c
c
      subroutine p18(xupw,yupw,n)
      implicit double precision(a-h,o-z)
      dimension xupw(n),yupw(n)
      common/p18c/ delta,pmu0,theta1
c
c     # find the upper wall of the rarefaction channel for problem 18
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
      h = 1.7d0/dfloat(n)
      xw0 = .1d0
      yw0 = .302d0
      width0 = .2d0
      theta0 = 0.d0
      x0 = xw0 + width0/dtan(theta0+dmu0)
      y0 = yw0 + width0
      xupw(1) = 0.d0
      yupw(1) = y0
      xupw(2) = x0
      yupw(2) = y0
c
      do 10 i=1,n-2
	 xw = xw0 + i*h
	 yw = p18yw(xw)
	 theta = p18th(xw)
	 dmu = p18mu(theta)
	 aratio = dsqrt((1.d0+0.4d0/(2.d0*dsin(dmu0)**2)) /
     &        (1.d0+0.4d0/(2.d0*dsin(dmu)**2)))
	 q = a0*aratio/dsin(dmu)
	 aa = a0*aratio
	 p = p0*((1.d0+0.4d0/(2.d0*dsin(dmu0)**2)) /
     &        (1.d0+0.4d0/(2.d0*dsin(dmu)**2)))**(1.4d0/0.4d0)
	 rho = 1.4d0*p/aa**2
	 width = rho0*width0*q0 / (rho*q*dcos(pi2-dmu))
	 x = xw + width*dcos(theta+dmu)
	 y = yw + width*dsin(theta+dmu)
	 xupw(i+2) = x
	 yupw(i+2) = y
  10     continue
      return
      end
