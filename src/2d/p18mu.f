c
c-----------------------------------------------------------------------
c
c     
c
      double precision function p18mu(theta)
      implicit double precision (a-h,o-z)
      external p18f
      common/p18c/ delta,pmu0,theta1
c
c     compute mu as a function of theta by solving eqn (6.168) in
c     Whitham.  Theta must be passed to f in common block.
c
      theta1 = theta
      p18mu = zeroin(0., 1.57d0, p18f, 1d-10)
      return
      end
