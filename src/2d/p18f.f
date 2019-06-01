c
c-----------------------------------------------------------------------
c
c
c
      double precision function p18f(amu)
      implicit double precision (a-h,o-z)
      common/p18c/ delta,pmu0,theta
      p18f = delta*datan(delta*dtan(amu)) - amu - pmu0 - theta
      return
      end
