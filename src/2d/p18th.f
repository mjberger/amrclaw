c
c-----------------------------------------------------------------------
c
c
c
      double precision function p18th(x)
      implicit double precision (a-h,o-z)
c
c     function specifying theta as a function of x along the lower wall
c     theta is the angle of the wall = angle of streamline
c
      if (x.le.0.1d0) slope = 0.d0
      if (x.gt.0.1d0 .and. x.lt.0.7d0) slope = -0.6d0*(x-0.1d0)
      if (x.ge.0.7d0) slope = -0.36d0
      p18th = datan(slope)
      return
      end
