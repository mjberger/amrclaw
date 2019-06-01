c
c-----------------------------------------------------------------------
c
c
c
      double precision function p18yw(x)
      implicit double precision (a-h,o-z)
c
c     function specifying y as a function of x along the lower wall
c
      if (x.le.0.1d0) y = 0.302d0
      if (x.gt.0.1d0 .and. x.lt.07d0) y = 0.302d0 - 0.3d0*(x-0.1d0)**2
      if (x.ge.0.7d0) y = .194d0-0.36d0*(x-0.7d0)
      p18yw = y
      return
      end
