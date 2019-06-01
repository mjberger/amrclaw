c
c-----------------------------------------------------------------------
c
c
c
      double precision function p18fxw(xw)
      implicit double precision (a-h,o-z)
      common/p18xwc/ x,y
c
c     # The slope of the characteristic between (x,y) and (xw,yw) is 
c     # tan(theta + mu)
c
      theta = p18th(xw)
      p18fxw = p18yw(xw)-y + dtan(theta+p18mu(theta)) *(x-xw)
      return
      end
