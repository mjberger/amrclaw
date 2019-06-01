
c
c
c------------------------------------------------------------------------
c
c
      function fs(s)
      implicit double precision (a-h,o-z)
      common/fscorn/ xc0,yc0,xc1,yc1
c   
c     # compute fbody at distance s between corners (xc0,yc0) and (xc1,yc1)
c
      x = xc0 + s*(xc1-xc0)
      y = yc0 + s*(yc1-yc0)
      fs = fbody(x,y)
      return
      end
