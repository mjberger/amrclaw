c
c-----------------------------------------------------------------------
c
c
      double precision function p18xw(x,y)
      implicit double precision (a-h,o-z)
      external p18fxw
      common/p18xwc/ x1,y1
c
c     # for arbitrary point (x,y), find the point xw so that the point
c     # (xw,yw) and the point (x,y) lie on the same charactersitic.
c
      x1 = x
      y1 = y
      p18xw = zero2(-1.d0, 2.d0, p18fxw, 1d-9)
      return
      end
