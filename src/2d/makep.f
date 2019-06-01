c
c
c ---------------------------------------------------------------
c
      subroutine makep(poly,i,j,xlow,ylow,dx,dy)

      implicit double precision (a-h,o-z)
      dimension poly(10,2)
c
c     # create a list of vertices of the polygon corresponding to grid
c     # cell (i,j).
c
      poly(1,1) = xlow + (i-1)*dx
      poly(1,2) = ylow + (j-1)*dy
      poly(2,1) = xlow + (i-1)*dx
      poly(2,2) = ylow + j*dy
      poly(3,1) = xlow + i*dx
      poly(3,2) = ylow + j*dy
      poly(4,1) = xlow + i*dx
      poly(4,2) = ylow + (j-1)*dy
      poly(5,1) = xlow + (i-1)*dx
      poly(5,2) = ylow + (j-1)*dy
c     poly(6,1) = -1
      poly(6,1) = -11
      return
      end
