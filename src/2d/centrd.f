c
c ---------------------------------------------------------------
c
      subroutine centrd(p,xc,yc,area)
      implicit double precision (a-h,o-z)
      dimension p(10,2)
      dimension xcen(10),ycen(10),subar(10)
c
c     # Compute the centroid (xc,yc) of the polygon p.  Also returns
c     # the area.
c
      epsar = 1d-7
      x0 = p(1,1)
      y0 = p(1,2)
      do 100 i=3,10
         if (p(i+1,1) .eq. -11) then
           if (i.eq.3) then
             write(6,*) 'degenerate polygon in centrd'
             stop
             endif
           go to 200
           endif
c
c        # Compute centroid of triangle with vertices 1, i-1, and i
c
         x1 = p(i-1,1)
         y1 = p(i-1,2)
         x2 = p(i,1)
         y2 = p(i,2)
         xcen(i) = (x1+x2+x0)/3.d0
         ycen(i) = (y1+y2+y0)/3.d0
c
c        # compute area:
c 
         subar(i) = 0.5d0 * ((y0+y1)*(x1-x0) + (y1+y2)*(x2-x1)
     &                     + (y2+y0)*(x0-x2))
c
  100    continue
      write(6,*) 'error in centrd:  too many vertices'
      stop
c
  200 continue
c
c     # take weighted average of triangle centroids to find centroid
c     # of polygon:
c
      area = 0.d0
      xc = 0.d0
      yc = 0.d0
      do 220 j=3,i-1
         area = area + subar(j)
         xc = xc + subar(j)*xcen(j)
         yc = yc + subar(j)*ycen(j)
  220    continue
c     if (area.eq.0.d0) then
      if (area.eq.1.d-16) then
          xc = x0
          yc = y0
        else
          xc = xc/area
          yc = yc/area
        endif
c     
      return
      end
