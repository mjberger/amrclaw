c
c ---------------------------------------------------------------
c
      subroutine weights(p,xc,yc,totxx,totxy,totyy,wt,area,points,
     .                   hx,hy,k)
      implicit double precision (a-h,o-z)
      dimension  p(10,2), wt(24)
      dimension  subar(20), points(24,2)
      dimension  txx(20),txy(20),tyy(20)
c
c     # compute integral over cell p of quadratic terms
c     # shifted by centroid.  have already computed 
c     # centroid (xc,yc) of the polygon p.  
c
      index = 1
      do 3 i = 1, 24
        wt(i) = 0.d0
 3    continue

      x0 = xc
      y0 = yc
      do 100 i=2,20
c        if (p(i+1,1) .eq. -1) then
         if (p(i,1) .eq. -11) then
           if (i.eq.2) then
             write(6,*) 'degenerate polygon in weights'
             stop
           endif
         go to 200
         endif

c        # Compute centroid of triangle with vertices 1, i-1, and i
c
         x1 = p(i-1,1)
         y1 = p(i-1,2)
         x2 = p(i,1)
         y2 = p(i,2)

c
c        # compute area:
c 
         subar(i) = 0.5d0 * ((y0+y1)*(x1-x0) + (y1+y2)*(x2-x1)
     &                     + (y2+y0)*(x0-x2))
c
c        # use 3rd order approximation of average of midpoints
c
         xmid = (x0+x1)/2.d0
         ymid = (y0+y1)/2.d0
         points(index,1)  = xmid
         points(index,2)  = ymid
         wt(index) = wt(index) + 1.d0/3.d0*subar(i)/area
         fnxx1 = (xmid-xc)*(xmid-xc)
         fnyy1 = (ymid-yc)*(ymid-yc)
         fnxy1 = (xmid-xc)*(ymid-yc)

         xmid = (x2+x1)/2.d0
         ymid = (y2+y1)/2.d0
         index = index + 1
         points(index,1)  = xmid
         points(index,2)  = ymid
         wt(index) = wt(index) + 1.d0/3.d0*subar(i)/area
         fnxx2 = (xmid-xc)*(xmid-xc)
         fnyy2 = (ymid-yc)*(ymid-yc)
         fnxy2 = (xmid-xc)*(ymid-yc)

         xmid = (x0+x2)/2.d0
         ymid = (y0+y2)/2.d0
         index = index + 1
         points(index,1)  = xmid
         points(index,2)  = ymid
         wt(index) = wt(index) + 1.d0/3.d0*subar(i)/area
         fnxx3 = (xmid-xc)*(xmid-xc)
         fnyy3 = (ymid-yc)*(ymid-yc)
         fnxy3 = (xmid-xc)*(ymid-yc)

         txx(i-1)  = .5d0*subar(i)*(fnxx1 + fnxx2 + fnxx3) / 3.d0
         txy(i-1)  =      subar(i)*(fnxy1 + fnxy2 + fnxy3) / 3.d0
         tyy(i-1)  = .5d0*subar(i)*(fnyy1 + fnyy2 + fnyy3) / 3.d0

c
  100    continue
      write(6,*) 'error in weights:  too many vertices'
      stop
c
  200 continue
c
c
c   compute final integrals of quadratic terms
      totxx   = 0.d0
      totxy   = 0.d0
      totyy   = 0.d0
      checkar = 0.d0
      do  230 j = 2, i-1
         totxx = totxx + txx(j-1)
         totxy = totxy + txy(j-1)
         totyy = totyy + tyy(j-1)
         checkar = checkar + subar(j)
 230   continue
      totxx = totxx / area
      totxy = totxy / area
      totyy = totyy / area
      if (dabs(checkar-area)/area .gt. 1.e-9) then
          write(6,*) 'areas differ: ', checkar,area
      endif
c     
      return
      end
