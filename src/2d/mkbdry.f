c
c--------------------------------------------------------------------
c
c
      subroutine mkbdry(k0,istart,jstart,irr,mitot,mjtot,xlow,ylow,
     &           hx,hy,lstgrd,mptr,gradThreshold)

      use amr_module
      implicit double precision (a-h,o-z)

      external fs
      logical   lfix, lsmall, debug
      dimension x(7222),y(7222)
      integer irr(mitot,mjtot), idx(8), idy(8)
      common /sxy/ s1,x0,y0,iedge
      common/fscorn/ xc0,yc0,xc1,yc1
      character*1 charPrint
      data  debug/.true./
      data idx /0,1,0,-1,0,1,0,-1/
      data idy /1,0,-1,0,1,0,-1,0/
c
c     # make boundary segment.  Sets up arrays irr, poly, ar, ix,
c     #    iy and nxtirr for irregular boundary cells on one segment.
c     # On input, 
c     #   k0     = pointer to last irregular grid location already used
c     #   (istart,jstart) = coordinates of grid cell in which boundary
c     #                     segment begins.  
c     #   mitot  = total number of grid cells in x-direction.
c     #   mjtot  = total number of grid cells in y-direction.
c     #   (xlow,ylow)  = coordinates of lower left corner of grid
c     #   dx,dy  = mesh spacing
c     #
c     # some local variables:
c     #   iedge   = edge of cell currently being examined. 
c     #               1 = left, 2 = top, 3 = right, 4 = bottom, 5 = left, etc.
c     #             repeats periodically for convenience.
c     #   (x0,y0) = first place boundary itersects this cell (as we move
c     #             around boundary with wall on right).
c     #   (x1,y1) = second place boundary itersects this cell 
c     #   (ixc*,iyc*) (for *=0 or 1) = index of a corner
c     #   (xc*,yc*) = coordinates of that corner
c     #   idx(j)  = change in ixc index when traversing edge j
c     #   idy(j)  = change in iyc index when traversing edge j
c     #
c     # basic approach:
c     #   Starting at (x0,y0) for one cell, use zeroin to find (x1,y1),
c     #   set up polygon information for this cell, and then (x1,y1) becomes
c     #   (x0,y0) for the adjacent cell.
c     # 
c
      if (debug) write(13,900) mptr
 900  format(' making boundary for grid ',i4)
      if (mitot .gt. 7222 .or. mjtot .gt. 7222) then
        write(6,*)" make mkbdry scratch arrays bigger "
        stop
      endif
      dx = hx
      dy = hy
      eps = 1d-4*dmax1(dx,dy)
c     tol = 1d-17
      tol = 1d-19
c     epside = 1.d-10
      epside = 1.d-16
      ixuse = istart
      iyuse = jstart
      ix0 = istart
      iy0 = jstart
      xup = xlow + mitot*dx
      yup = ylow + mjtot*dy
      do 5 i=1,mitot+1
         x(i) = xlow + dfloat(i-1)*dx
 5    continue
      do 6 j=1,mjtot+1
         y(j) = ylow + dfloat(j-1)*dy
 6    continue
      k = k0
c     
c     # find side of cell (istart,jstart) on which boundary starts
      ixc0 = ixuse
      iyc0 = iyuse
      f0 = fbody(x(ixc0),y(iyc0))
      do 10 iedge=1,4
         ixc1 = ixc0 + idx(iedge)
         iyc1 = iyc0 + idy(iedge)
         f1 = fbody(x(ixc1),y(iyc1))
         if (f0.lt.0. .and. f1.gt.0.) go to 20
         if (f1.eq.0.) write(6,601)
         f0 = f1
         ixc0 = ixc1
         iyc0 = iyc1
 10   continue
      write(6,*) 'error *** fbody does not change sign on boundary of'
      write(6,*) '          cell (istart,jstart) in mkbdry'
      stop
 20   continue
c     
c     # compute point on this side where boundary starts. 
c     
      xc0 = x(ixc0)
      yc0 = y(iyc0)
      xc1 = x(ixc1)
      yc1 = y(iyc1)
      z = zeroin(0.d0,1.d0,fs,tol)
      x0 = xc0 + z*(xc1-xc0)
      y0 = yc0 + z*(yc1-yc0)
c     
c     # beginning of main loop to march around boundary:
c     
      lfix = .false.
      lastk = 0
      do 200 ibcell=1,mitot*mjtot
c     
c     # get index k for storage of info for this irregular cell, unless
c     # lfix=true which indicates that previous cell was eliminated
c     # and hence the previous value of k can be used again:
         if (.not. lfix) then
            nextk = lstget(dummy)
            if (k .eq. k0) then
               nxtirr(k) = -nextk
            else
               nxtirr(k) = nextk
            endif
            k = nextk
         endif
c     
c       # compute coordinates of next corner past (x0,y0), going
c       # clockwise around the boundary.
         ixc0 = ixuse
         iyc0 = iyuse
         do 30 i=1,iedge
            ixc0 = ixc0 + idx(i)
            iyc0 = iyc0 + idy(i)
 30      continue
         xc0 = x(ixc0)
         yc0 = y(iyc0)
         f0 = fbody(xc0,yc0)
         if (f0.lt.0.) then
            if (lfix) then
               write(13,*) '> skipping cell ',ixuse,iyuse
               irr(ixuse,iyuse) = -1
               go to 110
            else
               write(6,*) 'error in mkbdry *** negative f0'
               stop
            endif
         endif
         if (f0.eq.0.) write(6,601)
c     
c     # start constructing polygon for this cell, starting with the
c     # current side, which goes from (x0,y0) to (xc0,yc0):
c     
         ixg(k) = ixuse
         iyg(k) = iyuse
         irr(ixuse,iyuse) = k
c     
         poly(1,1,k) = x0
         poly(1,2,k) = y0
         poly(2,1,k) = xc0
         poly(2,2,k) = yc0
c     
c     # march around the cell clockwise, computing corners and adding 
c     # vertices to polygon, until we come to the "short" side that
c     # is cut off by the boundary.
c     
         kside = 2
         do 40 i=iedge+1, iedge+3
            kside = kside + 1
            ixc1 = ixc0 + idx(i)
            iyc1 = iyc0 + idy(i)
            xc1 = x(ixc1)
            yc1 = y(iyc1)
            f1 = fbody(xc1,yc1)
            if (f1.lt.0.) go to 50
            if (f1.eq.0.) then
               write(6,601)
               stop
            endif
c     
c     # this is a full side.  Add next corner to polygon vertex list:
            poly(kside,1,k) = xc1
            poly(kside,2,k) = yc1
            ixc0 = ixc1
            iyc0 = iyc1
            xc0 = xc1
            yc0 = yc1
 40      continue
         write(6,*) 'error *** in do 40 loop'
         stop
c     
 50      continue
         iedge = i
         if (iedge.gt.4) iedge = iedge-4
         z = zeroin(0.d0,1.d0,fs,tol)
         x1 = xc0 + z*(xc1-xc0)
         y1 = yc0 + z*(yc1-yc0)
c     
c     # add the last two sides to polygon and set up half-space array:
         poly(kside,1,k) = x1
         poly(kside,2,k) = y1
         poly(kside+1,1,k) = x0
         poly(kside+1,2,k) = y0
c     poly(kside+2,1,k) = -1
         poly(kside+2,1,k) = -11
c     
c     # compute area of polygon:
c     
         ar(k) = 0.d0
         yk = poly(1,2,k)
         do 70 i=1,10
c     if (poly(i+1,1,k).eq.-1) go to 80
            if (poly(i+1,1,k).eq.-11) go to 80
            ar(k) = ar(k) + .5d0*((poly(i,2,k)-yk)+(poly(i+1,2,k)-yk))*
     &           (poly(i+1,1,k)-poly(i,1,k))
 70      continue
 80      continue
c     
c     # compute length of final (irregular) side to see if this
c     # polygon should be eliminated:
c     
         hside = dsqrt((x1-x0)**2 + (y1-y0)**2)
         if (hside/dx  .gt. sqrt(epside)) then
c     .        ar(k)/ar(lstgrd) .gt. 12.d0*dx*dx) then
            lfix = .false.
         else
c     # delete this polygon
            lfix = .true.
c     # replace (x1,y1) by nearest corner:
            ixc1 = idnint((x1-xlow)/dx) + 1
            iyc1 = idnint((y1-ylow)/dy) + 1
            x1 = x(ixc1)
            y1 = y(iyc1)
            if (lastk.ne.0) then
               write(13,*) '> fixing up.. old x,y:',poly(lside,1,lastk),
     &              poly(lside,2,lastk)
               write(13,*) '>                area:',ar(lastk)
            endif
c           if (ar(k) .gt. dx*dy/2.0d0) then
            if (ar(k) .gt. .9d0*dx*dy) then
c     # polygon is nearly full: replace by full regular cell
               irr(ixuse,iyuse) = lstgrd
               lsmall = .false.
               write(13,*) '> regular cell',ixuse,iyuse,'  with area',
     &              ar(k)
c     # patch up previous cell: move vertex to corner
               if (lastk.ne.0) then
                  poly(lside,1,lastk) = x1
                  poly(lside,2,lastk) = y1
               endif
            else
c     # area is nearly zero: mark as exterior
               irr(ixuse,iyuse) = -1
               lsmall = .true.
               write(13,*) '> exterior cell',ixuse,iyuse,'  with area',
     &              ar(k)
c     # patch up previous cell: delete vertex lside since 
c     # the corner is already vertex lside-1
               if (lastk.ne.0) then
                  poly(lside,1,lastk) = poly(lside+1,1,lastk)
                  poly(lside,2,lastk) = poly(lside+1,2,lastk)
c     poly(lside+1,1,lastk) = -1
                  poly(lside+1,1,lastk) = -11
                  if (lside.eq.3) then
                     write(6,*) 'error in mkbdry...two-sided polygon'
                     stop
                  endif
               endif
            endif
c     # recompute area and centroid of previous cell:
            if (lastk.ne.0) then
               call centrd(poly(1,1,lastk),xcirr(lastk),ycirr(lastk),
     &              ar(lastk))
c     ar(lastk) = 0.d0
c     yk = poly(1,2,lastk)
c     do 170 i=1,10
c     if (poly(i+1,1,lastk).eq.-1) go to 180
c     ar(lastk) = ar(lastk) + .5d0*((poly(i,2,lastk)-yk)
c     &             +(poly(i+1,2,lastk)-yk))*
c     &                  (poly(i+1,1,lastk)-poly(i,1,lastk))
c     170          continue
c     180         continue
               write(13,*) '>           new x,y:',poly(lside,1,lastk),
     &              poly(lside,2,lastk)
               write(13,*) '>              area:',ar(lastk)
            endif
            go to 110
         endif
c
c
c       # This cell was not eliminated... store k and beginning of last
c       # side for future reference in case next cell is eliminated:
         lastk = k
         lside = kside
c     
c     # compute centroid of polygon:
         call centrd(poly(1,1,k),xcirr(k),ycirr(k),ar(k))
c     
         if (debug) then
            hmax = dmax1(dx,dy)
            volFrac = ar(k)/ar(lstgrd)
            if (volFrac.lt.gradThreshold) then
               charPrint="*"
            else
               charPrint=" "
            endif
            write(13,121) ixuse,iyuse,k,ar(k),volFrac,charPrint
 121        format('ixuse,iyuse,k,area,volFrac:',3i5,2e15.7,a1)
            write(13,*) '   centroid:', xcirr(k),ycirr(k)
            do 61 i=1,10
c     if (poly(i,1,k).eq.-1) go to 62
               if (poly(i,1,k).eq.-11) go to 62
               if (i .eq. 1) then
                  write(13,901) poly(i,1,k),poly(i,2,k)
               else
                  rat = (dsqrt( (poly(i,1,k)-poly(i-1,1,k))**2+
     1                 (poly(i,2,k)-poly(i-1,2,k))**2 )) / hmax
                  write(13,901) poly(i,1,k),poly(i,2,k),rat
               endif
 901           format(2e25.15,4x,e15.7)
 61         continue
         endif
 62      continue
c     
 110     continue
c     
c     # determine next mesh cell boundary passes through:
c     
         go to (111,112,113,114) iedge
 111     ixuse = ixuse-1
         go to 120
 112     iyuse = iyuse+1
         go to 120
 113     ixuse = ixuse+1
         go to 120
 114     iyuse = iyuse-1
 120     continue
         iedge = iedge+2
c     
c     # if we just deleted a small cell then (x1,y1) was shifted to
c     # a corner and iedge should be incremented one more:
         if (lfix .and. lsmall) iedge = iedge+1
c     
         if (iedge.gt.4) iedge = iedge-4
c     
c     # check to see if we're back at the beginning of a loop:
         if (ixuse.eq.ix0 .and. iyuse.eq.iy0) go to 210
c     
c     # check to see if we're at one of the edges of the grid:
c     if (x1.lt.xlow+eps .or. x1.gt.xup-eps) go to 210
c     if (y1.lt.ylow+eps .or. y1.gt.yup-eps) go to 210
         if (ixuse.eq.0 .or. ixuse.gt.mitot .or. 
     &       iyuse.eq.0 .or. iyuse.gt.mjtot) 
     &        go to 210
c     
c     # otherwise, continue on to the adjacent cell
         x0 = x1
         y0 = y1
c     
 200  continue
      write(6,*) 'error in mkbdry **** too many iterations'
      write(6,*) 'ibcell,k:',ibcell,k
      stop
 210  continue
      if (lfix) k = lastk
      nxtirr(k) = 0
      k0 = k
 601  format('possible error in mkbdry... Rewrite fbody so it only', 
     &     /,'       returns positive or negative values, never 0')

      write(13,999) hx*hy
 999  format(" area of full cell is ",e25.15)
      return
      end
