c
c ------------------------------------------------------
c
      subroutine setirr(mitot,mjtot,mptr,quad,gradThreshold,irr)
c
      use amr_module
      implicit double precision (a-h,o-z)
      integer  irr(mitot,mjtot)
      logical  debug/.false./,done
      character ch
c
c lstgrd = starting index into the linked list of irregular info.
c    for grid mptr. In fact, the first entry is the "regular"
c    cell info, and stores the area of the regular cells.
c
c set the irregular cell array for the augmented grid as follows:
c -1 = cell is completely exterior to the domain
c  lstgrd = cell is a regular cell
c  k <> lstgrd but > 0 = cell is irregular, with k = index into the
c                     linked list for that cell's info.
c
c determine this by calling fbody = returns postiive if in the domain,
c  negative if outside. when a cell has a sign change, a boundary
c  segment has been found. must be careful about orientation -
c  exterior of the domain is on the right. when such a change
c  has been found,  call mkbdry to make the list for that
c  boundary piece. must also take care of completely closed
c  loops inside the grid.
c
       irr = iinfinity
       node(lstptr,mptr) = lstget(dummy)
       lstgrd          = node(lstptr,mptr)
       level           = node(nestlevel,mptr)
       hx              = hxposs(level)
       hy              = hyposs(level)
       ar(lstgrd)      = hx * hy
       nxtirr(lstgrd)  = 0
       done            = .false.
c
c initialize irr array. mkbdry will later modify one of the cells
c at a sign change. use the loc. of the cell center for this first pass
c
       xlow  = rnode(cornxlo,mptr) - nghost*hx
       ylow  = rnode(cornylo,mptr) - nghost*hy
       centx = xlow + .5d0*hx
       centy = ylow + .5d0*hy
       do 10 i = 1, mitot
       do 10 j = 1, mjtot 
          x = centx + (i-1)*hx
          y = centy + (j-1)*hy
          d = fbody(x,y)
          if (d .lt. 0.) then
             irr(i,j) = -1
          else
             irr(i,j) = lstgrd
          endif
 10    continue
c
c look for sign changes in irr as an indication of the starting
c or ending location of a boundary segment. only for starting
c segments, call mkbdry for the info. on that bndry segment
c as it affects this grid.  only look at the perimeter to
c find starting and ending segments. only other case is
c segment is completely contained in interior. must be more
c careful here and look at outside cell edge not center loc.
c this is an example why:
c            |   |   |   |        |___|___|___|___|___|
c            |__/|\__|___|        |___|___|__|=|__|___|
c                                 |___|___|__| |__|___|  
c                                            | |
c
c must catch boundary segments that only glancingly touch the
c fine grid, but don't cover cell center. once it's caught,
c mkbdry will fix it up in the interior, as in the ex. on the right.
c do the bottom, top, left, and right sides.
c
       kloc = lstgrd
       fpre = fbody(xlow,ylow)
       do 20 i = 1, mitot
          fpost = fbody(xlow+i*hx,ylow)
          if ((fpost .lt. 0) .and. (fpre .gt. 0)) 
     1         call mkbdry(kloc,i,1,irr,mitot,mjtot,xlow,ylow,
     2                     hx,hy,lstgrd,mptr,gradThreshold)
          fpre  = fpost
 20    continue
      yhigh = ylow + mjtot*hy
      fpre  = fbody(xlow,yhigh)
       do 30 i = 1, mitot
          fpost = fbody(xlow+i*hx,yhigh)
          if ((fpost .gt. 0) .and. (fpre .lt. 0)) 
     1      call mkbdry(kloc,i,mjtot,irr,mitot,mjtot,xlow,ylow,
     2                  hx,hy,lstgrd,mptr,gradThreshold)
          fpre = fpost
 30    continue
       fpre = fbody(xlow,ylow)
       do 40 j = 1, mjtot
          fpost = fbody(xlow,ylow+j*hy)
          if ((fpost .gt. 0) .and. (fpre .lt. 0)) 
     1         call mkbdry(kloc,1,j,irr,mitot,mjtot,xlow,ylow,
     2                     hx,hy,lstgrd,mptr,gradThreshold)
          fpre = fpost
 40    continue
       xhigh = xlow + mitot*hx
       fpre = fbody(xhigh,ylow)
       do 50 j = 1, mjtot
          fpost = fbody(xhigh,ylow+j*hy)
          if ((fpost .lt. 0) .and. (fpre .gt. 0)) 
     1      call mkbdry(kloc,mitot,j,irr,mitot,mjtot,xlow,ylow,
     2                  hx,hy,lstgrd,mptr,gradThreshold)
          fpre = fpost
 50    continue
c
c last job is to check for closed loops completely contained inside
c the grid, so they might not have intersected the boundary.
c also, double check that every boundary segment accounted for.
c
       do 60 iloops = 1, nloops
          x = xloops(iloops)
          y = yloops(iloops)
          if ((x .lt. xlow) .or. (x .gt. xhigh) .or. (y .lt. ylow) .or.
     1      (y .gt. yhigh)) go to 60
c         point inside the augmented grid. has it been acounted for yet?
          i = (x - xlow)/hx + 1.00d0
          j = (y - ylow)/hy + 1.00d0
          if ((irr(i,j) .eq. lstgrd) .or. (irr(i,j) .eq. -1))
     1    call mkbdry(kloc,i,j,irr,mitot,mjtot,xlow,ylow,
     2                hx,hy,lstgrd,mptr,gradThreshold)
 60    continue
c
c  for each cell with a boundary segment, calculate the weight
c  needed for quadratic reconstruction if nec.
c
cc     if (quad) then
cc     ###  compute weights all the time, even if not used
cc
c         do 70 i = 1, mitot
c         do 70 j = 1, mjtot
c           k = irr(i,j)
c           if (k .eq. -1) go to 70
c           if (k .ne. lstgrd) then
c              call weights(poly(1,1,k),xcirr(k),ycirr(k),totxx,totxy,
c     &                     totyy,wt(1,k),ar(k),points(1,1,k),hx,hy,k)
cc             save tots in rest of poly
c              poly( 8,1,k) = totxx
c              poly( 9,1,k) = totxy
c              poly(10,1,k) = totyy
c           else if (.not. done) then
c               call makep(poly(1,1,lstgrd),i,j,xlow,ylow,hx,hy)
c               xc = (poly(1,1,lstgrd) + poly(3,1,lstgrd))/2.d0
c               yc = (poly(1,2,lstgrd) + poly(2,2,lstgrd))/2.d0
c               call weights(poly(1,1,lstgrd),xc,yc,totxx,totxy,totyy,wt,
c     .                  ar(lstgrd),points(1,1,lstgrd),hx,hy,k)
c               poly(8,1,lstgrd)  = totxx
c               poly(9,1,lstgrd)  = totxy
c               poly(10,1,lstgrd) = totyy
c               done = .true.
c            endif
c 70      continue
cc     endif

       if (debug) then
           write(*,*)" setirr irr grid ", mptr
           do j = 1, mjtot
           do i = 1, mitot
             k = irr(i,j)
             if (k .eq. -1) then
                ch = '*'
             else if (k .eq. lstgrd) then
                ch = ' '
             else
                ch = '+'
             endif
             write(*,222)ch,i,j,irr(i,j)
 222        format(A1,2i4,i7)
           end do
           end do
       endif

       return
       end
