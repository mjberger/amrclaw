c
c -----------------------------------------------------------
c
       subroutine getYface(i,jcol,xface,yface,irr,mitot,mjtot,
     .                     xlow,ylow,dx,dy,lstgrd)

       use amr_module
       implicit double precision (a-h,o-z)
       dimension irr(mitot,mjtot)
       logical IS_CUT 

       IS_CUT(k) = ((k .ne. lstgrd) .and. (k .ne. -1))

c need midpoint of y face between 2 cells. 
c    if bottom is irregular need top face
c    (or if top is irregular need bottom face)
c    if at least one is regular, easy case
c    if one is solid, return face isnt used.
c
      if (irr(i,jcol) .eq. lstgrd .or. irr(i,jcol+1) .eq. lstgrd) then  ! face must be regular
             xface = xlow + (dfloat(i)-.5d0)*dx
             yface = ylow + dfloat(jcol)*dy

      else if (IS_CUT(irr(i,jcol)) .and. IS_CUT(irr(i,jcol+1))) then !use bottom cell to find top face
              k = irr(i,jcol)
              yreg =  ylow + jcol*dy
              do 10 kside = 1, 6
                if (poly(kside+2,1,k) .eq. -11) then
                   write(*,*)"missing side in getYface for ",i,jcol
c                  cut cell situation - just return something to not bomb
                   xface = poly(1,1,k)
                   yface = poly(1,2,k)
                   go to 99
                endif
                x1 = poly(kside,1,k)
                y1 = poly(kside,2,k)
                x2 = poly(kside+1,1,k)
                y2 = poly(kside+1,2,k)
                if (x1 .eq. x2) go to 10  ! not a y face
                if (y1 .ne. y2)  then
                  write(*,*)"screwy side in getYface for cell ",i,jcol
                  stop
                endif
                if (dabs(y1-yreg) .lt. dy/1000.)  then
                   xface = .5d0*(x1+x2)     ! midpoint of face
                   yface = y1
                   go to 99
                endif
 10           continue

      else  ! one or both must be solid, no need for real computation

          xface = 0.d0
          yface = 0.d0

      endif

 99   return
      end
