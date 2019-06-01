c
c --------------------------------------------------------------------------------
c
       subroutine cellsClose(ffluxlen,gfluxlen,mitot,mjtot,irr,lstgrd,
     .                       lwidth)
c
       use amr_module
       implicit double precision (a-h, o-z)

       dimension ffluxlen(mitot+1,mjtot+1), gfluxlen(mitot+1,mjtot+1)
       dimension irr(mitot,mjtot)


       do 10 j = lwidth+1, mjtot-lwidth
       do 10 i = lwidth+1, mitot-lwidth
          k = irr(i,j)
          if (k .eq. -1 .or. k .eq. lstgrd) go to 10
c
c         make sure cut cell closes
c
          sumx = ffluxlen(i+1,j) - ffluxlen(i,j) 
          sumy = gfluxlen(i,j+1) - gfluxlen(i,j)

c         next get length of irregular side
          do 20 kside = 1, 6
            if (poly(kside+2,1,k).eq.-11.) then
               x1 = poly(kside,1,k)
               y1 = poly(kside,2,k)
               x2 = poly(kside+1,1,k)
               y2 = poly(kside+1,2,k)
               go to 25 
            endif
 20      continue

 25      rlen = sqrt((y1-y2)**2 + (x1-x2)**2)
         alf = (y2-y1)/rlen
         beta = (x1-x2)/rlen
         sumx = sumx - alf*rlen
         sumy = sumy - beta*rlen
c        if (dabs(sumx).gt.1.d-12*dx .or. dabs(sumy).gt.1.d-12*dy) then
           write(*,900) i,j,sumx,sumy
 900       format(" cut cell ",2i5," doesnt close by ",2e25.15)
c        endif
c
 10      continue

       return
       end
