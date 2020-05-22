c
c ---------------------------------------------------------------------
c
       subroutine mnslopes(irr,mitot,mjtot,lwidth,dx,dy,xlow,ylow,
     &                    lstgrd,numHoods,ncount,areaMin,mioff,mjoff,
     &                    qMerge,nvar,qmx,qmy)

       use amr_module
       implicit double precision (a-h, o-z)
       include "cuserdt.i"

       dimension irr(mitot,mjtot)
       dimension qmx(nvar,mitot,mjtot), qmy(nvar,mitot,mjtot)
       dimension qMerge(nvar,mitot,mjtot), numHoods(mitot,mjtot)
       dimension nCount(mitot,mjtot)
       dimension mioff(mitot,mjtot),mjoff(mitot,mjtot)

       dimension a(5,5),b(5),rhs(5,nvar)
       dimension nborList(35,2)
       character c2

       logical IS_OUTSIDE, OUT_OF_RANGE
       logical quad, nolimiter,verbose
       common /order2/ ssw, quad, nolimiter

       IS_OUTSIDE(x,y) = (x .lt. xlower .or. x .gt. xupper .or.
     .                    y .lt. ylower .or. y .gt. yupper)

       OUT_OF_RANGE(i,j) = (i .lt. 1 .or. i .gt. mitot .or.
     .                      j .lt. 1 .or. j .gt. mjtot)


        do 20 j = lwidth, mjtot-lwidth+1
        do 20 i = lwidth, mitot-lwidth+1
            k = irr(i,j)
            if (k .eq. -1 .or. k .eq. lstgrd) cycle 
            if (ar(k) .gt. areaMin) cycle  ! wont need gradient
            call getCellCentroid(lstgrd,i,j,xc,yc,xlow,ylow,dx,dy,k)
            if (IS_OUTSIDE(xc,yc) .or. OUT_OF_RANGE(i,j)) cycle  

            nborCount = 0 ! if not enough cant use quad, drop to linear

            ! can reconstruct high order poly here, but for now do linear reconstruction
            ! solving least squares problem with all neighboring merged vals. 
            ! Should preprocess the matrix and save
            rhs = 0.d0 ! initialize for accumulation
            a = 0.d0
            x0 = xcentMerge(k)
            y0 = ycentMerge(k)
            ! mioff/mjoff precomputed for stability
            do 22 joff = -mjoff(i,j), mjoff(i,j)
            do 22 ioff = -mioff(i,j), mioff(i,j)
               if (ioff .eq. 0 .and. joff .eq. 0) go to 22 ! no eqn to solve
               if (OUT_OF_RANGE(i+ioff,j+joff)) cycle
               koff = irr(i+ioff,j+joff)
               if (koff .eq. -1) go to 22
               call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
     &                              xlow,ylow,dx,dy,koff)
               if (IS_OUTSIDE(xcn,ycn)) go to 22

               nborCount = nborCount + 1
               nborList(nborCount,1) = i+ioff
               nborList(nborCount,2) = j+joff
               deltax = xcentMerge(koff) - x0
               deltay = ycentMerge(koff) - y0
               if (koff .eq. lstgrd) then  ! tile centroid same as reg cell
                  call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
     &                                 xlow,ylow,dx,dy,koff)
               else
                  xcn = xcentMerge(koff)
                  ycn = ycentMerge(koff)
               endif
               deltax = xcn - x0
               deltay = ycn - y0
               a(1,1) = a(1,1) + deltax*deltax
               a(1,2) = a(1,2) + deltax*deltay
               a(2,2) = a(2,2) + deltay*deltay
               rhs(1,:) = rhs(1,:) + deltax * 
     .                    (qMerge(:,i+ioff,j+joff) - qMerge(:,i,j))
               rhs(2,:) = rhs(2,:) + deltay * 
     .                    (qMerge(:,i+ioff,j+joff) - qMerge(:,i,j))
               ! if (numMergeTerms .eq. 5) then 
               if (igradChoice .eq. 2) then  ! ptwise quadratic
                  a(1,3) = a(1,3) + deltax*deltax**2
                  a(1,4) = a(1,4) + deltax*deltax*deltay
                  a(1,5) = a(1,5) + deltax*deltay**2
                  a(2,3) = a(2,3) + deltay*deltax**2
                  a(2,4) = a(2,4) + deltay*deltax*deltay
                  a(2,5) = a(2,5) + deltay*deltay**2
                  a(3,3) = a(3,3) + deltax**2*deltax**2
                  a(3,4) = a(3,4) + deltax**2*deltax*deltay
                  a(3,5) = a(3,5) + deltax**2*deltay**2
                  a(4,4) = a(4,4) + deltax*deltay * deltax*deltay
                  a(4,5) = a(4,5) + deltax*deltay * deltay**2
                  a(5,5) = a(5,5) + deltay**4
                  rhs(3,:) = rhs(3,:) + deltax**2*
     .                      (qMerge(:,i+ioff,j+joff) - qMerge(:,i,j))
                  rhs(4,:) = rhs(4,:) + deltax*deltay*
     .                      (qMerge(:,i+ioff,j+joff) - qMerge(:,i,j))
                  rhs(5,:) = rhs(5,:) + deltay**2*
     .                      (qMerge(:,i+ioff,j+joff) - qMerge(:,i,j))
               endif
 22          continue

             if (igradChoice .eq. 1) then ! 1st order gradients, 2 unknowns
               c11 = sqrt(a(1,1))
               c12 = a(1,2) / c11
               c22 = sqrt(a(2,2) - c12**2)

               if (c22 .ne. 0.d0) then ! might have to test c11 too?
                  do m = 1, nvar
                    b(1) = rhs(1,m) / c11
                    b(2) = (rhs(2,m) - c12*b(1))/c22
                    qmy(m,i,j) = b(2) / c22
                    qmx(m,i,j) = (b(1) - c12*qmy(m,i,j))/c11
                  end do
               else
                 write(*,*)"found c22=0 for cell ",i,j
                 ! set gradient to zero?
               endif

             ! solve a^t*a  * grad = a^t rhs.  First form cholesky factors.
             ! will have to add robustness checks
             else ! larger matrix, larger cholesky, 
               c11 = sqrt(a(1,1))
               c12 = a(1,2) / c11   
               c13 = a(1,3) / c11   
               c14 = a(1,4) / c11   
               c15 = a(1,5) / c11   

               c22 = sqrt(a(2,2)-c12**2)
               c23 = (a(2,3)-c12*c13)/c22
               c24 = (a(2,4)-c12*c14)/c22
               c25 = (a(2,5)-c12*c15)/c22

               c33 = sqrt(a(3,3)-c13**2-c23**2)
               c34 = (a(3,4)- c13*c14-c23*c24)/c33
               c35 = (a(3,5)- c13*c15-c23*c25)/c33

               c44 = sqrt(a(4,4)-c14**2-c24**2-c34**2)
               c45 = (a(4,5)-c14*c15-c24*c25-c34*c35)/c44
               c55 = sqrt(a(5,5)-c15**2-c25**2-c35**2-c45**2)
                  
             ! solving c^t b = rhs, then c b = gradients
             do m = 1, nvar
                b(1) = rhs(1,m)/c11
                b(2) = (rhs(2,m)-c12*b(1))/c22
                b(3) = (rhs(3,m)-c13*b(1)-c23*b(2))/c33
                b(4) = (rhs(4,m)-c14*b(1)-c24*b(2)-c34*b(3))/c44
                b(5) = (rhs(5,m)-c15*b(1)-c25*b(2)-c35*b(3) -
     &                             c45*b(4))/c55

                w5 = b(5)/c55
                w4 = (b(4)-c45*w5)/c44
                w3 = (b(3)-c35*w5-c34*w4)/c33
                w2 = (b(2)-c25*w5-c24*w4-c23*w3)/c22
                w1 = (b(1)-c15*w5-c14*w4-c13*w3-c12*w2)/c11

                qmy(m,i,j) = w2
                qmx(m,i,j) = w1
             end do

             endif
 20     continue


       return
       end
