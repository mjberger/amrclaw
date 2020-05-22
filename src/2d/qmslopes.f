      subroutine qmslopes(irr,mitot,mjtot,lwidth,dx,dy,xlow,ylow,
     . lstgrd, numHoods, ncount,areaMin,
     & mioff, mjoff,
     . qMerge, 
     . nvar,
     . qmx,qmy)

! !   INPUTS:
! !   REQUIRED GRID DATA
! !
! !   irr,mitot,mjtot,lwidth,hx,hy,xlow,ylow,lstgrd
! !
! !   REQUIRED MERGING DATA
! !   volMerge, xcmerge, ycmerge, numHoods,
! !
! !   mioff, mjoff: RECONSTRUCTION NEIGHBORHOODS ON MERGING NEIGHS
! !
! !   qmshifts, qMerge
! !
! !   nvar: number of equations
! !
! !   OUTPUT
! !   qmx,qmy: second order merging slopes


      use amr_module
      implicit double precision (a-h, o-z)

      dimension irr(mitot,mjtot), numHoods(mitot,mjtot)
      dimension rtri(3), stri(3), wtri(3)
      dimension rquad(4), squad(4), wquad(4)

      dimension qMerge(nvar,mitot,mjtot)
      dimension qmx(nvar,mitot,mjtot), qmy(nvar,mitot,mjtot)

      !dimension mi(irrsize,10),mj(irrsize,10)
      dimension ncount(mitot,mjtot)

      dimension a(9,9), rhs(9,nvar), b(9), f(9), w(9,nvar),G(9,9)

      dimension mioff(mitot,mjtot), mjoff(mitot,mjtot)

       logical IS_OUTSIDE
       logical OUT_OF_RANGE

       IS_OUTSIDE(x,y) = (x .lt. xlower .or. x .gt. xupper .or.
     .                    y .lt. ylower .or. y .gt. yupper)

       OUT_OF_RANGE(i,j) = (i .lt. 1 .or. i .gt. mitot .or.
     .                      j .lt. 1 .or. j .gt. mjtot)

      ! QUADRATURE RULE ON TRIANGLES
      ! This quadrature rule integrates quadratics exactly
      rtri = (/ 1.d0/6.d0, 2.d0/3.d0,1.d0/6.d0 /)
      stri = (/ 1.d0/6.d0, 1.d0/6.d0,2.d0/3.d0 /)
      wtri = (/ 1.d0/3.d0, 1.d0/3.d0,1.d0/3.d0 /)
      ntriquad = 3

      ! QUADRATURE RULE ON CARTESIAN CELLS
      ! This quadrature rule integrates cubics exactly
      rquad = (/-dsqrt(1.d0/3.d0),  dsqrt(1.d0/3.d0),
     .           dsqrt(1.d0/3.d0), -dsqrt(1.d0/3.d0)/)
      squad = (/ -dsqrt(1.d0/3.d0), -dsqrt(1.d0/3.d0),
     .            dsqrt(1.d0/3.d0),  dsqrt(1.d0/3.d0) /)
      wquad = (/ 1.d0/4.d0, 1.d0/4.d0, 1.d0/4.d0, 1.d0/4.d0 /)
      nquadquad = 4


      qmx = 0.d0
      qmy = 0.d0

      !jst = max(2,lwidth/2)
      !jend = min(mjtot-1,mjtot-lwidth/2)
      !ist = max(2,lwidth/2)
      !iend = min(mitot-1,mitot-lwidth/2)
      ! this only works if the first "real" cells only needs info
      ! from no more than 2 away. Otherwise need to make larger num ghost cell
      ! need to add a test
      ist = lwidth
      jst = lwidth
      iend = mitot-lwidth+1
      jend = mjtot-lwidth+1

        do 820 j = jst, jend
        do 820 i = ist, iend
             k = irr(i,j)
             if (k .eq. -1 .or. k .eq. lstgrd) cycle
             if (ar(k) .gt. areaMin) cycle
             call getCellCentroid(lstgrd,i,j,xc,yc,xlow,ylow,dx,dy,k)
             if (IS_OUTSIDE(xc,yc) .or. OUT_OF_RANGE(i,j)) cycle

             rhs = 0.d0 ! initialize for accumulation
             a = 0.d0
             x0 = xcentMerge(k)
             y0 = ycentMerge(k)


            do 834 joff = -mjoff(i,j), mjoff(i,j)
            do 834 ioff = -mioff(i,j), mioff(i,j)
                 if (ioff .eq. 0 .and. joff .eq. 0) go to 834
                 if (OUT_OF_RANGE(i+ioff,j+joff)) go to 834
                 koff = irr(i+ioff,j+joff)
                 if (koff .eq. -1) goto 834
                 call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
     &                                xlow,ylow,dx,dy,koff)
                 if (IS_OUTSIDE(xcn,ycn)) go to 834

                 if(koff .eq. lstgrd) then
                     xcoff = xlow + (i+ioff-0.5d0)*dx
                     ycoff = ylow + (j+joff-0.5d0)*dy
                     deltax = xcoff - x0
                     deltay = ycoff - y0
                 else
                     deltax = xcentMerge(koff) - x0
                     deltay = ycentMerge(koff) - y0
                 endif

                 f(1) = deltax/dx
                 f(2) = deltay/dy

                 f(3) = -qmshifts(1,k)
                 f(4) = -qmshifts(2,k)
                 f(5) = -qmshifts(3,k)

                 if (koff .eq. lstgrd) then ! set up data on whole cell
                    ! AG orig
                    !nc(koff)  = 1
                    !mi(koff, 1) = i+ioff
                    !mj(koff, 1) = j+joff
                    !volMerge(koff) = dx * dy / numHoods(i+ioff, j+joff)
                    cvm = dx * dy / numHoods(i+ioff, j+joff)
                 else
                    cvm = volMerge(koff) ! off vol merge
                 endif


                 ! compute weighted inner product of monomials on this neighborhood
                 do 897 ic = 1, ncount(i+ioff,j+joff)+1 
                    ! AG orig code
                    !icurr = mi(koff,ic)
                    !jcurr = mj(koff,ic)
                    if (ic .eq. 1) then !nhood is itself
                       icurr = i+ioff
                       jcurr = j+joff
                       kcurr = irr(icurr,jcurr)
                    else if (ic .eq. 2) then  !nbor had its own nhood
                       kcurr = irr(i+ioff,j+joff)
                       icurr = i+ioff + svi(kcurr)
                       jcurr = j+joff + svj(kcurr)
                       kcurr = irr(icurr,jcurr)
                    else
                       write(*,*)"shouldnt have ic = ",i3," in qmslopes"
                    endif

                    nhc = numHoods(icurr, jcurr) ! num hoods current
                    if (kcurr .eq. lstgrd) then
                      itri = 2
                      poly(1,1,kcurr) = xlow + (dfloat(icurr)-1.d0)*dx
                      poly(1,2,kcurr) = ylow + (dfloat(jcurr)-1.d0)*dy

                      poly(2,1,kcurr) = xlow + (dfloat(icurr)-0.d0)*dx
                      poly(2,2,kcurr) = ylow + (dfloat(jcurr)-1.d0)*dy

                      poly(3,1,kcurr) = xlow + (dfloat(icurr)-0.d0)*dx
                      poly(3,2,kcurr) = ylow + (dfloat(jcurr)-0.d0)*dy

                      poly(4,1,kcurr) = xlow + (dfloat(icurr)-1.d0)*dx
                      poly(4,2,kcurr) = ylow + (dfloat(jcurr)-0.d0)*dy
                    else
                     ivert = 1
                      do 890 while (poly(ivert+1,1,kcurr) .ne. -11.)
                        ivert = ivert + 1
  890                   continue
                       itri = ivert - 3
                    endif

      ! computing the inner product on each triangle of neighborhood member
                  idx1 = 1
                  do 891 it = 1, itri ! for each  triangle
                    idx2 = it + 1
                    idx3 = it + 2

                    x1 = poly(idx1,1,kcurr)
                    y1 = poly(idx1,2,kcurr)

                    x2 = poly(idx2,1,kcurr)
                    y2 = poly(idx2,2,kcurr)

                    x3 = poly(idx3,1,kcurr)
                    y3 = poly(idx3,2,kcurr)

                    artri = triangle_area(x1, x2, x3, y1, y2, y3)

                    do 892 itq = 1,ntriquad
                        xval = x1 * rtri(itq) + x2 * stri(itq)
     .                     +  x3 * (1.d0-rtri(itq)-stri(itq))
                        yval = y1 * rtri(itq) + y2 * stri(itq)
     .                     +  y3 * (1.d0-rtri(itq)-stri(itq))

                        f(3) = f(3) + (artri/nhc/cvm) * wtri(itq) *
     .                          ( (xval-x0)**2 ) / (dx**2)

                        f(4) = f(4) + (artri/nhc/cvm) * wtri(itq) *
     .                          ( (xval-x0)/dx )*( (yval-y0)/dy )

                        f(5) = f(5) + (artri/nhc/cvm) * wtri(itq) *
     .                          ( (yval-y0)**2 ) / (dy**2)
 892        continue ! for each quadrature point on each triangle
 891        continue ! for each triangle
 897        continue ! for each subelement of the current neighborhood
        ! integration is complete, let's accumulate

            do ii = 1, 5
              do jj = 1, 5
                a(ii,jj) = a(ii,jj) + f(ii)*f(jj)
              enddo

              rhs(ii,:) = rhs(ii,:)
     .               + f(ii) * (qMerge(:,i+ioff,j+joff) - qMerge(:,i,j))
            enddo

 8334       continue


 834        continue



      ! SOLVE THE NORMAL EQUATIONS USING THE CHOLESKY FACTORIZATION.
            call cholesky(9, 5, a, G)
            do mm = 1, nvar
               call solve(9,5,G,rhs(:,mm))
               qmy(mm,i,j)  =  rhs(2,mm)/dy
               qmx(mm,i,j)  =  rhs(1,mm)/dx
            end do
c           write(*,*)i,j,qmx(1,i,j),0.d0,qmy(1,i,j),2.*y0
c           write(*,900)i,j,qmx(1,i,j),qmy(1,i,j)
 900        format(2i4,4e15.7)

 820    continue ! iterate over merging tiles on entire grid

      return
      end
