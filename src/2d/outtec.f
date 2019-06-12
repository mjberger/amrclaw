c
c -------------------------------------------------------------
c
      subroutine outtec(q,nvar,mptr,irr,
     1                  mitot,mjtot,lstgrd,
     2                  dx,dy,xlow,ylow,time,ncount,numHoods)
c
      use amr_module
      implicit double precision (a-h,o-z)

      dimension q(nvar,mitot,mjtot) 
      ! use temporary array for qp to avoid converting back
      dimension qp(nvar,mitot,mjtot)  
      integer ncount(mitot,mjtot), numHoods(mitot,mjtot) 
      integer irr(mitot,mjtot) 
      dimension qx(nvar,mitot,mjtot),qy(nvar,mitot,mjtot)
      dimension valprim(4)
      dimension exactsoln(1)
      common /userdt/ cflcart,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold,pwconst
      common /order2/ ssw, quad, nolimiter
      logical  quad, nolimiter
      logical pwconst
      logical isAllSolid,checkIfAllSolid
      character ch
c
c     first count to see if any non solid cells on grid
c     tecplot doesnt like zeroes
      isAllSolid = checkIfAllSolid(irr,mitot,mjtot,nghost)
      call countCellType(irr,mitot,mjtot,nghost,numSolid,numCut,
     &                         numFull,lstgrd)
      if (isAllSolid) then
        if (numSolid .ne. (mitot-2*nghost)*(mjtot-2*nghost)) then
          write(*,*)"count off in outtec"
          stop
        else
          return
        endif
      endif


 8    xlowb = xlow - nghost*dx
      ylowb = ylow - nghost*dy

      ! output primitive variables, not conserved
      if (iprob .ne. 20)  call vctoprm(q,qp,mitot,mjtot,nvar)
c
c     ### call for exterior bcs at each stage so can use slopes
c    ## NOTE THAT BNDRY CELLS FROM OTHER GRIDS NOT SET
            xhigh = xlowb + mitot*dx
            yhigh = ylowb + mjtot*dy
            call pphysbdlin(xlowb,xhigh,ylowb,yhigh,level,mitot,mjtot,
     &                   nvar,qp,time,dx,dy,qx,qy,irr,lstgrd)
 
       qx = 0.d0
       qy = 0.d0

c  set pwconst true for piecewise constant plots, set to false for slopes in tec output
c     pwconst =  .true.
c     pwconst =  .false.
      if (pwconst) go to 9

      if (ssw .ne. 0.d0) then
        call slopes(qp,qx,qy,mitot,mjtot,irr,lstgrd,nghost,dx,dy,
     &               xlowb,ylowb,nvar)
        call qslopes(qp,qx,qy,mitot,mjtot,irr,lstgrd,nghost,dx,dy,
     &                 xlowb,ylowb,mptr,nvar)
      endif

 9    continue

c  count needed for unstructured tec format (so dont have to look up new format)
      nCellsinPlane = 0  
      do 10 i = nghost+1, mitot-nghost
      do 10 j = nghost+1, mjtot-nghost
c      do 10 i = nghost-1, mitot-nghost+2
c      do 10 j = nghost-1, mjtot-nghost+2
         if (irr(i,j) .ne. -1) then
            nCellsinPlane = nCellsinPlane+1
         endif
 10   continue
c
      if (iprob .eq. 20) then
         write(14,101) 4*nCellsinPlane,nCellsinPlane
 101     format('VARIABLES = x,y,Rho,U,V,Pressure,Xcent,Ycent,Err',/,
     .          'Zone T="Cut",N =',i8,' E= ',i8,' F=FEPOINT')
      else
         write(14,103) 4*nCellsinPlane,nCellsinPlane
 103     format('VARIABLES = x,y,Rho,U,V,Pressure,Xcent,Ycent,',
     .                      'ncount,numHoods,i,j,k,volFrac,perim',/,
     .          'Zone T="Cut",N =',i10,' E= ',i10,' F=FEPOINT')
      endif

c  only output real rows and cols, no ghost cells 
c
      do 15 i = nghost+1, mitot-nghost
      do 15 j = nghost+1, mjtot-nghost
c       do 15 i = nghost-1,mitot-nghost+2
c       do 15 j = nghost-1,mjtot-nghost+2
         kirr = irr(i,j)
         if (kirr .eq. -1) go to 15

c        # test again for kirr -1 in case want to view in tecplot
         if (kirr .eq. lstgrd .or. kirr .eq. -1) then
            xcen = xlowb + (dfloat(i)-.5d0)* dx
            ycen = ylowb + (dfloat(j)-.5d0)* dy
            ch = ' '
         else 
            xcen = xcirr(kirr)
            ycen = ycirr(kirr)
            ch = '+'
         endif

         volFrac = ar(kirr)/ar(lstgrd)
         perim =  0.d0 ! maybe perim will indicate reason for instability?
         if (kirr .eq. lstgrd) then 
           perim = 2.d0*(dx+dy)
         else if (kirr .ne. -1) then
           do kside = 1,6
             if (poly(kside+1,1,kirr).eq.-11) go to 62
             hside=dsqrt((poly(kside,1,kirr)-poly(kside+1,1,kirr))**2 +
     &                   (poly(kside,2,kirr)-poly(kside+1,2,kirr))**2)
             perim = perim + hside
           end do
         endif

62       continue

         xcorner = xlowb + (dfloat(i)-1.)* dx    
         ycorner = ylowb + (dfloat(j)-1.)* dy
         do itimes = 1,4
            if (itimes .eq. 1) then ! get all 4 corners of mesh, pw constant sol
               xc = xcorner
               yc = ycorner
            else if (itimes .eq. 2) then
               xc = xcorner
               yc = ycorner + dy
            else if (itimes .eq. 3) then
               xc = xcorner + dx
               yc = ycorner + dy
            else if (itimes .eq. 4) then
               xc = xcorner + dx
               yc = ycorner 
            endif
c
c  reconstruct to corners to can contour through disjoint dataset
c
         do ivar = 1, nvar
            valprim(ivar) = qp(ivar,i,j)+(xc-xcen)*qx(ivar,i,j) +
     .                                   (yc-ycen)*qy(ivar,i,j)
         end do
c        if (iprob .eq. 20)then
c           call p20fn(xcen,ycen,exactsoln,time)
c            errprim =  q(i,j,1) - exactsoln(1)
c        endif

          if (iprob .eq. 20) then  ! output error and soln, nvar=1
             write(14,102) xc,yc,(valprim(ivar),ivar=1,nvar),
     &                  xcen,ycen,errprim
          else if (iprob .eq. 16) then  ! only output soln
c for debugging output error instead
             rhoex = ycen - .1d0*xcen + .5d0
c            rhoex = 1.d0
c            uncomment next line to output error
c            stuffed into w field
             valprim(3) = valprim(1) - rhoex
c            uncomment next line to output density
c            valprim(1) = valprim(1) 
c            valprim(2) = valprim(2)
c            valprim(3) = valprim(3)
c            valprim(4) = valprim(4)
             write(14,1022) xc,yc,(valprim(ivar),ivar=1,nvar),
     &                  xcen,ycen,ncount(i,j),numHoods(i,j),i,j,
     &                  kirr,volFrac,perim
345          format(a1,2i3,2f6.3,4e15.7)

          else if (iprob .eq. 19) then  
            if (kirr .eq. lstgrd) then
               call makep(poly(1,1,kirr),i,j,xlowb,ylowb,dx,dy)
            endif
            call p19tru(xcen,ycen,rhot,ut,vt,pt,poly(1,1,kirr),kirr)
            valprim(1) = valprim(1) - rhot
            valprim(2) = valprim(2) - ut
            valprim(3) = valprim(3) - vt
            valprim(4) = valprim(4) - pt
             write(14,1022) xc,yc,(valprim(ivar),ivar=1,nvar),
     &                  xcen,ycen
          else
            write(14,1022) xc,yc,(valprim(ivar),ivar=1,nvar),
     &                  xcen,ycen,ncount(i,j),numHoods(i,j),i,j,
     &                  kirr,volFrac,perim
          endif
 102      format(9e18.9)
 1022      format(8e25.15,5i8,2e10.2)
        end do

 15    continue
c
c write mesh
c
       ico = 0
       do i = 1, nCellsinPlane
          write(14,104) ico+1,ico+2,ico+3,ico+4
 104      format(4i10)
          ico = ico + 4
       end do

      return
      end
c
c -----------------------------------------------------------
c
      logical function checkIfAllSolid(irr,mitot,mjtot,nghost)

      dimension irr(mitot,mjtot)

      checkIfAllSolid = .true.

      do 5 j = nghost+1, mjtot-nghost
      do 5 i = nghost+1, mitot-nghost
         if (irr(i,j) .ne. -1) then
           checkIfAllSolid = .false.
           go to 99
         endif
 5    continue

 99   return 
      end
c
c -----------------------------------------------------------
c
      subroutine countCellType(irr,mitot,mjtot,nghost,numSolid,numCut,
     &                         numFull,lstgrd)
      
      dimension irr(mitot,mjtot)

      numSolid = 0
      numCut = 0
      numFull = 0

      do 10 j = nghost+1,mjtot-nghost
      do 10 i = nghost+1,mitot-nghost
         if (irr(i,j) .eq. -1) then
            numSolid = numSolid + 1
         else if (irr(i,j) .ne. lstgrd) then
            numCut = numCut + 1
         else
            numFull = numFull + 1
         endif
 10   continue

      return
      end
