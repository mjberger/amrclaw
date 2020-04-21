      subroutine mymethod(q,qold,mitot,mjtot,lwidth,
     &                  dtn,dtnewn,
     &                  dx,dy,nvar,xlow,ylow,mptr,maux,aux,irr,
     &                  lstgrd,ncount,numHoods,vtime,istage,time)

      use amr_module
      implicit double precision (a-h,o-z)

      dimension q(nvar,mitot,mjtot), qold(nvar,mitot,mjtot)
      dimension f(nvar,mitot,mjtot),g(nvar,mitot,mjtot)
      dimension qx(nvar,mitot,mjtot),qy(nvar,mitot,mjtot)
      dimension aux(maux,mitot,mjtot)
      dimension ur(nvar,max(mitot,mjtot)),ul(nvar,max(mitot,mjtot))
      dimension res(nvar,mitot,mjtot)
      dimension ff(nvar,max(mitot,mjtot))
      dimension ffluxlen(mitot+1,mjtot+1),gfluxlen(mitot+1,mjtot+1)

      integer   irr(mitot,mjtot),ncount(mitot,mjtot)
      integer   numHoods(mitot,mjtot)

      common   /RKmethod/ coeff(5),mstage
      common   /order2/ ssw,quad,nolimiter
      include "cuserdt.i"

      logical   debug, vtime 
      dimension firreg(nvar,-1:irrsize)
      dimension resid(nvar)
      dimension fakeState(nvar)
      character ch

      integer    xrp, yrp
      data       debug/.false./
      data       xrp/1/, yrp/0/
      data       pi/3.14159265357989d0/

c
c
c xrp (= 1) to solve x riemann problem, yrp(= 0) for y riemann problem
c
c       Modified 2/3/89 to solve riemann problems one row at a time
c     # 
c
c cell (i,j) owns the fluxes to the left and bottom
c
      msize = max(mitot,mjtot)
      if ( msize .lt. max(mitot,mjtot) ) then
          write(6,*) 'Not enough memory allocated for rowwise flux'
          write(6,*) 'Calculations. Allocated Size ',msize
          write(6,*) 'Space required ',max(mitot,mjtot)
          write(6,*) 'Remember also to allocate the same additional'
          write(6,*) 'space in subroutine vrm '
          stop
       endif

       if (2*mstage .gt. lwidth) then
          write(6,*) 'Need twice as many ghost cells as stages'
          write(6,*) "or run with single large grid and not parallel"
          write(6,*) "number of ghost cells: ",lwidth
          write(6,*) "number of RK stages:   ",mstage
          stop
       endif

       ! if this grid patch is all solid cells skip it
       ! this method only counts interior cells, not ghost cells.
       call countCellType(irr,mitot,mjtot,lwidth,numSolid,numCut,
     .                    numFull,lstgrd)
       nx = mitot-2*lwidth
       ny = mjtot-2*lwidth
       if (numSolid .eq. nx*ny) then
c        write(*,*) "grid ",mptr," all solid" 
         dtnewn = rinfinity     ! output var, need to return someting
         return
       endif
c      if (numFull .eq. nx*ny) write(*,*) "grid ",mptr," all flow" 
c      if (numCut+numSolid+numFull .ne. nx*ny) then
c         write(*,*)"count doesn't work for grid ",mptr
c      endif

       ! in primitive variables
       fakeState(1) = 1.4d0
       fakeState(2) = 0.d0
       fakeState(3) = 0.d0
       fakeState(4) = 1.0d0

c
      ar(-1) = 1.d0   ! fake value for solid area to avoid zero divides
      dtn = .07d0
c
c     if (debug) then
c        ! for debugging, won't work for multiple grids at same level
c        totmass1 = gridconck(q,irr,mitot,mjtot,lwidth,nvar)
c        write(*,909) totmass1
c909     format("           from method initial mass is ",e30.16)
c     endif
c
c     # initialize fluxes:
      firreg(:,-1)  = 0.d0
      f = 0.d0
      g = 0.d0
      qx = 0.d0
      qy = 0.d0

c need routine to set face lengths and  midpoints
         call getirrlen(irr,mitot,mjtot,dtn,dx,dy,lstgrd,
     &               mptr,nvar,ffluxlen,gfluxlen)
c
c  turn on for debugging, or could move to regridding section of code
c        call cellsClose(ffluxlen,gfluxlen,mitot,mjtot,irr,lstgrd,
c    .                    lwidth)
 
c   ::::: 1 rk stage with linear reconstruction follows ::::::
c
c  store primitive variables in f for now
c
      if (iprob .ne. 20) call vctoprm(q,q,mitot,mjtot,nvar)
c     now using bc2amr in conserved variables. copyingn internally
c     from different grids for each stage
c     ### call for exterior bcs at each stage so can use slopes
c           xhigh= xlow + mitot*dx
c           yhigh = ylow + mjtot*dy
c           call pphysbdlin(xlow,xhigh,ylow,yhigh,level,mitot,mjtot,
c    &                      nvar,q,time,dx,dy,qx,qy,irr,lstgrd)
c        endif
      if (ssw .ne. 0.d0) then   ! recalc slopes at each stage
         call slopes(q,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &               xlow,ylow,nvar,mptr)
         call qslopes(q,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &                 xlow,ylow,mptr,nvar)
      endif

c
c  loop through rows of q calculating fluxes one row at time
c  vertical riemann problem first
c
      do 800 jcol = 2, mjtot-2 
c
         do 511 i = lwidth-2, mitot-lwidth+3
            call getYface(i,jcol,xface,yface,irr,mitot,mjtot,
     .               xlow,ylow,dx,dy,lstgrd)
            call getCellCentroid(lstgrd,i,jcol+1,xcentp,ycentp,
     .                      xlow,ylow,dx,dy,irr(i,jcol+1))
            call getCellCentroid(lstgrd,i,jcol,xcent,ycent,
     .                      xlow,ylow,dx,dy,irr(i,jcol))
            do 511 m = 1, nvar
             if (gfluxlen(i,jcol+1) .ne. 0.d0) then  ! real face
                ur(m,i) = q(m,i,jcol+1) + (xface-xcentp)*qx(m,i,jcol+1)
     .                                  + (yface-ycentp)*qy(m,i,jcol+1)
                ul(m,i) = q(m,i,jcol)   + (xface-xcent)*qx(m,i,jcol)
     .                                  + (yface-ycent)*qy(m,i,jcol)
             else
c              ur(m,i) = q(m,i,jcol+1)  ! dont change vals
c              ul(m,i) = q(m,i,jcol)  ! which might lead to bad rp
               ur(m,i) = fakeState(m)
               ul(m,i) = fakeState(m)
             endif
  511    continue
c
c store fluxes in ff vector, copy into permanent flux array
c
         call vrm(ur,ul,ff,lwidth-2,mitot-lwidth+3,yrp,msize,mptr)
c
         do 720 i = lwidth-2, mitot-lwidth+3
         do 720 m = 1, nvar
            g(m,i,jcol+1) = ff(m,i)
  720    continue
c
c
        
 800  continue
c
c
c    Horizontal riemann problems next
c
      do 900 irow = 2, mitot-2
c
         do 611 j = lwidth-2, mjtot-lwidth+3
            call getXface(irow,j,xface,yface,irr,mitot,mjtot,
     .               xlow,ylow,dx,dy,lstgrd)
            call getCellCentroid(lstgrd,irow+1,j,xcentp,ycentp,
     .                      xlow,ylow,dx,dy,irr(irow+1,j))
            call getCellCentroid(lstgrd,irow,j,xcent,ycent,
     .                      xlow,ylow,dx,dy,irr(irow,j))
            do 611 m = 1, nvar
             if (ffluxlen(irow+1,j) .ne. 0.) then ! real face
                ur(m,j) = q(m,irow+1,j) + (xface-xcentp)*qx(m,irow+1,j)
     .                                  + (yface-ycentp)*qy(m,irow+1,j)
                ul(m,j) = q(m,irow,j)   + (xface-xcent)*qx(m,irow,j)
     .                                  + (yface-ycent)*qy(m,irow,j)
             else
c              ur(m,j) = q(m,irow+1,j)
c              ul(m,j) = q(m,irow,j)
               ur(m,j) = fakeState(m)
               ul(m,j) = fakeState(m)
             endif
  611   continue
  
c
c store fluxes in ff 
c
         call vrm(ur,ul,ff,lwidth-2,mjtot-lwidth+3,xrp,msize,mptr)
c
         do 721  m = 1, nvar
         do 721 j = lwidth-2, mjtot-lwidth+3
            f(m,irow+1,j) = ff(m,j)
  721    continue
c
 900  continue
c
c irregflux computes the cut cell bndry flux. since no flow
c  through bndry use eval pressure there.
          call irregFlux(q,firreg,irr,mitot,mjtot,dx,dy,lstgrd,
     .                   xlow,ylow,mptr,qx,qy,lwidth,nvar)
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c  multiply fluxes by mesh spacing. 
c  this zeros out unused fluxes so solid vals dont get updated.
c
         do 580 i = 1, mitot
         do 580 j = 1, mjtot
            f(:,i,j) = f(:,i,j) * ffluxlen(i,j)
            g(:,i,j) = g(:,i,j) * gfluxlen(i,j)
 580     continue
c 
        ! need conserved vars for second stage and for ghost cell
        ! filling for multiple stages
        if (iprob .ne. 20) call vprmtoc(q,mitot,mjtot,nvar)  
c
c      # finite volume update
c
         c      = coeff(istage)
         ar(-1) = 1.d0   ! prevent zero divides for solid cells
         do 917 j = 1+istage, mjtot-istage
         do 917 i = 1+istage, mitot-istage

         k = irr(i,j)
         do m = 1, nvar
            resid(m) = (f(m,i+1,j)-f(m,i,j)+g(m,i,j+1)-g(m,i,j))
     &                 - firreg(m,k)
            res(m,i,j) = resid(m) ! for debugging save both
            if (istage .eq. 1) then
               q(m,i,j) = qold(m,i,j) - dtn/ar(k)*resid(m)
            else
               q(m,i,j) = q(m,i,j) - dtn/ar(k)*resid(m)  ! second stage builds on first
            endif
         end do
 917     continue

c       call checkPhys(q,irr,mitot,mjtot,mptr,istage,
c    .                  lstgrd,'from my_method')
c       call checkPhysInt(q,mitot,mjtot,mptr,istage,
c    .                    lwidth,'from my_method')
c       write(*,*)"done with physCheckInt for grid ",mptr

c     write(*,*)"calling symcheck after update stage ",istage,
c    .          " before srd"
c     call symcheck(q,irr,mitot,mjtot,nvar,lwidth)

c  postprocess for stability of cut cells. c
c  do it in conserved variables for conservation purposes, (but maybe prim better?)
c
         if (numCut .gt. 0)  then
            if (ismp .eq. 1) then
              call srd_cellMerge(q,nvar,irr,mitot,mjtot,qx,qy,lstgrd,
     .                           dx,dy,lwidth,xlow,ylow,istage,
     .                           ncount,numHoods,mptr,ffluxlen,gfluxlen)
c           else if (ismp .eq. 2) then
c              call drd_cellMerge(q,nvar,irr,mitot,mjtot,qx,qy,lstgrd,
c    .                            dx,dy,lwidth,xlow,ylow,istage,
c    .                            ncount,numHoods,mptr)
            endif
         endif

         if (istage .eq. 2) then  
            q = 0.5d0*(q + qold)
            rhodifmax = 0.d0
            do j = lwidth+1, mjtot-lwidth
            do i = lwidth+1, mitot-lwidth
              rhodifmax =  max(rhodifmax,abs(qold(1,i,j)-q(1,i,j)))
            end do
            end do
            write(*,*)" rhodifmax ",rhodifmax
         endif

c
c
c     :: for multistage method need to adjust for fluxes at each stage
c     :: note that not conservative if using amr unless sum fluxes at each stage around boundary

c
c      :: output adjusted mass for multistage scheme
c
       if (debug) then
          totmass2 = gridconck(q,irr,mitot,mjtot,lwidth,nvar)  ! final mass on grid
c         time step adjustment factor depends on stage.
c         continue to accumulate for 2 stage method and change factor
c         :::: next formula ONLY works for 1 and 2 stage methods
       endif
c
c     # output irregular fluxes for debugging purposes:
      if (debug) then
         write(11,*)"grid ",mptr," flux:",k,i,j
         k = lstgrd
         do 810 ik = 1, irrsize
            k = iabs(nxtirr(k))
            if (k.eq.0) go to 822
               write(11,820) k,ixg(k),iyg(k),firreg(1,k)
               do 810 m = 2, nvar
                  write(11,821) firreg(m,k)
  810             continue
  822    continue
      endif
  820 format(3i4,2d22.10)
  821 format(12x,2d22.10)
c
c     # output fluxes for debugging purposes:
      if (debug) then
         write(11,*)"Inviscid: f             g          "
         do 830 i = lwidth+1, mitot-1
            do 830 j = lwidth+1, mjtot-1
               write(11,831) i,j,f(1,i,j),g(1,i,j)
               do 830 m = 2, nvar
                  write(11,832) f(m,i,j),g(m,i,j)
  830             continue
  831          format(2i4,4d20.10)
  832          format(8x, 4d20.10)
      endif

c
c     estimate new time step for next round. even if not used will give cfl 
      if (vtime .and. istage .eq. mstage) then
c        speedmax = 0.d0
c        do 140 j = lwidth+1, mjtot-lwidth
c        do 140 i = lwidth+1, mitot-lwidth
c           k = irr(i,j)
c           if (k .ne. -1) then
c               p = gamma1* (q(4,i,j)- .5d0* (q(2,i,j)**2 +
c    &                       q(3,i,j)**2)/q(1,i,j))
c               c2 = gamma*p/q(1,i,j)
c               if (c2 .le. 0.d0) go to 140
c               c = dsqrt(c2)
c               u = q(2,i,j)/q(1,i,j)
c               v = q(3,i,j)/q(1,i,j)
c               ::: dt using muscl cfl limit
c                !speed = dmax1(dabs(u)+c,dabs(v)+c)  !use max norm for muscl
c                speed = dsqrt(u*u+v*v) +c  
c                speedmax = max(speedmax,speed) 
c                !dtx = dx / (dabs(u)+c)
c                !dty = dy / (dabs(v)+c)
c                !delT = MIN(dtx,dty)
c           endif
c140     continue
c        effh = dx*dy/(2.d0*dx+2.d0*dy)
c        dtnewn = cflcart* effh / speedmax
          call estdt(q,irr,mitot,mjtot,nvar,dx,dy,dtnewn,lwidth,
     &               aux,naux,cfl)
      endif

c     symmetry check across y = 2 for debugging
c     write(*,*)"calling symcheck after srd  stage ",istage
c     call symcheck(q,irr,mitot,mjtot,nvar,lwidth)

      return
      end
c
c --------------------------------------------------------------------
c
      double precision function gridconck(q,irr,mitot,mjtot,lwidth,nvar)
c
      use amr_module, only: ar 
      implicit double precision (a-h,o-z)

      dimension q(nvar,mitot,mjtot), irr(mitot,mjtot)

      totmass = 0.d0

      do 10 j = lwidth+1, mjtot-lwidth
      do 10 i = lwidth+1, mitot-lwidth

         k = irr(i,j)
         if (k .eq. -1) cycle 
         totmass = totmass + q(1,i,j)*ar(k)

 10   continue

      gridconck = totmass

      return
      end
c
c --------------------------------------------------------------------
c
      subroutine symcheck(q,irr,mitot,mjtot,nvar,nghost)

      implicit real*8 (a-h,o-z)

      dimension q(nvar,mitot,mjtot)
      dimension irr(mitot,mjtot)

      jend = (mjtot-nghost)/2
      tol = 1.d-10

      do i = nghost+1,mitot-nghost
      do j = nghost+1,jend
         if (irr(i,j) .eq. -1) cycle
         j2 = mjtot - j + 1
         dif1 = abs(q(1,i,j) - q(1,i,j2))
         dif4 = abs(q(4,i,j) - q(4,i,j2))
         dif2 = abs(q(2,i,j) - q(2,i,j2))
         dif3 = abs(q(3,i,j) + q(3,i,j2))
         if (dif1 .gt. tol .or. dif2 .gt. tol .or.
     .       dif3 .gt. tol .or. dif4 .gt. tol) then
            write(*,100) i,j,dif1,dif2,dif3,dif4
 100        format(2i5," difs ",5e12.5)
         endif
      end do
      end do
      return
      end
