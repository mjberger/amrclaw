      subroutine mymethod(q,fm,fp,gm,gp,mitot,mjtot,lwidth,dtn,dtnewn,
     &                  dx,dy,nvar,xlow,ylow,mptr,maux,aux,irr,
     &                  lstgrd,ncount,numHoods,vtime)

      use amr_module
      implicit double precision (a-h,o-z)

      dimension q(nvar,mitot,mjtot)
      dimension f(nvar,mitot,mjtot),g(nvar,mitot,mjtot)
      dimension fm(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)
      dimension fp(nvar,mitot,mjtot),gp(nvar,mitot,mjtot)
      dimension qx(nvar,mitot,mjtot),qy(nvar,mitot,mjtot)
      dimension ur(nvar,max(mitot,mjtot)),ul(nvar,max(mitot,mjtot))
      dimension qtemp(nvar,mitot,mjtot)
      dimension res(nvar,mitot,mjtot)
      dimension ff(nvar,max(mitot,mjtot))
      dimension ffluxlen(mitot+1,mjtot+1),gfluxlen(mitot+1,mjtot+1)

      integer   irr(mitot,mjtot),ncount(mitot,mjtot)
      integer   numHoods(mitot,mjtot)

      common   /order2/ ssw,quad,nolimiter
      common   /userdt/ cflcart,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                  ismp,gradThreshold

      logical   debug, vtime  
      dimension firreg(nvar,-1:irrsize)
      dimension resid(nvar)
      dimension fakeState(nvar)
      character ch

      integer    xrp, yrp
      data       debug/.true./
      data       xrp/1/, yrp/0/
      data       pi/3.14159265357989d0/

c      dimension coeff(3)
c      data       mstage/3/
c      data       coeff/.21d0,.5d0,1.d0/  

       dimension coeff(2)
       data mstage/2/
       data      coeff/0.5d0,1.d0/

c     dimension coeff(1)
c     data      mstage/1/
c     data      coeff/1.d0/

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

       fakeState(1) = 1.d0
       fakeState(2) = 0.d0
       fakeState(3) = 0.d0
       fakeState(4) = 1.0d0

c
      ix1 = lwidth + 1
      ixn = mitot - lwidth
      iy1 = lwidth + 1
      iyn = mjtot - lwidth
      ar(-1) = 1.d0   ! fake value for solid area to avoid zero divides
c
      if (debug) then
         ! for debugging, won't work for multiple grids at same level
         totmass1 = gridconck(q,irr,mitot,mjtot,lwidth,nvar)
         write(*,909) totmass1
 909     format("           from method initial mass is ",e30.16)
      endif

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
c  could turn on for debugging, or move to regridding section of code
c        call cellsClose(ffluxlen,gfluxlen,mitot,mjtot,irr,lstgrd,
c    .                    lwidth)
 
c
      istage = 1        

 12   continue

c   :::::   rk with linear reconstruction follows ::::::
c
c  store primitive variables in f for now
c
      if (iprob .ne. 20) call vctoprm(q,q,mitot,mjtot,nvar)
c     ### call for exterior bcs at each stage so can use slopes
            xhigh= xlow + mitot*dx
            yhigh = ylow + mjtot*dy
         if (istage .eq. 1 .or. istage .eq. 2) then
            call pphysbdlin(xlow,xhigh,ylow,yhigh,level,mitot,mjtot,
     &                      nvar,q,time,dx,dy,qx,qy,irr,lstgrd)
         endif
      if (ssw .ne. 0.d0) then   ! recalc slopes at each stage
         call slopes(q,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &               xlow,ylow,nvar)
         call qslopes(q,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &                 xlow,ylow,mptr,nvar)
      endif

      if (istage .eq. 1) then  ! moved to here so that bcs in q from new bc routine not phsybd
c        copy q to qtemp for multi-stage RK
         qtemp = q
c        need to convert back to conserved vars for updating below
         call vprmtoc(qtemp,mitot,mjtot,nvar)
       endif
c
c
c  loop through rows of q calculating fluxes one row at time
c  vertical riemann problem first
c
      do 800 jcol = lwidth-2, mjtot-lwidth+3
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
         call vrm(ur,ul,ff,lwidth-2,mitot-lwidth+3,yrp,msize)
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
      do 900 irow = lwidth-2, mitot-lwidth+3
c
         xright = xlow + (dfloat(irow)-.5d0)*dx
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
         call vrm(ur,ul,ff,lwidth-2,mjtot-lwidth+3,xrp,msize)
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
         do 580 i = 2, mitot
         do 580 j = 2, mjtot
            f(:,i,j) = f(:,i,j) * ffluxlen(i,j)
            g(:,i,j) = g(:,i,j) * gfluxlen(i,j)
 580     continue
c 
        if (istage .eq. mstage .and. iprob .ne. 20) then
             call vprmtoc(q,mitot,mjtot,nvar)  ! need conserved vars for last stage only
        endif
c
c
c      # finite volume update
c
         c      = coeff(istage)
         ar(-1) = 1.d0   ! prevent zero divides for solid cells
         do 917 j = iy1-lwidth+istage, iyn+lwidth-istage
         do 917 i = ix1-lwidth+istage, ixn+lwidth-istage
         do 917 m = 1, nvar
            k = irr(i,j)
            resid(m) = (f(m,i+1,j)-f(m,i,j)+g(m,i,j+1)-g(m,i,j))
     &                 - firreg(m,k)
            res(m,i,j) = resid(m)
            if (istage .eq. 1) then
               q(m,i,j) = qtemp(m,i,j) - dtn/ar(k)*resid(m)
            else
c              q(m,i,j) = .5*(qtemp(m,i,j) + q(m,i,j)   ! this way is final answer
c    &                        - dtn/ar(k)*resid(m))
               q(m,i,j) = q(m,i,j) - dtn/ar(k)*resid(m)  ! second stage alone
            endif
 917     continue

c  postprocess for stability of cut cells. c
c  do it in conserved variables for conservation purposes, (but maybe prim better?)
c
         if (ismp .eq. 1) then
            call srd_cellMerge(q,nvar,irr,mitot,mjtot,qx,qy,lstgrd,dx,
     .                  dy,lwidth,xlow,ylow,istage,ncount,numHoods,mptr)
         else if (ismp .eq. 2) then
            call drd_cellMerge(q,nvar,irr,mitot,mjtot,qx,qy,lstgrd,dx,
     .                  dy,lwidth,xlow,ylow,istage,ncount,numHoods)
         endif

         if (istage .eq. 2) then  ! if did above update the 2nd way
            q = 0.5d0*(q + qtemp)
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
         if (istage .lt. mstage ) then
            istage = istage + 1
            go to 12
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
      if (vtime) then
         arreg = dx*dy  ! get regular cell info  
         rlen = dsqrt(arreg)
         dt3 = 1.d10  ! initialize
         do 140 j = lwidth+1, mjtot-lwidth
         do 140 i = lwidth+1, mitot-lwidth
            k = irr(i,j)
            if (k .ne. -1) then
                p = gamma1* (q(4,i,j)- .5d0* (q(2,i,j)**2 +
     &                       q(3,i,j)**2)/q(1,i,j))
                c2 = gamma*p/q(1,i,j)
                if (c2 .le. 0.d0) go to 140
                c = dsqrt(c2)
                u = q(2,i,j)/q(1,i,j)
                v = q(3,i,j)/q(1,i,j)
c               ::: dt using muscl cfl limit
                 speed = dmax1(dabs(u)+c,dabs(v)+c)  !use max norm for muscl
                 dtx = dx / (dabs(u)+c)
                 dty = dy / (dabs(v)+c)
                 delT = MIN(dtx,dty)
                 dt3 =  min(cfl* delT,dt3)
            endif
 140     continue
         dtnewn = dt3     ! output var, need to return someting
      else
         dtnewn = dtn     ! output var, need to return someting
      endif

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
      ar(-1)  = 0.d0

      do 10 j = lwidth+1, mjtot-lwidth
      do 10 i = lwidth+1, mitot-lwidth

         k = irr(i,j)
         totmass = totmass + q(1,i,j)*ar(k)

 10   continue

      gridconck = totmass

      return
      end
