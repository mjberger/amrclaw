c
      subroutine mymethod(q,qold,mitot,mjtot,lwidth,
     &                       dtn,dtnewn,
     &                       dx,dy,nvar,xlow,ylow,mptr,maux,aux,irr,
     &                       lstgrd,ncount,numHoods,vtime,istage,time)



      use amr_module
      implicit  double precision (a-h,o-z)

      dimension q(nvar,mitot,mjtot), qold(nvar,mitot,mjtot)
      dimension aux(maux,mitot,mjtot)
      dimension f(nvar,mitot,mjtot),g(nvar,mitot,mjtot)
      dimension qx(nvar,mitot,mjtot),qy(nvar,mitot,mjtot)
      dimension corr(nvar)

      dimension divcoef(mitot,mjtot)
      dimension alldivf(nvar,mitot,mjtot),alldivg(nvar,mitot,mjtot)
      dimension ffluxlen(mitot+1,mjtot+1),gfluxlen(mitot+1,mjtot+1)
      dimension xfaceMidpt(2,mitot+1,mjtot+1)
      dimension yfaceMidpt(2,mitot+1,mjtot+1)

      dimension ur(nvar,max(mitot,mjtot)),ul(nvar,max(mitot,mjtot))
      dimension ff(nvar,max(mitot,mjtot),2)
      integer   irr(mitot,mjtot)
      dimension numHoods(mitot,mjtot),ncount(mitot,mjtot)
      dimension  firreg(nvar,-1:irrsize)

      include "cuserdt.i"
      common /moredata/ sloc

      common  /order2/ ssw,quad,nolimiter
      logical    debug, vtime, tranvd/.true./
      logical   divterm/.true./
      !logical   divterm/.false./

      integer    xrp, yrp, istage
      data       debug/.false./
      data       xrp/1/, yrp/0/
c
c
c xrp (= 1) to solve x riemann problem, yrp(= 0) for y riemann problem
c
c       Modified 2/3/89 to solve riemann problems one row at a time
c     # 
c     # tranvd = true if we want to include transverse derivatives
c
c cell (i,j) owns the fluxes to the left and bottom
c 
      msize = max(mitot,mjtot)

c     ! if this grid patch is all solid cells skip it
      ! this method only counts interior cells, not ghost cells.
       call countCellType(irr,mitot,mjtot,lwidth,numSolid,numCut,
     .                    numFull,lstgrd)

      nx = mitot-2*lwidth
      ny = mjtot-2*lwidth
      if (numSolid .eq. nx*ny) then
c       write(*,*) "grid ",mptr," all solid" 
        dtnewn = rinfinity     ! output var, need to return someting
        return
      endif
c     if (numFull .eq. nx*ny) write(*,*) "grid ",mptr," all flow" 
c     if (numCut+numSolid+numFull .ne. nx*ny) then
c        write(*,*)"count doesn't work for grid ",mptr
c     endif

      ar(-1) = 1.d0
c     # initialize fluxes
      firreg(:,-1) = 0.d0
      firreg(:,lstgrd) = 0.d0
      f = 0.d0
      g = 0.d0
      qx = 0.d0
      qy = 0.d0
c need routine to set face lengths and  midpoints
       call getirrlen(irr,mitot,mjtot,dtn,dx,dy,lstgrd,
     &                mptr,nvar,ffluxlen,gfluxlen)
      call getFaceMidpts(xfaceMidpt,yfaceMidpt,irr,
     &                   dx,dy,xlow,ylow,lstgrd,mitot,mjtot)

c
c
c   :::::   second order Muscl scheme follows ::::::
c
c  work in primitive variables 
c
      call vctoprm(q,q,mitot,mjtot,nvar)
      if (ssw .ne. 0.d0) then
         call slopes(q,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,dx,dy,
     &               xlow,ylow,nvar,mptr)
      if (numCut .gt. 0) 
     &    call qslopes(q,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,
     &                 dx,dy,xlow,ylow,mptr,nvar)
      endif
c
c vertical rp gives transverse fluxes for midstate predction
c for very first row.
c     prepare for Riemann solve for transverse derivative

c  Still need to adjust for irregular in these cols too
      do 510 i = 1, mitot
            ur(:,i) = q(:,i,2)
            ul(:,i) = q(:,i,1) 
 510  continue
c
c store fluxes in ff 
c     ff(1,1,jtog) contains upper row of fluxes
c     ff(1,1,itog) contains lower row of fluxes
c
      jtog = 1
      itog = 2
      call vrm(ur,ul,ff(1,1,jtog),1,mitot,yrp,msize)
c
c  loop through rows of q calculating fluxes one row at time
c
      do 800 jcol = 2, mjtot-1
c
c previous top row of transverse fluxes becomes bottom row
c
         jtemp = jtog
         jtog  = itog
         itog  = jtemp
c
c vertical rp gives transverse fluxes for midstate prediction for next row
c
         do 511 i = 1, mitot
            ur(:,i) = q(:,i,jcol+1)
            ul(:,i) = q(:,i,jcol) 
            kl = irr(i,jcol)
            kr = irr(i,jcol+1)
            !if (kl .eq. -1 .or. kr .eq. -1) cycle 
            ! if one is solid then that edge is not a real edge 
            if (gfluxlen(i,jcol+1) .eq. 0.d0) cycle
            !  for edges next to irregular cells, reconstruct to edge midpoint
            if (kr .eq. lstgrd)then
                xcenr = xlow + (i-0.5d0)*dx
                ycenr = ylow + (jcol+0.5d0)*dy
            else
                xcenr = xcirr(kr)
                ycenr = ycirr(kr)
            endif
            if (kl .eq. lstgrd) then
               xcenl = xlow + (i-0.5d0)*dx
               ycenl = ylow + (jcol-0.5d0)*dy
            else
               xcenl = xcirr(kl)
               ycenl = ycirr(kl)
            endif
            if (kr .eq. lstgrd .or. kl .eq. lstgrd) then   
            ur(:,i) = ur(:,i) + (yfaceMidpt(1,i,jcol+1)-xcenr)*
     &                qx(:,i,jcol+1)+(yfaceMidpt(2,i,jcol+1)-ycenr)*
     &                qy(:,i,jcol+1)
           
            ul(:,i) = ul(:,i) + (yfaceMidpt(1,i,jcol+1)-xcenl)*
     &                qx(:,i,jcol)+ (yfaceMidpt(2,i,jcol+1)-ycenl)*
     &                qy(:,i,jcol)
            endif
  511    continue
c
c store fluxes in ff 
c
        call vrm(ur,ul,ff(1,1,jtog),1,mitot,yrp,msize)

c
c overwrite ur,ul with half step predictor
c
         call vpredx(q,qx,irr,ur,ul,mitot,mjtot,dtn,dx,dy,jcol,nvar,
     &               msize,lstgrd,qy,xlow,ylow,xfaceMidpt,ffluxlen)
         call rowprmtoc(ul,nvar,msize,1,mitot)
         call rowprmtoc(ur,nvar,msize,1,mitot)
c
c  add transverse derivatives from rp above to get midstates
c  fix up irregular cells in second pass 
c
         hdtdy = .5d0*dtn/dy
         hdt = .5d0 * dtn
         if (tranvd) then  
         do 525 i = 2, mitot-1
c        computing ul/r states for edge i,jcol
c        right state uses transverse deriv from
c        i,jcol+-1, left state from i-1,jcol+-1
         if (ffluxlen(i,jcol) .eq. 0.d0) cycle ! not a real edge
         ! need top and bottom g flux edges to be real to do transverse deriv
!!         if (gfluxlen(i,jcol+1) .gt. 0.d0 .and. 
!!     &       gfluxlen(i,jcol)   .gt. 0.d0) the
!!        if (irr(i,jcol+1) .eq. lstgrd .or. 
!!     &      irr(i,jcol)   .eq. lstgrd) then
          if (irr(i,jcol) .eq. lstgrd) then
             ur(:,i) = ur(:,i) - hdtdy*(ff(:,i,jtog)-ff(:,i,itog))
         else if (irr(i,jcol).ne.lstgrd) then 
            rho = q(1,i,jcol)
            u = q(2,i,jcol)
            v = q(3,i,jcol)
            p = q(4,i,jcol)
            corr(1) = v*qy(1,i,jcol) + rho*qy(3,i,jcol)
            corr(2) = v*qy(2,i,jcol)
            corr(3) = v*qy(3,i,jcol) + qy(4,i,jcol)/rho
            corr(4) = gamma*p*qy(3,i,jcol) + v*qy(4,i,jcol)
            call convert2prim(ur(:,i),nvar)
            ur(:,i) = ur(:,i) - hdt *  corr(:)
            call convert2cons(ur(:,i),nvar)
          endif
!!          if (gfluxlen(i-1,jcol+1) .gt. 0.d0 .and.
!!     &        gfluxlen(i-1,jcol)   .gt. 0.d0) then
!!        if (irr(i-1,jcol+1) .eq. lstgrd .and.
!!   &        irr(i-1,jcol)   .eq. lstgrd) then
          if (irr(i-1,jcol) .eq. lstgrd) then
            ul(:,i) = ul(:,i) - hdtdy*(ff(:,i-1,jtog)-ff(:,i-1,itog))
          else if (irr(i-1,jcol) .ne. lstgrd) then
            rho = q(1,i-1,jcol)
            u = q(2,i-1,jcol)
            v = q(3,i-1,jcol)
            p = q(4,i-1,jcol)
            corr(1) = v*qy(1,i-1,jcol) + rho*qy(3,i-1,jcol)
            corr(2) = v*qy(2,i-1,jcol)
            corr(3) = v*qy(3,i-1,jcol) + qy(4,i-1,jcol)/rho
            corr(4) = gamma*p*qy(3,i-1,jcol) + v*qy(4,i-1,jcol) 
            call convert2prim(ul(:,i),nvar)
            ul(:,i) = ul(:,i) - hdt *  corr(:)
            call convert2cons(ul(:,i),nvar)     
          endif
  525    continue
         endif

c
c  calculate fluxes for the row  - into ff, then store into f.
c
         call vrmc(ur,ul,ff(1,1,itog),1,mitot,xrp,msize)
         do 720 i = 1, mitot
            f(:,i,jcol) = ff(:,i,itog)
  720    continue
c
 800  continue
c
c
c    All done with rows. Now do column problems
c
c vertical rp gives transverse fluxes for midstate prediction
c for very first column.
c
         do 610 j = 1, mjtot
            ur(:,j) = q(:,2,j)
            ul(:,j) = q(:,1,j)  
  610  continue
c
c store fluxes in ff 
c     ff(1,1,jtog) contains front column of fluxes
c     ff(1,1,itog) contains back column of fluxes
c
      jtog = 1
      itog = 2
      call vrm(ur,ul,ff(1,1,jtog),1,mjtot,xrp,msize)
c
c  loop through columns of q calculating fluxes one column at time
c
      do 900 irow = 2, mitot-1
c
c previous front column of transverse fluxes becomes back column
         jtemp = jtog
         jtog  = itog
         itog  = jtemp
c
c vertical rp gives transverse fluxes for midstate prediction for next column
c
         do 611 j = 1, mjtot
            ur(:,j) = q(:,irow+1,j)
            ul(:,j) = q(:,irow,j)
            kl = irr(irow,j)
            kr = irr(irow+1,j)
            if (ffluxlen(irow+1,j) .eq. 0.d0) cycle ! not a real edge

            ! for irregular cell edges reconstruct to edge midpoint
            if (kr .eq. lstgrd)then
                xcenr = xlow + (irow+0.5d0)*dx
                ycenr = ylow + (j-0.5d0)*dy
            else
                xcenr = xcirr(kr)
                ycenr = ycirr(kr)
            endif
            if (kl .eq. lstgrd) then
               xcenl = xlow + (irow-0.5d0)*dx
               ycenl = ylow + (j-0.5d0)*dy
            else
               xcenl = xcirr(kl)
               ycenl = ycirr(kl)
            endif
c           if (kr.eq.lstgrd .or. kl.eq.lstgrd) then
            ur(:,j) = ur(:,j) + (xfaceMidpt(2,irow+1,j)-ycenr)*
     &                qy(:,irow+1,j)+(xfaceMidpt(1,irow+1,j)-xcenr)*
     &                qx(:,irow+1,j)
            ul(:,j) = ul(:,j) + (xfaceMidpt(2,irow+1,j)-ycenl)*
     &                qy(:,irow,j)+ (xfaceMidpt(1,irow+1,j)-xcenl)*
     &                qx(:,irow,j)
c             endif
  611  continue
c
c store fluxes in ff 
c
         call vrm(ur,ul,ff(1,1,jtog),1,mjtot,xrp,msize)
c
c overwrite ur,ul with half step predictor
c
         call vpredy(q,qy,irr,ur,ul,mitot,mjtot,dtn,dx,dy,irow,nvar,
     &               msize,lstgrd,qx,xlow,ylow,yfaceMidpt,gfluxlen)
         call rowprmtoc(ul,nvar,msize,1,mjtot)
         call rowprmtoc(ur,nvar,msize,1,mjtot)
c
c  add transverse derivatives from rp above to get midstates
c  again - no slope but yes transverse deriv. at irreg. cells
c
         hdtdx = .5d0*dtn/dx
         hdt = .5d0*dtn
         if (tranvd) then
         do 626 j = 2, mjtot-1
c        making l/r states for flux at edge irow,j
c        right (top) state uses transverse deriv from
c        (irow+-1,j) - (irow+-1,j-1)
         if (gfluxlen(irow,j) .eq. 0.d0) cycle ! not a real edge
c        else at least we know cells on both sides of edge exist so use
c        transverse deriv (but might be inaccurate and need correction)
         ! need top and bottom g flux edges to be real to do transverse deriv
!!         if (ffluxlen(irow+1,j) .gt. 0.d0 .and. 
!!     &       ffluxlen(irow,j)   .gt. 0.d0) then
!!         if (irr(irow+1,j) .eq. lstgrd .or.
!!     &       irr(irow,j)   .eq. lstgrd) then
           if (irr(irow,j) .eq. lstgrd) then
             ur(:,j) = ur(:,j) - hdtdx * (ff(:,j,jtog)-ff(:,j,itog))
         ! correct df/dy for different in ymidpts at x faces
            else if (irr(irow,j) .ne. lstgrd) then
                rho = q(1,irow,j)
                u = q(2,irow,j)
                v = q(3,irow,j)
                p = q(4,irow,j)
                corr(1) = u*qx(1,irow,j) + rho*qx(2,irow,j)
                corr(2) = u*qx(2,irow,j) +qx(4,irow,j)/rho
                corr(3) = u*qx(3,irow,j)
                corr(4) = gamma*p*qx(2,irow,j) + u*qx(4,irow,j)
                call convert2prim(ur(:,j),nvar)
                ur(:,j) = ur(:,j) - hdt *  corr(:)
                call convert2cons(ur(:,j),nvar)
          endif
!!          if (ffluxlen(irow+1,j-1) .gt. 0.d0 .and. 
!!     &       ffluxlen(irow,j-1)   .gt. 0.d0) then
!!          if (irr(irow+1,j-1) .eq. lstgrd .or. 
!!     &       irr(irow,j-1)   .eq. lstgrd) then  ! must be full f face
            if (irr(irow,j-1) .eq. lstgrd) then
             ul(:,j) = ul(:,j) - hdtdx * (ff(:,j-1,jtog)-ff(:,j-1,itog))
         ! correct df/dy for different in ymidpts at x faces
             else if (irr(irow,j-1) .ne. lstgrd) then
                rho = q(1,irow,j-1)
                u = q(2,irow,j-1)
                v = q(3,irow,j-1)
                p = q(4,irow,j-1)
! corr is df/dy at centroid
                corr(1) = u*qx(1,irow,j-1) + rho*qx(2,irow,j-1)
                corr(2) = u*qx(2,irow,j-1) +qx(4,irow,j-1)/rho
                corr(3) = u*qx(3,irow,j-1)
                corr(4) = gamma*p*qx(2,irow,j-1) + u*qx(4,irow,j-1)
                call convert2prim(ul(:,j),nvar)
                ul(:,j) = ul(:,j) - hdt *  corr(:)
                call convert2cons(ul(:,j),nvar)
          endif
 626     continue
         endif
c
c  calculate fluxes for the column into ff, then store into g
c
         call vrmc(ur,ul,ff(1,1,itog),1,mjtot,yrp,msize)
         do 721 j = 1, mjtot
            g(:,irow,j) = ff(:,j,itog)
  721    continue
c
 900  continue
c
c    All done with columns 
c
c irregflux computes the cut cell bndry flux. since no flow
c  through bndry use eval pressure there.
      if (numCut .gt. 0)
     &  call irregFlux(q,firreg,irr,mitot,mjtot,dx,dy,lstgrd,
     &                 xlow,ylow,mptr,qx,qy,lwidth,nvar)
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c  multiply fluxes by mesh spacing. 
c  this zeros out unused fluxes so solid vals dont get updated.
c
c    still need to modify for irregular cells
c
      if (divterm) then
c          q in primitive vars, qold still in conserved vars
c          divcorn routine works with prims
c          but result applied to conserved
c          call orig_divcorn(q,divcoef,mitot,mjtot,dx,dy,mptr,nvar)
         call divcorn(q,qold,divcoef,mitot,mjtot,dx,dy,mptr,nvar,lstgrd,
     &                  xlow,ylow,irr,nghost,qx,qy,
     &                  alldivf,alldivg)
           do 580 j = 2, mjtot-1
           do 580 i = 2, mitot-1
           if (irr(i,j) .eq. -1) cycle

           if (irr(i-1,j) .ne. -1)
     &         f(:,i,j) = f(:,i,j) - alldivf(:,i,j)
c     &         f(:,i,j) = f(:,i,j)-.5d0*
c     &         (divcoef(i,j)+divcoef(i,j+1))*
c     &         (qold(:,i,j)-qold(:,i-1,j))

           if (irr(i,j-1) .ne. -1)
     &          g(:,i,j) = g(:,i,j) - alldivg(:,i,j)
c     &         g(:,i,j) = g(:,i,j)-.5d0*
c     &          (divcoef(i,j)+divcoef(i+1,j))*
c     &          (qold(:,i,j)-qold(:,i,j-1))
 580       continue
      endif

      do 585 j = 1, mjtot
      do 585 i = 1, mitot
         f(:,i,j) = f(:,i,j) * ffluxlen(i,j)
         g(:,i,j) = g(:,i,j) * gfluxlen(i,j)
 585  continue

c
c     # update q values by differencing fluxes. q in prim, qold in cons, overwrite q.
c
      ar(-1) = 1.d0
      do 917 j = lwidth-1, mjtot-lwidth+1
      do 917 i = lwidth-1, mitot-lwidth+1
         k = irr(i,j)
         if (k .ne. -1) then  ! dont bother updating solid cells 
           do m = 1, nvar
              q(m,i,j) = qold(m,i,j) - dtn/ar(k)* (f(m,i+1,j) - f(m,i,j)
     &                   + g(m,i,j+1) - g(m,i,j) - firreg(m,k))
           end do
         endif
 917   continue

       call checkPhys(q,irr,mitot,mjtot,mptr,istage,
     &                lstgrd,'from my_method',
     &                lwidth,mitot-lwidth+1,lwidth,mjtot-lwidth+1)

c      if (numCut .gt. 0)  then
c         if (ismp .eq. 1) then
c           call srd_cellMerge(q,nvar,irr,mitot,mjtot,qx,qy,lstgrd,
c    .                         dx,dy,lwidth,xlow,ylow,istage,
c    .                         ncount,numHoods,mptr,ffluxlen,gfluxlen)
c         endif
c      endif

c     # output fluxes for debugging purposes:
      if (debug) then
          write(11,*)' fluxes for grid ',mptr
          k = lstgrd
          do 810 ik = 1, 3933
             k = iabs(nxtirr(k))
             if (k.eq.0) go to 822
                 write(11,820) k,ixg(k),iyg(k),firreg(1,k)
                 do 810 m = 2, 4
                    write(11,821) firreg(m,k)
  810               continue
  822     continue

          write(11,*)"Inviscid: f             g          "
          do 830 i = lwidth+1, mitot-1
          do 830 j = lwidth+1, mjtot-1
             write(11,831) i,j,f(1,i,j),g(1,i,j)
             do 830 m = 2, nvar
                write(11,832) f(m,i,j),g(m,i,j)
  830           continue
       endif
  820 format(3i4,d16.6)
  821 format(12x,d16.6)
  831 format(2i4,4d20.10)
  832 format(8x, 4d20.10)
c
c compute max speed for variable time stepping
c
c     if (vtime) then
c        call estdt(q,irr,mitot,mjtot,nvar,dx,dy,dtnewn,lwidth,
c    &              aux,naux,cfl)
c     endif
c
      return
      end
c
c -------------------------------------------------------------
c
      subroutine rowprmtoc(u,nvar,msize,ixmin,ixmax)

c     need special routine to convert only a row, with different dimensions than entire grid patch

      implicit double precision (a-h,o-z)
      dimension u(nvar,msize)
      include "cuserdt.i"


       do 20 k = ixmin,ixmax
         rho  = u(1,k)
         uvel = u(2,k)
         vvel = u(3,k)
         pr   = u(4,k)
         ei   = pr/(rho*gamma1)
         en = ei + 0.5d0*(uvel*uvel+vvel*vvel)
         u(2,k) = uvel *rho
         u(3,k) = vvel *rho
         u(4,k) = en*rho
 20    continue
c     
       return
       end

c
c -------------------------------------------------------------
c     
       subroutine vpredx(q,qx,irr,ur,ul,mitot,mjtot,dt,dx,dy,
     &                   j,nvar,msize,lstgrd,qy,xlow,ylow,
     &                   xfaceMidpt,ffluxlen)

       use amr_module, only: xcirr,ycirr
       implicit double precision(a-h,o-z)
       include "cuserdt.i"

       dimension q(nvar,mitot,mjtot)
       dimension qx(nvar,mitot,mjtot),qy(nvar,mitot,mjtot)
       dimension irr(mitot,mjtot)
       dimension ur(nvar,msize),ul(nvar,msize)
       dimension xfaceMidpt(2,mitot+1,mjtot+1)
       dimension ffluxlen(mitot+1,mjtot+1)

c
c predictor step for variables at midtime at cell edges for muscl
c working on row number j
c cell i,j own left interface.  when computing interface to right
c the states come from cells  i,j and i+1,j and it is intfc i+1.

c
       !this version has the dx in the slopes so take it out of hdtdx
       !hdtdx = .5d0*dt/dx
       hdt = .5d0*dt  
          do 20 i = 1, mitot-1
             edgeLength = ffluxlen(i+1,j)
             k = irr(i,j)

             rho = q(1,i,j)
             u   = q(2,i,j)
             v   = q(3,i,j)
             pr  = q(4,i,j) 
             cs2 = gamma*pr / rho
             cs  = dsqrt(cs2)
             bt1 = hdt*(qx(4,i,j)/cs-rho*qx(2,i,j))
             bt2 = hdt*cs*(qx(1,i,j)-qx(4,i,j)/cs2)
             bt3 = hdt*cs*qx(3,i,j)
             if (u-cs .le. 0.d0) bt1 = 0.d0
             if (u .le. 0.d0) then
                bt2 = 0.d0
                bt3 = 0.d0
             endif  
             if (k .eq. lstgrd .or. k .eq. -1) then
                xdif = .5d0*dx
             else 
                xdif = xfaceMidpt(1,i+1,j) - xcirr(k)
             endif
             off = xdif  - dmax1(u+cs,0.d0)*hdt
c
             ul(1,i+1) = bt1 + bt2 + rho + off*qx(1,i,j)
             ul(2,i+1) = -bt1*cs/rho + u + off*qx(2,i,j)
             ul(3,i+1) = bt3 + v + off*qx(3,i,j)
             ul(4,i+1) = bt1*cs2 + pr + off*qx(4,i,j)

             if (k .ne. lstgrd .and. k .ne. -1) then
              if (edgeLength .gt. 0.d0) then
                ydif = xfaceMidpt(2,i+1,j) - ycirr(k)
                ul(:,i+1) = ul(:,i+1) + ydif*qy(:,i,j)
              else
                ul(:,i+1) = q(:,i,j)
              endif
            endif

             k = irr(i+1,j)

             rho = q(1,i+1,j)
             u   = q(2,i+1,j)
             v   = q(3,i+1,j)
             pr =  q(4,i+1,j) 
             cs2 = gamma*pr / rho
             cs  = dsqrt(cs2)
             bt2 = hdt*cs*(qx(4,i+1,j)/cs2-qx(1,i+1,j))
             bt3 = -1.d0*hdt*cs*qx(3,i+1,j)
             bt4 = -1.d0*hdt*(qx(2,i+1,j)*rho+qx(4,i+1,j)/cs)
             if (u .gt. 0.d0) then
                bt2 = 0.d0
                bt3 = 0.d0
             endif
             if (u+cs .ge. 0.d0) bt4 = 0.d0

             if (k .eq. lstgrd .or. k .eq. -1) then
                xdif = .5d0*dx
             else 
                xdif =  xcirr(k) - xfaceMidpt(1,i+1,j)
             endif
             off = xdif + dmin1(u-cs,0.d0)*hdt
c
             ur(1,i+1) = bt2 + bt4 + rho - off*qx(1,i+1,j)
             ur(2,i+1) = bt4*cs/rho + u - off*qx(2,i+1,j)
             ur(3,i+1) = bt3 + v - off*qx(3,i+1,j)
             ur(4,i+1) = bt4*cs2 + pr - off*qx(4,i+1,j)

             if (k .ne. lstgrd .and. k .ne. -1) then
               if ( edgeLength .gt. 0.d0) then
                ydif = xfaceMidpt(2,i+1,j) - ycirr(k)
                ur(:,i+1) = ur(:,i+1) + ydif*qy(:,i+1,j)
               else
                ur(:,i+1) = q(:,i+1,j)
             endif
           endif
               
 20       continue
c
       return
       end
c
c -------------------------------------------------------------
c 
       subroutine vpredy(q,qy,irr,ut,ub,mitot,mjtot,dt,dx,dy,
     &                   i,nvar,msize,lstgrd,qx,xlow,ylow,yfaceMidpt,
     &                   gfluxlen)

       use amr_module, only: xcirr,ycirr     
       implicit double precision(a-h,o-z)

       include "cuserdt.i"
       dimension q(nvar,mitot,mjtot)
       dimension qy(nvar,mitot,mjtot),qx(nvar,mitot,mjtot)
       dimension irr(mitot,mjtot)
       dimension ub(nvar,msize),ut(nvar,msize)
       dimension yfaceMidpt(2,mitot+1,mjtot+1)
       dimension gfluxlen(mitot+1,mjtot+1)
c
c predictor step for variables at midtime at cell edges for muscl
c        working on ith column
c
       !hdtdy = .5d0*dt/dy
       ! this version has the dy in the sloeps so take it out of
       ! hdtdy.  Need to adjust when include cut cells
       hdt = 0.5d0*dt ! better variable for new stuff
       do 80 j = 1, mjtot-1
             edgeLength = gfluxlen(i,j+1)
             k = irr(i,j)

             rho = q(1,i,j)
             u   = q(2,i,j)
             v   = q(3,i,j)
             pr  = q(4,i,j) 

             cs2 = gamma*pr/rho
             cs  = dsqrt(cs2)
             bt1 = hdt*(qy(4,i,j)/cs-rho*qy(3,i,j))
             bt2 = hdt*cs*(qy(1,i,j)-qy(4,i,j)/cs2)
             bt3 = -1.d0*hdt*cs*qy(2,i,j)
             if (v-cs .le. 0.d0) bt1 = 0.d0
             if (v .le. 0.d0) then
                bt2 = 0.d0
                bt3 = 0.d0
             endif
             if (k .eq. lstgrd .or. k .eq. -1) then
                ydif = 0.5*dy
             else if (k .ne. -1) then
                ydif = yfaceMidpt(2,i,j+1) - ycirr(k)
             endif
             off = ydif - dmax1(v+cs,0.d0)*hdt
c
             ub(1,j+1) = bt1 + bt2 + rho + off*qy(1,i,j)
             ub(2,j+1) = -bt3 + u + off*qy(2,i,j)
             ub(3,j+1) = -bt1*cs/rho + v + off*qy(3,i,j)
             ub(4,j+1) = bt1*cs2 + pr + off*qy(4,i,j)

             if (k .ne. lstgrd .and. k.ne. -1) then
               if (edgeLength .gt. 0.d0) then
                xdif = yfaceMidpt(1,i,j+1) - xcirr(k)
                ub(:,j+1) = ub(:,j+1) + xdif*qx(:,i,j)
               else
                ub(:,j+1) = q(:,i,j)  ! not real RP
               endif
             endif
c
             k = irr(i,j+1)

             rho = q(1,i,j+1)
             u   = q(2,i,j+1)
             v   = q(3,i,j+1)
             pr =  q(4,i,j+1) 
             cs2 = gamma*pr/rho
             cs  = dsqrt(cs2)
             bt2 =  hdt*cs*(qy(4,i,j+1)/cs2-qy(1,i,j+1))
             bt3 =  hdt*cs*qy(2,i,j+1)
             bt4 =  -1.d0*hdt*(qy(3,i,j+1)*rho+qy(4,i,j+1)/cs)
             if (v .ge. 0.d0) then
                bt2 = 0.d0
                bt3 = 0.d0
             endif
             if (v+cs .ge. 0.d0) bt4 = 0.d0
             if (k .eq. lstgrd .or. k .eq. -1) then
               ydif = 0.5d0*dy
             else 
               ydif = ycirr(k) - yfaceMidpt(2,i,j+1)
             endif
             off = ydif + dmin1(v-cs,0.d0)*hdt
c
             ut(1,j+1) = bt2 + bt4 + rho - off*qy(1,i,j+1)
             ut(2,j+1) = -bt3 + u - off*qy(2,i,j+1)
             ut(3,j+1) = bt4*cs/rho + v - off*qy(3,i,j+1)
             ut(4,j+1) = bt4*cs2 + pr - off*qy(4,i,j+1)

             if (k .ne. lstgrd .and. k .ne. -1) then
               if (edgeLength .gt. 0.d0) then
                xdif = yfaceMidpt(1,i,j+1) - xcirr(k)
                ut(:,j+1) = ut(:,j+1) + xdif*qx(:,i,j+1)
              else
                ut(:,j+1) = q(:,i,j+1)
             endif
            endif

 80    continue
c
       return
       end
c
c ------------------------------------------------------------------
c
       subroutine orig_divcorn(q,div,mitot,mjtot,dx,dy,mptr,nvar)

       implicit double precision (a-h,o-z)
       dimension q(nvar,mitot,mjtot),div(mitot,mjtot)
c      dimension u(mitot,mjtot),v(mitot,mjtot)
c      parameter (qm=.2d0)

c :::::::::::::::::::::: DIVCORN :::::::::::::::::::::::::::::::::::::;
c  compute divergence to use in artificial viscosity in method
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

c  ::: convert to velocities
c  THIS VERISON ALREADY IN PRIMITIVE VARS
c      do 10 j = 2, mjtot-1
c      do 10 i = 2, mitot-1
c         u(i,j) = q(i,j,2)/q(i,j,1)
c         v(i,j) = q(i,j,3)/q(i,j,1)
c         div(i,j)       = 0.d0
c10    continue

       div       = 0.d0

c
c :::  cell i,j owns lower left corner divergence
c
       do 20 j = 2, mjtot-1
       do 20 i = 2, mitot-1

          ux = .5d0*(q(2,i,j)-q(2,i-1,j)+q(2,i,j-1)-q(2,i-1,j-1))
          vy = .5d0*(q(3,i,j)-q(3,i,j-1)+q(3,i-1,j)-q(3,i-1,j-1))

          div(i,j) = qm*dmax1(0.d0,-(ux + vy))
          if (dabs(div(i,j)) .lt. 1.e-25) div(i,j) = 0.d0
 20    continue

c
       return
       end
c
c ----------------------------------------------------------------------------
c
        subroutine getFaceMidpts(xfaceMidpt,yfaceMidpt,irr,
     &                   dx,dy,xlow,ylow,lstgrd,mitot,mjtot)

      use amr_module
      implicit double precision (a-h, o-z)

      dimension xfaceMidpt(2,mitot+1,mjtot+1)
      dimension yfaceMidpt(2,mitot+1,mjtot+1)
      dimension irr(mitot,mjtot)

c     compute midpoints (x,y) for vertical (f or x-faces) and
c     horizontal (g or y-faces) for the whole grid.

c     first do as if everything regular, then overwrite for cut cells
c     solid cells not important and will leave
c     index 1 is for x, 2 for y

      do j = 1, mjtot
      do i = 1, mitot+1
         xfaceMidpt(1,i,j) = xlow + (i-1)*dx
         xfaceMidpt(2,i,j) = ylow + (j-0.5d0)*dy
      end do
      end do

      do j = 1, mjtot+1
      do i = 1, mitot
         yfaceMidpt(1,i,j) = xlow + (i-0.5d0)*dx
         yfaceMidpt(2,i,j) = ylow + (j-1)*dy
      end do
      end do

c     now go around for each cut cell
      do j = 1, mjtot
      do i = 1, mitot
         k = irr(i,j)
         if (k .eq. lstgrd .or. k .eq. -1) cycle
c        
         do kside = 1,6
            if (poly(kside+2,1,k) .eq. -11) exit
            x1 = poly(kside,1,k)
            y1 = poly(kside,2,k)
            x2 = poly(kside+1,1,k)
            y2 = poly(kside+1,2,k)
            if (y1 .eq. y2) then !  y face
               ixadj = i 
               iyadj = j + isig(x2-x1)
               if (iyadj .lt. j) then ! set left edge owned by cell i,j
                  yfaceMidpt(1,i,j) = 0.5d0*(x1+x2)
               else               
                yfaceMidpt(1,i,iyadj) = 0.5d0*(x1+x2)
               endif
               
            else if (x1 .eq. x2) then  ! x face
               ixadj = i + isig(y1-y2)
               iyadj = j
               if (ixadj .lt. i) then ! set left side owned by cell ij
                  xfaceMidpt(2,i,j) = 0.5d0*(y1+y2)
               else
                  xfaceMidpt(2,ixadj,j) = 0.5d0*(y1+y2)
               endif   
            else
              write(*,*)"screwy face for cell ",i,j
            endif
         end do
      end do
      end do

      return
      end
c
c -------------------------------------------------------------------
c
      subroutine convert2cons(corr,nvar)

      implicit double precision (a-h,o-z)

      include "cuserdt.i"
      dimension corr(nvar)

       rho = corr(1)
       if (rho .eq. 0.d0) return  !or assert all entries 0?
       u = corr(2)
       v = corr(3)
       p = corr(4)
       ei = p/(rho*gamma1)
       en = ei + .5d0*(u*u+v*v)
       corr(2) = corr(2) * rho
       corr(3) = corr(3) * rho
       corr(4) = en*rho

       return
       end
c
c -------------------------------------------------------------------
c
      subroutine convert2prim(state,nvar)

      implicit double precision (a-h,o-z)

       include "cuserdt.i"
       dimension state(nvar)

       rho = state(1)
       if (rho .eq. 0.d0) return  !or assert all entries 0?
       u = state(2)/rho
       v = state(3)/rho
       p = gamma1*(state(4)-0.5d0*rho*(u*u+v*v))
       state(2) = u
       state(3) = v
       state(4) = p

       return
       end
c
c ------------------------------------------------------------------
c
       subroutine divcorn(q,qold,div,mitot,mjtot,dx,dy,mptr,nvar,lstgrd,
     .                    xlow,ylow,irr,nghost,qx,qy,
     .                    alldivf,alldivg)

       use amr_module, only:  xcirr,ycirr
       implicit double precision (a-h,o-z)
       dimension q(nvar,mitot,mjtot), qold(nvar,mitot,mjtot)
       dimension div(mitot,mjtot)
       dimension alldivf(nvar,mitot,mjtot), alldivg(nvar,mitot,mjtot)
       dimension qx(nvar,mitot,mjtot),qy(nvar,mitot,mjtot)
       
       dimension irr(mitot,mjtot)
       dimension qplus(nvar),qminus(nvar)
c      dimension u(mitot,mjtot),v(mitot,mjtot)
       logical ALL4REG, ALL5REG, all9reg

c      parameter (qm=.2d0)
       parameter (qm=.30d0)
c      parameter (qm=.250d0)

c :::::::::::::::::::::: DIVCORN :::::::::::::::::::::::::::::::::::::;1438
c  compute divergence to use in artificial viscosity in method
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

c  q already in primitive variables
c  qold in conserved variables

       div     = 0.d0
       alldivf = 0.d0
       alldivg = 0.d0
c
c :::  cell i,j owns lower left corner divergence
c
       do 20 j = 2, mjtot-1
       do 20 i = 2, mitot-1

        if (irr(i,j) .eq. -1) cycle

        if (irr(i,j) .eq. lstgrd) then 
c          at least we know  nboring cells exist
           ux = q(2,i+1,j) - q(2,i-1,j)
           vy = q(3,i,j+1) - q(3,i,j-1)
           quse = qm
           div(i,j) = quse*dmax1(0.d0,-(ux + vy))
           if (dabs(div(i,j)) .lt. 1.e-25) div(i,j) = 0.d0
        endif

 20    continue

       do 30 j = 2, mjtot-1
       do 30 i = 2, mitot-1
          if (irr(i,j) .eq. -1) cycle

          coeff = div(i,j)
          do joff = -1, 1
          do ioff = -1, 1
             coeff = max(coeff,div(i+ioff,j+joff))
          end do
          end do
             alldivf(:,i,j) = coeff*(qold(:,i,j)-qold(:,i-1,j))
c          endif
         
c          if (irr(i,j-1) .ne. -1) then
             alldivg(:,i,j) = coeff*(qold(:,i,j) - qold(:,i,j-1))
c          endif

 30    continue

c
       return
       end
c
c -----------------------------------------------------------------------
c
       subroutine getMidpt(i,j,x,y,xlow,ylow,dx,dy,irr,lstgrd,
     &                     mitot,mjtot)

        use amr_module, only:  xcirr,ycirr
       implicit double precision (a-h,o-z)
       dimension irr(mitot,mjtot)

       k = irr(i,j)
       if (k .eq. lstgrd) then
         x = xlow + (i-0.5d0)*dx
         y = ylow + (j-0.5d0)*dy
      else if (k .eq. -1) then
        write(*,*)"should not find solid cell in getMdpt; index",i,j
        stop
      else
       x = xcirr(k)
       y = ycirr(k)
      endif

      return
      end
       
c
c ----------------------------------------------------------------------
c
       subroutine xfacerecenter(qrec,q,i,j,qx,qy,lstgrd,xlow,ylow,dx,dy,
     &                     irr,mitot,mjtot,xface,yface,nvar)

       use amr_module, only:  xcirr,ycirr
       implicit double precision (a-h,o-z)
       dimension irr(mitot,mjtot),qrec(nvar)
       dimension q(nvar,mitot,mjtot)
       dimension qx(nvar,mitot,mjtot),qy(nvar,mitot,mjtot)

       k = irr(i,j)
       if (k .eq. lstgrd) then
         xcell = xlow + (i-0.5d0)*dx
         ycell = ylow + (j-0.5d0)*dy
      else if (k .eq. -1) then
        write(*,*)"should not find solid cell in getMdpt; index",i,j
        stop
      else
       xcell = xcirr(k)
       ycell = ycirr(k)
      endif

c     recenter to y face val of x face so can difference
      qrec = q(:,i,j) + (yface-ycell)*qy(:,i,j)

      return  
      end
c
c ----------------------------------------------------------------------
c
       subroutine yfacerecenter(qrec,q,i,j,qx,qy,lstgrd,xlow,ylow,dx,dy,
     &                     irr,mitot,mjtot,xface,yface,nvar)

       use amr_module, only:  xcirr,ycirr
       implicit double precision (a-h,o-z)
       dimension irr(mitot,mjtot),qrec(nvar)
       dimension q(nvar,mitot,mjtot)
       dimension qx(nvar,mitot,mjtot),qy(nvar,mitot,mjtot)

       k = irr(i,j)
       if (k .eq. lstgrd) then
         xcell = xlow + (i-0.5d0)*dx
         ycell = ylow + (j-0.5d0)*dy
      else if (k .eq. -1) then
        write(*,*)"should not find solid cell in getMdpt; index",i,j
        stop
      else
       xcell = xcirr(k)
       ycell = ycirr(k)
      endif

c     recenter to y face val of x face so can difference
      qrec = q(:,i,j) + (xface-xcell)*qx(:,i,j)

      return  
      end
