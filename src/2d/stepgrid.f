c
c -------------------------------------------------------------
c
!> Take a time step on a single grid **mptr** and overwrite solution array **q**. 
!! A modified version of the clawpack routine step2 is used.
!! 
!! Return new solution **q** as well as fluxes in fm,fp and gm,gp.
!! Patch has room for ghost cells (mbc of them) around the grid.
!! Everything is the enlarged size (**mitot** by **mjtot**).
!! 
!! \param[in] mbc number of ghost cells  (= lwidth)
!! \param[in] mptr grid number (for debugging)
!! \param[in] xlow left edge of enlarged grid (including ghost cells).
!! \param[in] ylow lower edge of enlarged grid (including ghost cells).
!! \param[in] dt incoming time step
!! \param[in] dx mesh size in x direction for this grid
!! \param[in] dx mesh size in y direction for this grid
!! \param[in] irr irregular cell array
!! \param[in] lstgrd  indicates which cell index is regular on this grid
!! \param[in,out] q solution array
!! \param[out] dtnew  return suggested new time step for this grid's soln.
!! \param[out] fm fluxes on the left side of each vertical edge
!! \param[out] fp fluxes on the right side of each vertical edge
!! \param[out] gm fluxes on the lower side of each horizontal edge
!! \param[out] gp fluxes on the upper side of each horizontal edge
      subroutine stepgrid(q,fm,fp,gm,gp,mitot,mjtot,mbc,dt,dtnew,dx,dy,
     &                    nvar,xlow,ylow,time,mptr,maux,aux,irr,lstgrd,
     &                    ncount,numHoods,vtime)
c
c          
c ::::::::::::::::::: STEPGRID   ::::::::::::::::::::::::::::::::::::
c  was stepgrid, now in a  for cut cell version
c  calls method ( not step2 ) to do the real work
c take a time step on a single grid. overwrite solution array q. 
c A modified version of the clawpack routine step2 is used.
c
c return fluxes in fm,fp and gm,gp.
c patch has room for ghost cells (mbc of them) around the grid.
c everything is the enlarged size (mitot by mjtot).
c
c mbc       = number of ghost cells  (= lwidth)
c mptr      = grid number  (for debugging)
c xlow,ylow = lower left corner of enlarged grid (including ghost cells).
c dt         = incoming time step
c dx,dy      = mesh widths for this grid
c dtnew      = return suggested new time step for this grid's soln.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      use amr_module
      implicit double precision (a-h,o-z)
      external rpn2,rpt2

      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

      dimension q(nvar,mitot,mjtot)
      dimension fp(nvar,mitot,mjtot),gp(nvar,mitot,mjtot)
      dimension fm(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)
      dimension aux(maux,mitot,mjtot)
      integer   irr(mitot,mjtot), ncount(mitot,mjtot)
      integer   numHoods(mitot,mjtot)


      logical    debug,  dump, vtime
      data       debug/.false./,  dump/.false./

c
c     # set tcom = time.  This is in the common block comxyt that could
c     # be included in the Riemann solver, for example, if t is explicitly
c     # needed there.

      write(*,*)"stepgrid 0 vtime ",vtime
      tcom = time

      if (dump) then
         write(outunit,*) "dumping grid ",mptr," at time ",time
         do i = 1, mitot
         do j = 1, mjtot
            write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar) 
c    .                  ,(aux(ivar,i,j),ivar=1,maux)
 545        format(2i4,5e15.7)
         end do
         end do
      endif
c
      meqn   = nvar
      mx = mitot - 2*mbc
      my = mjtot - 2*mbc
      maxm = max(mx,my)       !# size for 1d scratch array
      mbig = maxm
      xlowmbc = xlow + mbc*dx
      ylowmbc = ylow + mbc*dy

c     # method(2:7) and mthlim
c     #    are set in the amr2ez file (read by amr)
c
      method(1) = 0
c
c  in case anything need adjusting in aux arrays
      call b4step2(mbc,mx,my,nvar,q,
     &             xlowmbc,ylowmbc,dx,dy,time,dt,maux,aux)
c
c
c     # take one step (2 stages)  on the conservation law:
c
      write(*,*)"stepgrid 1 vtime ",vtime
      call mymethod(q,fm,fp,gm,gp,mitot,mjtot,mbc,dt,dtnew,dx,dy,
     &            nvar,xlow,ylow,mptr,maux,aux,irr,lstgrd,
     &            ncount,numHoods,vtime)
      write(*,*)"stepgrid after mymethod dt,dtnew ",dt,dtnew
c
c
c!$OMP  CRITICAL (cflm)

c        cfl_level = dmax1(cfl_level,cflgrid)

c!$OMP END CRITICAL (cflm)

c
c
      if (method(5).eq.1) then
c        # with source term:   use Godunov splitting
         call src2(nvar,mbc,mx,my,xlowmbc,ylowmbc,dx,dy,
     &             q,maux,aux,time,dt)
         endif
c
c
c
c     # output fluxes for debugging purposes:
      if (debug) then
         write(dbugunit,*)" fluxes for grid ",mptr
c        do 830 j = mbc+1, mjtot-1
            do 830 i = mbc+1, mitot-1
         do 830 j = mbc+1, mjtot-1
               write(dbugunit,831) i,j,fm(1,i,j),fp(1,i,j),
     .                             gm(1,i,j),gp(1,i,j)
               do 830 m = 2, meqn
                  write(dbugunit,832) fm(m,i,j),fp(m,i,j),
     .            gm(m,i,j),gp(m,i,j)
  831          format(2i4,4d16.6)
  832          format(8x,4d16.6)
  830    continue
      endif

c
c
c For variable time stepping, use max speed seen on this grid to 
c choose the allowable new time step dtnew.  This will later be 
c compared to values seen on other grids.
c
c      cflgrid = cfl  ! not really should fix this but just want dtnew
c      if (cflgrid .gt. 0.d0) then
c          dtnew = dt*cfl/cflgrid
c        else
c          # velocities are all zero on this grid so there's no 
c          # time step restriction coming from this grid.
c           dtnew = rinfinity
c         endif

c     # give a warning if Courant number too large...
c
c     if (cflgrid .gt. cflv1) then
c           write(*,810) cflgrid
c           write(outunit,810) cflgrid, cflv1
c 810       format('*** WARNING *** Courant number  =', d12.4,
c    &              '  is larger than input cfl_max = ', d12.4)
c           endif
c
      if (dump) then
         write(outunit,*) "dumping grid ",mptr," after stepgrid"
         do i = 1, mitot
         do j = 1, mjtot
            write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar)
         end do
         end do
      endif
      write(*,*)"stepgrid at end dt,dtnew ",dt,dtnew
      return
      end


