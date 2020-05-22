c
!> Integrate all grids at the input **level** by one step of its delta(t)
!!
!! this includes:  
!! - setting the ghost cells 
!! - advancing the solution on the grid
!! - adjusting fluxes for flux conservation step later
c --------------------------------------------------------------
c
      subroutine advanc (level,nvar,dtlevnew,vtime,naux)
c
      use amr_module
      implicit double precision (a-h,o-z)
      common /RKmethod/ coeff(5),mstage


      logical    vtime
      integer omp_get_thread_num, omp_get_max_threads
      integer mythread/0/, maxthreads/1/
      integer listgrids(numgrids(level))
      integer clock_start, clock_finish, clock_rate
      integer clock_startStepgrid,clock_startBound,clock_finishBound
      real(kind=8) cpu_start, cpu_finish
      real(kind=8) cpu_startBound, cpu_finishBound
      real(kind=8) cpu_startStepgrid, cpu_finishStepgrid

c     maxgr is maximum number of grids  many things are
c     dimensioned at, so this is overall. only 1d array
c     though so should suffice. problem is
c     not being able to dimension at maxthreads


c
c  ::::::::::::::; ADVANC :::::::::::::::::::::::::::::::::::::::::::
c  integrate all grids at the input  'level' by one step of its delta(t)
c  this includes:  setting the ghost cells 
c                  advancing the solution on the grid
c                  adjusting fluxes for flux conservation step later
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c get start time for more detailed timing by level
      call system_clock(clock_start,clock_rate)
      call cpu_time(cpu_start)

c  mjb adapted june 8, 2019 for multistage RK methods.
c  loop over all grids, each stage, and copy ghost cells
c  to avoid needing so many
c
c  mjb adapted May 18, 2020 to do SRD call after all grdis updated and new
c ghost cell values copied in

      do istage = 1, mstage

      hx   = hxposs(level)
      hy   = hyposs(level)
      delt = possk(level)
c
      call system_clock(clock_startBound,clock_rate)
      call cpu_time(cpu_startBound)


c     maxthreads initialized to 1 above in case no openmp
!$    maxthreads = omp_get_max_threads()

c We want to do this regardless of the threading type
!$OMP PARALLEL DO PRIVATE(j,locnew, locaux, mptr,nx,ny,mitot,
!$OMP&                    mjtot,time,levSt),
!$OMP&            SHARED(level, nvar, naux, alloc, intrat, delt,
!$OMP&                   listOfGrids,listStart,nghost,
!$OMP&                   node,rnode,numgrids,listgrids,istage),
!$OMP&            SCHEDULE (dynamic,1)
!$OMP&            DEFAULT(none)
      do j = 1, numgrids(level)
         levSt = listStart(level)
         mptr   = listOfGrids(levSt+j-1)
         mitot  = node(ndihi,mptr) - node(ndilo,mptr) + 1 + 2*nghost
         mjtot  =  node(ndjhi,mptr) - node(ndjlo,mptr) + 1 + 2*nghost
         locnew = node(store1,mptr)
         locaux = node(storeaux,mptr)
         time   = rnode(timemult,mptr)
c     
          call bound(time,nvar,nghost,alloc(locnew),mitot,mjtot,mptr,
     1               alloc(locaux),naux)
       end do
!$OMP END PARALLEL DO
      call system_clock(clock_finishBound,clock_rate)
      call cpu_time(cpu_finishBound)
      timeBound = timeBound + clock_finishBound - clock_startBound
      timeBoundCPU=timeBoundCPU+cpu_finishBound-cpu_startBound
      
c
      dtlevnew = rinfinity
c     cfl_level = 0.d0    !# to keep track of max cfl seen on each level

c 
      call system_clock(clock_startStepgrid,clock_rate)
      call cpu_time(cpu_startStepgrid)


!$OMP PARALLEL DO PRIVATE(j,mptr,nx,ny,mitot,mjtot)  
!$OMP&            PRIVATE(mythread,dtnew)
!$OMP&            SHARED(rvol,rvoll,level,nvar,mxnest,alloc,intrat)
!$OMP&            SHARED(nghost,intratx,intraty,hx,hy,naux,listsp)
!$OMP&            SHARED(node,rnode,dtlevnew,numgrids,listgrids)
!$OMP&            SHARED(istage,mstage,ar,ixg,iyg,nxtirr)
!$OMP&            SHARED(listOfGrids,listStart,levSt,vtime,delt)
!$OMP&            SCHEDULE (DYNAMIC,1)
!$OMP&            DEFAULT(none)
      do j = 1, numgrids(level)
         levSt = listStart(level)
         mptr = listOfGrids(levSt+j-1)
         mitot  = node(ndihi,mptr) - node(ndilo,mptr) + 1 + 2*nghost
         mjtot  = node(ndjhi,mptr) - node(ndjlo,mptr) + 1 + 2*nghost
c
         call par_advanc(mptr,mitot,mjtot,nvar,naux,dtnew,vtime,
     &                   istage,mstage)

c        update for oneeach stage so that bc's set properly
c        for next  stage. will take it off when later
c
         rnode(timemult,mptr) = rnode(timemult,mptr)+delt
      end do
!$OMP END PARALLEL DO
c
c     next set of boundary conditions
      call system_clock(clock_startBound,clock_rate)
      call cpu_time(cpu_startBound)
c redo ghost cells (one level only, no refinement) before SRD
c
!$OMP PARALLEL DO PRIVATE(j,locnew, locaux, mptr,nx,ny,mitot,
!$OMP&                    mjtot,time,levSt),
!$OMP&            SHARED(level, nvar, naux, alloc, intrat, delt,
!$OMP&                   listOfGrids,listStart,nghost,
!$OMP&                   node,rnode,numgrids,listgrids,istage),
!$OMP&            SCHEDULE (dynamic,1)
!$OMP&            DEFAULT(none)
      do j = 1, numgrids(level)
         levSt = listStart(level)
         mptr   = listOfGrids(levSt+j-1)
         mitot  = node(ndihi,mptr) - node(ndilo,mptr) + 1 + 2*nghost
         mjtot  = node(ndjhi,mptr) - node(ndjlo,mptr) + 1 + 2*nghost
         locnew = node(store1,mptr)
         locaux = node(storeaux,mptr)
         time   = rnode(timemult,mptr)
c     
          call bound(time,nvar,nghost,alloc(locnew),mitot,mjtot,mptr,
     1               alloc(locaux),naux)

       end do
!$OMP END PARALLEL DO

      call system_clock(clock_finishBound,clock_rate)
      call cpu_time(cpu_finishBound)
      timeBound = timeBound + clock_finishBound - clock_startBound
      timeBoundCPU=timeBoundCPU+cpu_finishBound-cpu_startBound

!$OMP PARALLEL DO PRIVATE(j,locnew, locaux, mptr,nx,ny,lstgrd,xlow,ylow,
!$OMP&                    mitot,mjtot,time,levSt,locirr,locncount,
!$OMP&                    locnumHoods),
!$OMP&            SHARED(level, nvar, naux, alloc, intrat, delt,dtnew,
!$OMP&                   listOfGrids,listStart,nghost,istage,mstage,
!$OMP&                   dtlevnew,node,rnode,numgrids,listgrids,hx,hy),
!$OMP&            SCHEDULE (dynamic,1)
!$OMP&            DEFAULT(none)
      do j = 1, numgrids(level)
         levSt = listStart(level)
         mptr   = listOfGrids(levSt+j-1)
         mitot  = node(ndihi,mptr) - node(ndilo,mptr) + 1 + 2*nghost
         mjtot  = node(ndjhi,mptr) - node(ndjlo,mptr) + 1 + 2*nghost
         locnew = node(store1,mptr)
         locaux = node(storeaux,mptr)
         time   = rnode(timemult,mptr)
         locirr = node(permstore,mptr)
         locncount = locirr + mitot*mjtot
         locnumHoods = locncount + mitot*mjtot
         lstgrd = node(lstptr,mptr)
         xlow = rnode(cornxlo,mptr) - nghost*hx
         ylow = rnode(cornylo,mptr) - nghost*hy
         call srd_cellMerge(alloc(locnew),nvar,alloc(locirr),mitot,
     &                     mjtot,lstgrd,hx,hy,nghost,xlow,ylow,istage,
     &                     alloc(locncount),alloc(locnumHoods),mptr)
         if (istage .eq. mstage) then ! set new time step
            call estdt(alloc(locnew),alloc(locirr),mitot,mjtot,nvar,
     &                 hx,hy,dtnew,nghost,alloc(locaux),naux)
!$OMP CRITICAL (newdt)
          dtlevnew = dmin1(dtlevnew,dtnew)
!$OMP END CRITICAL (newdt)    
            ! grid has been updated in time for each stage. subtract off so
            ! only one update. Done this way to be able to tell that grids
            ! updated in case of refinement, or
            rnode(timemult,mptr)  = rnode(timemult,mptr)-(mstage-1)*delt
         endif
      end do
!$OMP END PARALLEL DO

      end do  ! end loop over each stage

c   --------- done with everything at this level --------------
c     final clock updates
      call system_clock(clock_finish,clock_rate)
      call cpu_time(cpu_finish)
      ! both stages count towards number of cell updates
      tvoll(level) = tvoll(level) + clock_finish - clock_start
      tvollCPU(level) = tvollCPU(level) + cpu_finish - cpu_start
      timeStepgrid = timeStepgrid +clock_finish-clock_startStepgrid
      timeStepgridCPU=timeStepgridCPU+cpu_finish-cpu_startStepgrid      
c

      return
      end
c
c --------------------------------------------------------------
c
!> Integrate grid **mptr**. grids are done in parallel.
      subroutine par_advanc (mptr,mitot,mjtot,nvar,naux,dtnew,
     .                       vtime,istage,mstage)
c
      use amr_module
      use gauges_module, only: update_gauges, num_gauges
      implicit double precision (a-h,o-z)


      integer omp_get_thread_num, omp_get_max_threads
      integer mythread/0/, maxthreads/1/

c     double precision fp(nvar,mitot,mjtot),fm(nvar,mitot,mjtot)
c     double precision gp(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)

      logical vtime


c
c  :::::::::::::: PAR_ADVANC :::::::::::::::::::::::::::::::::::::::::::
c  integrate this grid. grids are done in parallel.
c  extra subr. used to allow for stack based allocation of
c  flux arrays. They are only needed temporarily. If used alloc
c  array for them it has too long a lendim, makes too big
c  a checkpoint file, and is a big critical section.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      level = node(nestlevel,mptr)
      hx    = hxposs(level)
      hy    = hyposs(level)
      delt  = possk(level)
      nx    = node(ndihi,mptr) - node(ndilo,mptr) + 1
      ny    = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
      time  = rnode(timemult,mptr)

!$    mythread = omp_get_thread_num()

      locold = node(store2, mptr)
      locnew = node(store1, mptr)

c
c  copy old soln. values into  next time step's soln. values
c  since integrator will overwrite it. only for grids not at
c  the finest level. finest level grids do not maintain copies
c  of old and new time solution values.
c
         ! do it all the time now, needed for multistage RK
         ! even if no refinement.  But only do it on
         ! first stage, other locold holds qold for RK update
         !if (level .lt. mxnest) then
         if (istage .eq. 1) then
             ntot   = mitot * mjtot * nvar
cdir$ ivdep
             do 10 i = 1, ntot
 10            alloc(locold + i - 1) = alloc(locnew + i - 1)
         endif
c
         xlow = rnode(cornxlo,mptr) - nghost*hx
         ylow = rnode(cornylo,mptr) - nghost*hy

!$OMP CRITICAL(rv)
      rvol = rvol + nx * ny
      rvoll(level) = rvoll(level) + nx * ny
!$OMP END CRITICAL(rv)


         locaux = node(storeaux,mptr)
         locirr = node(permstore,mptr)
         locncount = locirr + mitot*mjtot
         locnumHoods = locncount + mitot*mjtot
c
         if (node(ffluxptr,mptr) .ne. 0) then
            lenbc  = 2*(nx/intratx(level-1)+ny/intraty(level-1))
            locsvf = node(ffluxptr,mptr)
            locsvq = locsvf + nvar*lenbc
            locx1d = locsvq + nvar*lenbc
         endif

c        # See if the grid about to be advanced has gauge data to output.
c        # This corresponds to previous time step, but output done
c        # now to make linear interpolation easier, since grid
c        # now has boundary conditions filled in.

c     should change the way print_gauges does io - right now is critical section
c     no more,  each gauge has own array.

      if (num_gauges > 0) then
           call update_gauges(alloc(locnew:locnew+nvar*mitot*mjtot),
     .                       alloc(locaux:locaux+nvar*mitot*mjtot),
     .                       xlow,ylow,nvar,mitot,mjtot,naux,mptr)
         endif

c
         if (dimensional_split .eq. 0) then
c           # Unsplit method
c        call stepgrid(alloc(locnew),alloc(locold),fm,fp,gm,gp,
c    2                 mitot,mjtot,nghost,
c    3                 delt,dtnew,hx,hy,nvar,
c    4                 xlow,ylow,time,mptr,naux,alloc(locaux),
c    5                 alloc(locirr),node(lstptr,mptr),
c    6                 alloc(locncount),alloc(locnumHoods),vtime,istage)
         call mymethod(alloc(locnew),alloc(locold),mitot,mjtot,nghost,
     1                 delt,dtnew,hx,hy,nvar,xlow,ylow,mptr,naux,
     2                 alloc(locaux),alloc(locirr),node(lstptr,mptr),
     3                 alloc(locncount),alloc(locnumHoods),vtime,
     4                 istage,time)
         else if (dimensional_split .eq. 1) then
c           # Godunov splitting
            write(6,*)"this option not supported"
            stop
         else 
c           # should never get here due to check in amr2
            write(6,*) '*** Strang splitting not supported'
            stop
         endif

c     if (node(cfluxptr,mptr) .ne. 0)
c    2   call fluxsv(mptr,fm,fp,gm,gp,
c    3               alloc(node(cfluxptr,mptr)),mitot,mjtot,
c    4               nvar,listsp(level),delt,hx,hy)
c        if (node(ffluxptr,mptr) .ne. 0) then
c        lenbc = 2*(nx/intratx(level-1)+ny/intraty(level-1))
c           locsvf = node(ffluxptr,mptr)
c        call fluxad(fm,fp,gm,gp,
c    2               alloc(locsvf),mptr,mitot,mjtot,nvar,
c    4                  lenbc,intratx(level-1),intraty(level-1),
c    5               nghost,delt,hx,hy)
c        endif
c
c        write(outunit,969) mythread,delt, dtnew
c969     format(" thread ",i4," updated by ",e15.7, " new dt ",e15.7)
         ! note above that for the second stage time incremetned
c
      return
      end
