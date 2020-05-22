      subroutine setprob

      use amr_module
      implicit real*8 (a-h,o-z)
      character(len=25) fname
c
      logical nolimiter, quad
      common /RKmethod/ coeff(5), mstage
      common /order2/ ssw, quad, nolimiter
      include "cuserdt.i"
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)

      read(7,*) mstage
      write(*,*) "will use evolution scheme with ", mstage," stages"

      coeff = 0.d0 
      if (mstage .eq. 1) then
         coeff(1) = 1.d0
      else if (mstage .eq. 2) then
         coeff(1) = 0.5d0
         coeff(2) = 1.0d0
      endif


      ! check for higher order if the structure in my method holds
      ! or need more scratch storage etc
      read(7,*) ismp
      write(*,*)"Using stabilization ismp = ",ismp
      write(outunit,*)"Using stabilization ismp = ",ismp

      read(7,*) pwconst
      write(*,*)"Plotting output with pwconst = ",pwconst
      write(outunit,*)"Plotting output with pwconst = ",pwconst

      read(7,*) max1d
      write(*,*)"Using grid patches of size  max1d = ",max1d
      write(outunit,*)"Using grid patches of size  max1d = ",max1d

      read(7,*) nloops
      write(*,*) "This geometry has ", nloops," loops"
      write(outunit,*) "This geometry has ", nloops," loops"

      do n= 1, nloops
        read(7,*) xloops(n)
        read(7,*) yloops(n)
      end do

      read(7,*) ghost_ccg
      write(*,*)"Use of ghost cells in gradients = ",ghost_ccg
      write(outunit,*)"Use of ghost cells in gradients = ",ghost_ccg

      read(7,*) limitTile
      write(*,*)"Tile limiter: 1 (BJ), 2 (LP) = ",limitTile
      write(outunit,*)"Tile limiter: 1 (BJ), 2 (LP) = ",limitTile

      read(7,*) lpChoice
      write(*,*)"lpChoice (if used): 1 (standard), 2 (relaxed) = ",
     &           lpChoice
      write(outunit,*)"lpChoice (if used): 1 (standard), 2 (relaxed)= ",
     &           lpChoice


      !read(7,*) nTerms
      !write(*,*) "Use ",nTerms," terms in cell gradient"
      !write(outunit,*) "Use ", nTerms," terms in cell gradient"

      !read(7,*) numMergeTerms
      !write(*,*) "Use ",numMergeTerms," terms in tile gradient"
      !write(outunit,*) "Use ", numMergeTerms," terms in tile gradient"
      !write(*,*)" 2 terms = 1st order gradient, 5 = 2nd order"
      !write(outunit,*)" 2 terms = 1st order gradient, 5 = 2nd order"

      read(7,*) igradChoice
      write(*,*)" Use ",igradCoice," choice for gradient"
      write(*,*)" 0=none, 1=1st order, 2=ptwise quad.,3=cell avg.quad."
      write(outunit,*)" Use ",igradCoice," choice for gradient"
      write(outunit,*)" 0 = none, 1=1st order, 2=ptwise quad., ",
     &                " 3 = cell avg.quad."
      if (igradChoice .eq. 0) then
         if (ssw .ne. 0) then
            write(*,*)" ssw = ",ssw," resetting to 0 to match ",
     &                " igradChoice =  ", igradChoice
            ssw = 0.d0
            nolimiter = .true. ! might as well since no gradients
         endif
      endif

      iprob = 19
      write(*,*)"Setprob is setting iprob = ",iprob
      write(outunit,*)"Setprob is setting iprob = ",iprob

      gamma = 1.4d0
      gamma1 = gamma - 1.d0
      xprob = xupper
      yprob = yupper
      cflcart = cfl  ! have to go through and only use one
      gradThreshold = 1.d-4


      return
      end
