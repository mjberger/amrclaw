      subroutine setprob

      use amr_module
      implicit real*8 (a-h,o-z)
      character(len=25) fname
c
      logical pwconst
      common /RKmethod/ coeff(5), mstage
      common /userdt/ cflcart,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold,pwconst
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)

      read(7,*) mstage
      write(*,*) "will use RK scheme with ", mstage," stages"

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

      read(7,*) pwconst
      write(*,*)"Plotting output with pwconst = ",pwconst

      read(7,*) max1d
      write(*,*)"Using grid patches of size  max1d = ",max1d

      read(7,*) nloops
      write(*,*) "This geometry has ", nloops," loops"

      do n= 1, nloops
        read(7,*) xloops(n)
        read(7,*) yloops(n)
      end do


      iprob = 21
      write(*,*)"Setprob is setting iprob = ",iprob

      gamma = 1.4d0
      gamma1 = gamma - 1.d0
      xprob = xupper
      yprob = yupper
      cflcart = cfl  ! have to go through and only use one


      return
      end
