      subroutine setprob

      use amr_module
      implicit real*8 (a-h,o-z)
      character(len=25) fname
c
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)

      read(7,*) nloops
      write(*,*) " This geometry has ", nloops," loops"

      do n= 1, nloops
        read(7,*) xloops(n)
        read(7,*) yloops(n)
      end do


      iprob = 16
      gamma = 1.4d0
      gamma1 = gamma - 1.d0
      xprob = xupper
      yprob = yupper

      write(*,*)"setprob iprob ",iprob
      return
      end
