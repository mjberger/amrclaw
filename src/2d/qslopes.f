c
c ---------------------------------------------------------------------
c
       subroutine qslopes(qp,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,
     &                    hx,hy,xlow,ylow,mptr,nvar) 
      use amr_module
      implicit double precision(a-h,o-z)
      include "cuserdt.i"

      dimension qp(nvar,mitot,mjtot),qx(nvar,mitot,mjtot),
     &          qy(nvar,mitot,mjtot),irr(mitot,mjtot)

      ! driver routine to call the correct gradient routine
      ! (use for both cell gradients and merge nhood gradients)
      ! igradChoice =     
      !  0 = no gradient
      !  1 = first order accurate gradients  (2 terms in least squares)
      !  2 = pointwise quadratic gradients (5 terms in least squares)
      !  3 = cellwise average  quadratic gradients (5 terms in least squares)


      if (igradChoice .eq. 1 .or. igradChoice .eq. 2) then
         call lp_qslopesPtQuad2(qp,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,
     &                   hx,hy,xlow,ylow,mptr,nvar)
      else if (igradChoice .eq. 3) then
         call lp_qslopesWithGhostAvgQuadratic(qp,qx,qy,mitot,mjtot,irr,
     &                   lstgrd,lwidth,hx,hy,xlow,ylow,mptr,nvar)
      else if (igradChoice .eq. 0) then
         write(*,*)"should not be here  should have set ssw = 0"
      else
         write(*,*)"unrecognized gradient choice",igradsChoice
         write(*,*)"should test in setprob and stop there"
         stop
      endif

      return
      end
