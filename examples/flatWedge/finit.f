c
c -------------------------------------------------------------
c
        subroutine finit(val,nvar,mitot,mjtot,irr,
     1               cornx,corny,hx,hy,lstgrd)
c
       implicit double precision (a-h,o-z)
       dimension val(mitot,mjtot,nvar), irr(mitot,mjtot)

      include  "cnodal.i"
      include  "cuser.i"
      include  "cdisc.i"
      include  "cdom.i"

      if (iprob .eq. 1) then
c     # shock/cloud initial conditions
c
        xshock = 2.*xupper 
 	do 30 i = 1, mitot
 	   xcen = cornx + (i-.5)*hx
 	   if (xcen .le. xshock) then
 		do 10 j = 1, mjtot
 		      val(i,j,1) = rhoshk
 		      val(i,j,2) = rhoshk*ushk
 		      val(i,j,3) = 0.d0
 		      val(i,j,4) = rhoshk*(eshk+.5d0*ushk**2)
 10             continue
            else
 		do 20 j = 1, mjtot
 		   ycen = corny + (j-.5)*hy
c		   dist = sqrt((xcen-xcircle)**2 + (ycen-ycircle)**2)
c		   if (dist - radius .ge. 0.d0) then
 		      rhotmp = rhoamb
 		      utmp   = uamb
 		      etmp   = eamb
c		   else
c		      rhotmp = rhocld
c		      utmp   = ucld
c 		      etmp   = ecld
c		   endif
 		   val(i,j,1) = rhotmp
 		   val(i,j,2) = rhotmp*utmp
 		   val(i,j,3) = 0.d0
 		   val(i,j,4) = rhotmp*(etmp + .5*utmp**2)
 
 20              continue
             endif
 30     continue

      else if (iprob .eq. 3) then

	do 130 i = 1, mitot
	   xcen = cornx + (i-.5)*hx
	   if (xcen .gt. xshock) then
	      do 110 j = 1, mjtot
		val(i,j,1) = rhoshk
		val(i,j,2) = rhoshk*ushk
		val(i,j,3) = 0.d0
		val(i,j,4) = rhoshk*(eshk+.5d0*ushk**2)
 110          continue
           else
	      do 120 j = 1, mjtot
		 ycen = corny + (j-.5)*hy
		 dist = dsqrt((xcen-xcircle)**2 + (ycen-ycircle)**2)
		 theta = datan2(ycen-ycircle,xcen-xcircle)
		 upert = eps*(dist/rc)*dexp(alpha*(1.d0-(dist/rc)**2))
     .                    *dsin(theta)
		 vpert =-eps*(dist/rc)*dexp(alpha*(1.d0-(dist/rc)**2))
     .                    *dcos(theta)
		 tpert = -((gamma-1.d0)*eps*eps/(4.d0*alpha*gamma))*
     .                   dexp(2.d0*alpha*(1.d0-(dist/rc)**2))
		 spert = 0.d0

		 ttmp    = tamb+tpert
		 rhotmp  = ttmp**(1.d0/(gamma-1.d0))
		 ptmp    = ttmp*rhotmp

		 utmp    = uamb + upert
		 vtmp    =        vpert

		 etmp    = ptmp/(gamma-1.d0) + 
     .                     .5d0*rhotmp*(utmp**2 +vtmp**2)

		 val(i,j,1) = rhotmp
		 val(i,j,2) = rhotmp*utmp
		 val(i,j,3) = rhotmp*vtmp
		 val(i,j,4) = etmp
 120         continue
          endif
 130   continue

      else if (iprob .eq. 2) then
c
c     2d riemann problem initial cond. (see Glaz et al SISSC 14, 1993)
c
	xmid = .5d0 * (xupper+xlower)
	ymid = .5d0 * (yupper+ylower)
	do 50 i = 1, mitot
	xcen = cornx + (i-.5)*hx
	do 50 j = 1, mjtot
	ycen = corny + (j-.5)*hy
	    if (xcen .le. xmid) then
	       if (ycen .le. ymid) then
c                  # 3rd quadrant data
c                  #configuration 3
c		   rho =  0.137992831541219d0
c		   pr  = 0.029032258064516d0
c		   u   =  1.206045378311055d0
c		   v   =  1.206045378311055d0
c                  #configuration B
c                  rho = 1.d0
c                  u   = -.75d0
c                  v   =  .5d0
c                  pr  = 1.d0
c                  #configuration A
                   rho = 1.d0
                   u   =  .75d0
                   v   =  .5d0
                   pr  = 1.d0
	       else
c                  # 2nd quadrant data
c                  #configuration 3
c		   rho = 0.532258064516129d0
c		   pr  = 0.3d0
c		   u   =  1.206045378311055d0
c		   v   =  0.d0
c                  #configuration B
c                  rho = 2.d0
c                  u   =  .75d0
c                  v   =  .5d0
c                  pr  = 1.d0
c                  #configuration A
                   rho = 2.d0
                   u   = -.75d0
                   v   =  .5d0
                   pr  = 1.d0
	       endif
	    else if (ycen .le. ymid) then
c                  # 4th quadrant data
c                  #configuration 3
c		   rho = 0.532258064516129d0
c		   pr  = 0.3d0
c		   u   = 0.d0
c		   v   =  1.206045378311055d0
c                  #configuration B
c                  rho = 3.d0
c                  u   = -.75d0
c                  v   = -.5d0
c                  pr  = 1.d0
c                  #configuration A
                   rho = 3.d0
                   u   =  .75d0
                   v   = -.5d0
                   pr  = 1.d0
	    else
c                  # 1st quadrant data
c                  #configuration 3
c		   rho = 1.5d0
c		   pr  = 1.5d0
c		   u   = 0.d0
c		   v   = 0.d0
c                  #configuration B
c                  rho = 1.d0
c                  u   =  .75d0
c                  v   = -.5d0
c                  pr  = 1.d0
c                  #configuration A
                   rho = 1.d0
                   u   = -.75d0
                   v   = -.5d0
                   pr  = 1.d0
	    endif

	    val(i,j,1) = rho
	    val(i,j,2) = rho * u
	    val(i,j,3) = rho * v
	    val(i,j,4) = pr/gamma1 + .5d0*rho*(u*u + v*v)

 50      continue

      else if (iprob .eq. 5) then

c     # mach 10 30 degree ramp reflection. shock starts at x= xshock at wall
c     # parameters set in setprob.
c     # 2nd ramp problem: mach 1.025 at 18 degrees.

         do 60 i = 1, mitot
	   xcen = cornx + (i-.5)*hx

	   do 60 j = 1, mjtot
	     ycen = corny + (j-.5)*hy

             call cellave(xcen-hx/2.d0,ycen-hy/2.d0,hx,hy,wl)
             rho = (1.d0-wl)*rhoamb + wl*rhoshk
             u =   (1.d0-wl)*uamb   + wl*ushk
             v =   (1.d0-wl)*vamb   + wl*vshk
             p =   (1.d0-wl)*pamb   + wl*pshk
         
	     val(i,j,1) = rho
	     val(i,j,2) = rho * u
	     val(i,j,3) = rho * v
	     val(i,j,4) = p/gamma1 + .5d0*rho*(u*u+v*v)

 60      continue


      endif

      return
      end
