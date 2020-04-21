c
c ------------------------------------------------------------------
c
      subroutine physbd(xleft,xright,ybot,ytop,
     1                  level,nrow,ncol,
     2                  nvar,val,time,hx,hy)
 
c
c :::::::::: PHYSBD ::::::::::::::::::::::::::::::::::::::::::::::;
c
c     # dummy physbd term routine for problems with no physical boundary
c     # conditions, (e.g. if periodic in both directions)

c :::::::::: PHYSBD ::::::::::::::::::::::::::::::::::::::::::::::;
c
c     Take a grid patch with mesh widths hx,hy, of dimensions nrow by
c     ncol,  and set the values of any piece of
c     of the patch which extends outside the physical domain
c     using the boundary conditions.
c
c     The corners of the grid patch are at
c        (xleft,ybot)  --  lower left corner
c        (xright,ytop) --  upper right corner
c
c     The physical domain itself is a rectangle bounded by
c        (xlower,ylower)  -- lower left corner
c        (xupper,yupper)  -- upper right corner
c
c     the picture is the following:
c
c               _____________________ (xupper,yupper)
c              |                     |
c          _________ (xright,ytop)   |
c          |   |    |                |
c          |   |    |                |
c          |   |    |                |
c          |___|____|                |
c (xleft,ybot) |                     |
c              |                     |
c              |_____________________|
c   (xlower,ylower)
c
c
c     Any cells that lie outside the physical domain are ghost cells whose
c     values should be set in this routine.  This is tested for by comparing
c     xleft with xlower to see if values need to be set at the left, as in
c     the figure above, and similarly at the other boundaries.
c
c     Patches are guaranteed to have at least 1 row of cells filled
c     with interior values so it is possible to  extrapolate.
c     Fix trimbd if you want more than 1 row pre-set.
c
c     Make sure the order the boundaries are specified is correct
c     so that diagonal corner cells are also properly taken care of.
c
c     Periodic boundaries are set before calling this routine, so you
c     can safely extrapolate there.  Don't overwrite them!
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;


      implicit double precision (a-h,o-z)

      include  "cuser.i"
      include  "cdisc.i"
      include  "cdom.i"

      dimension val(nrow,ncol,nvar)

      hxmarg = hx*.01
      hymarg = hy*.01
      pi     = atan(1.d0)*4.d0

      if (xperdom .and. yperdom) return

      if (iprob .eq. 1) then

c     left boundary - inflow
      if (xleft .lt. xlower-hxmarg) then
            nxl = (xlower+hxmarg-xleft)/hx
	    do 10 i = 1, nxl
 	    do 10 j = 1, ncol
 	      val(i,j,1) = rhoshk
 	      val(i,j,2) = rhoshk * ushk
 	      val(i,j,3) = 0.d0
 	      val(i,j,4) = rhoshk*(eshk+.5d0*ushk**2)
 10         continue
      endif

c     right boundary. - extrapolate
       if (xright .gt. xupper+hxmarg) then
         nxr = (xright - xupper + hxmarg)/hx
         nxr = nxr - 1
         ibeg = max0(nrow-nxr, 1)
 	     do 11 i    = ibeg, nrow
 	     do 11 j    = 1, ncol
 	     do 11 ivar = 1, nvar
 	       val(i,j,ivar) = val(2*ibeg-1-i,j,ivar) 
 11           continue
      endif

c     bottom boundary - reflecting
      if (ybot .lt. ylower-hymarg) then
        nyb = (ylower+hymarg-ybot)/hy
	   do 12 j=1,nyb     
	   do 12 i=1,nrow
	      val(i,j,1) =  val(i,2*nyb+1-j,1)
	      val(i,j,2) =  val(i,2*nyb+1-j,2)
	      val(i,j,3) = -val(i,2*nyb+1-j,3)
	      val(i,j,4) =  val(i,2*nyb+1-j,4)
  12       continue
c
      endif
 
c     top boundary - reflecting
c     top boundary - extrpolate
      if (ytop .gt. yupper-hymarg) then
	  nyt = (ytop - yupper + hymarg)/hy
	  jbeg = max0(ncol-nyt+1, 1)
	     do 13 j=jbeg,ncol
	     do 13 i=1,nrow
		val(i,j,1) =  val(i,2*jbeg-1-j,1)
		val(i,j,2) =  val(i,2*jbeg-1-j,2)
 		val(i,j,3) = -val(i,2*jbeg-1-j,3)
		val(i,j,4) =  val(i,2*jbeg-1-j,4)
  13         continue
      endif

c end of bc's for problem 1 
      endif
c
c
      if (iprob .eq. 2) then

c     bottom boundary  - extrapolate
      if (ybot .lt. ylower-hymarg) then
        nyb = (ylower+hymarg-ybot)/hy
           do 20 j = 1,nyb
           do 20 i = 1,nrow
              val(i,j,1) = val(i,nyb+1,1)
              val(i,j,2) = val(i,nyb+1,2)
              val(i,j,3) = val(i,nyb+1,3)
              val(i,j,4) = val(i,nyb+1,4)
20         continue
      endif

c     top boundary 
      if (ytop .gt. yupper-hymarg) then
        nyt = (ytop - yupper + hymarg)/hy
        jbeg = max0(ncol-nyt+1, 1)
           do 21 j= jbeg,ncol
           do 21 i    = 1, nrow
              val(i,j,1) = val(i,jbeg-1,1)
              val(i,j,2) = val(i,jbeg-1,2)
              val(i,j,3) = val(i,jbeg-1,3)
              val(i,j,4) = val(i,jbeg-1,4)
21         continue
       endif

c     right boundary 
      if (xright .gt. xupper-hxmarg) then
        nxr = (xright - xupper + hxmarg)/hx
        nxr = nxr - 1
        ibeg = max0(nrow-nxr, 1)
           do 22 i = ibeg, nrow
           do 22 j = 1, ncol
             val(i,j,1) = val(ibeg-1,j,1)
             val(i,j,2) = val(ibeg-1,j,2)
             val(i,j,3) = val(ibeg-1,j,3)
             val(i,j,4) = val(ibeg-1,j,4)
22         continue
      endif

c     left boundary - 
      if (xleft .lt. xlower-hxmarg) then
        nxl = (xlower+hxmarg-xleft)/hx
           do 23 i = 1,nxl
           do 23 j = 1,ncol
              val(i,j,1) = val(nxl+1,j,1)
              val(i,j,2) = val(nxl+1,j,2)
              val(i,j,3) = val(nxl+1,j,3)
              val(i,j,4) = val(nxl+1,j,4)
23         continue
        endif

c end of bc's for problem 2 
      endif


      if (iprob .eq. 5) then

c     right boundary. 
      if (xright .gt. xupper+hxmarg) then
        nxr = (xright - xupper + hxmarg)/hx
        nxr = nxr - 1
c       # extrapolate
        ibeg = max0(nrow-nxr, 1)
	 do 50 i    = ibeg, nrow
	 do 50 j    = 1, ncol
	 do 50 ivar = 1, nvar
	   val(i,j,ivar) = val(2*ibeg-1-i,j,ivar) 
50       continue
        endif


c     bottom boundary. 
      if (ybot .lt. ylower-hymarg) then
        nyb = (ylower+hymarg-ybot)/hy
	 do 51 i = 1, nrow
	 do 51 j = 1, nyb
	   xcen   = xleft + dfloat(i-.5d0)*hx
	   ycen   = ybot  + dfloat(j-.5d0)*hy
	   if (xcen .le. disp) then
	     rho = rhoshk
	     u = ushk
	     v = vshk
	     e = eshk
	   else
	     rho = val(i,2*nyb+1-j,1)
	     u =   val(i,2*nyb+1-j,2) / rho
	     v = - val(i,2*nyb+1-j,3) / rho
	     e = (val(i,2*nyb+1-j,4) -.5d0*rho*(u*u+v*v))/rho
	   endif
	   val(i,j,1) = rho
	   val(i,j,2) = rho*u
	   val(i,j,3) = rho*v
	   val(i,j,4) = rho*(e+.5d0*(u*u+v*v))
 51       continue
        endif
 
c     top boundary.  
      if (ytop .gt. yupper-hymarg) then
	  nyt = (ytop - yupper + hymarg)/hy
	  jbeg = max0(ncol-nyt+1, 1)
c	  topinit = (yupper-ylower)/1.732 + disp 
c	  xshock = 10.d0*time/sin(pi/3.d0) + topinit

c	  topinit = (yupper-ylower)*sin(pi/10.d0) + disp 
c	  xshock = 1.025d0*time*cos(pi/10.d0) + topinit
c	  x0 = xshock
c         y0 = yupper
          
	  x0 = disp + 1.025d0*time/cos(pi/10.d0)
	  y0 = ylower

	  do 52 j = jbeg, ncol
	    ycen   =  ybot + dfloat(j-.5d0)*hy

	     do 52 i = 1, nrow
	      xcen   = xleft + dfloat(i-.5d0)*hx

	      call cellave(xcen-hx/2.d0,ycen-hy/2.d0,hx,hy,wl)
	      rho = (1.d0-wl)*rhoamb + wl*rhoshk
	      u   = (1.d0-wl)*uamb + wl*ushk
	      v   = (1.d0-wl)*vamb + wl*vshk
	      p   = (1.d0-wl)*pamb + wl*pshk

	      val(i,j,1) =  rho
	      val(i,j,2) =  rho * u
	      val(i,j,3) =  rho * v
	      val(i,j,4) =  p/gamma1 + .5d0*rho*(u*u+v*v)
 52       continue
      endif

c     left boundary
      if (xleft .lt. xlower-hxmarg) then
        nxl = (xlower+hxmarg-xleft)/hx
	    do 53 i = 1, nxl
	    do 53 j = 1, ncol
	      val(i,j,1) = rhoshk
	      val(i,j,2) = rhoshk * ushk
	      val(i,j,3) = rhoshk * vshk
	      val(i,j,4) = rhoshk*(eshk+.5d0*(ushk**2+vshk**2))
 53        continue
        endif
c
c end of bc's for problem 5 
      endif


      if (iprob .eq. 6) then
c
c no problem 6 right now, but leave for anything needing
c all reflecting boundary conditions.

c     bottom boundary - reflecting
      if (ybot .lt. ylower-hymarg) then
        nyb = (ylower+hymarg-ybot)/hy
	   do 62 j=1,nyb     
	   do 62 i=1,nrow
	      val(i,j,1) =  val(i,2*nyb+1-j,1)
	      val(i,j,2) =  val(i,2*nyb+1-j,2)
	      val(i,j,3) = -val(i,2*nyb+1-j,3)
	      val(i,j,4) =  val(i,2*nyb+1-j,4)
  62       continue
c
      endif
 
c     top boundary - reflecting
      if (ytop .gt. yupper-hymarg) then
	  nyt = (ytop - yupper + hymarg)/hy
	  jbeg = max0(ncol-nyt+1, 1)
	     do 63 j=jbeg,ncol
	     do 63 i=1,nrow
		val(i,j,1) =  val(i,2*jbeg-1-j,1)
		val(i,j,2) =  val(i,2*jbeg-1-j,2)
		val(i,j,3) = -val(i,2*jbeg-1-j,3)
		val(i,j,4) =  val(i,2*jbeg-1-j,4)
  63         continue
      endif

c     left boundary - reflecting
      if (xleft .lt. xlower-hxmarg) then
        nxl = (xlower+hxmarg-xleft)/hx
           do 64 i = 1,nxl
           do 64 j = 1,ncol
              val(i,j,1) =  val(2*nxl+1-i,j,1)
              val(i,j,2) = -val(2*nxl+1-i,j,2)
              val(i,j,3) =  val(2*nxl+1-i,j,3)
              val(i,j,4) =  val(2*nxl+1-i,j,4)
64         continue
        endif

c     right boundary - reflecting
      if (xright .gt. xupper+hxmarg) then
        nxr  = (xright - xupper + hxmarg)/hx
        nxr  = nxr - 1
        ibeg = max0(nrow-nxr, 1)
           do 65 i = ibeg,nrow
           do 65 j = 1,ncol
              val(i,j,1) =  val(2*ibeg-1-i,j,1)
              val(i,j,2) = -val(2*ibeg-1-i,j,2)
              val(i,j,3) =  val(2*ibeg-1-i,j,3)
              val(i,j,4) =  val(2*ibeg-1-i,j,4)
65         continue
       endif

c
c end of bc's for problem 6 
c
      endif


      return
      end
