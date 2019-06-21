c
c ------------------------------------------------------------------
c
      subroutine pphysbdlin(xleft,xright,ybot,ytop,level,nrow,ncol,nvar,
     1                  val,time,hx,hy,qx,qy,irr,lstgrd)
 
c     This routine takes an (enlarged) grid (or grid patch)
c     with mesh widths hx,hy, and sets the values of any piece of
c     of the patch which extends outside the physical domain using the
c     values given by the boundary conditions. 
c
c 
c
      use amr_module
      implicit double precision (a-h,o-z)

      common /userdt/ cflcart,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     &                ismp,gradThreshold
      dimension val(nvar,nrow,ncol)
      dimension qx(nvar,nrow,ncol), qy(nvar,nrow,ncol)
      dimension irr(nrow,ncol)
 
 

      pi = 3.14159265358979d0
      sloc = .39d0
      alpha = 30.d0*pi/180.d0
      rho = 8.0d0   
      u = 8.25d0*cos(alpha)
      v = -8.25d0*sin(alpha)
      pr = 116.5d0
c     rho = 1.4d0
c     pi = 3.14159265358979d0
c     u = cos(pi/6.d0)
c     v = sin(pi/6.d0)
c     pr = 1.d0
c
      hxmarg = hx*.01
      hymarg = hy*.01

c     left boundary (inflow)
 
      if (xleft .lt. -hxmarg) then
 
        nxl = (hxmarg-xleft)/hx
 
           do 400 i = 1,nxl
           do 400 j = 1,ncol
               val(1,i,j) = rho
               val(2,i,j) = u
               val(3,i,j) = v
               val(4,i,j) = pr
400        continue
c
      endif
 
c     top boundary (extrap from interior)
 
      if (ytop .gt. yprob+hymarg) then
 
        nyt = (ytop - yprob + hymarg)/hy
        jbeg = max0(ncol-nyt+1, 1)
 
           do 100 j= jbeg,ncol
           do 100 i    = 1, nrow
             ycen = ybot  + (dfloat(j)-.5d0)* hy
             xcen = xleft  + (dfloat(i)-.5d0)* hx
             sloc_top = sloc + (ycen+20.d0*time)/sqrt(3.d0)
             if (xcen < sloc_top) then
                 rhot = 8.d0
                 ut = 8.25d0*cos(pi/6.d0)
                 vt = -8.25d0*sin(pi/6.d0)
	         pt =  116.5d0
             else
	         rhot = 1.4d0
                 ut = 0.d0
                 vt = 0.d0
	         pt = 1.d0
	     endif
             val(1,i,j) = rhot
             val(2,i,j) = ut
             val(3,i,j) = vt
             val(4,i,j) = pt

100        continue
c
      endif
 
c     right boundary. 
      if (xright .gt. xprob+hxmarg) then
 
        nxr = (xright - xprob + hxmarg)/hx
        nxr = nxr - 1
 
        ibeg = max0(nrow-nxr, 1)

c start extrap at bottom of grid, not including ghost cells
           if (ybot .lt. -hymarg) then
              jbeg = (hymarg-ybot)/hy + 1
           else
              jbeg = 1
           endif

           do 300 j = jbeg, ncol
           do 300 i = ibeg, nrow
            ycen = ybot  + (dfloat(j)-.5d0)* hy
            xcen = xleft + (dfloat(i)-.5d0)* hx
               val(1,i,j) = val(1,ibeg-1,j) 
               val(2,i,j) = val(2,ibeg-1,j)
               val(3,i,j) = val(3,ibeg-1,j)
               val(4,i,j) = val(4,ibeg-1,j)
300        continue
c
      endif
 
 
c     bottom boundary. 
 
      if (ybot .lt. -hymarg) then
        nyb = (hymarg-ybot)/hy
           do 200 j = 1,nyb
           do 200 i = 1,nrow   ! the first lwidth set in left bc above
               val(1,i,j) = val(1,i,nyb+1)
               val(2,i,j) = val(2,i,nyb+1)
               val(3,i,j) = -val(3,i,nyb+1) ! reflecting
               val(4,i,j) = val(4,i,nyb+1)
200        continue
c
      endif
c
 99   return
      end
