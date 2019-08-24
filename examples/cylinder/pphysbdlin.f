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
 
 
      rho = 8.0d0   ! 5.7143d0*1.40d0
      u = 8.25d0
      v = .0d0
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
             val(1,i,j) = val(1,i,jbeg-1)
             val(2,i,j) = val(2,i,jbeg-1)
             val(3,i,j) = val(3,i,jbeg-1)
             val(4,i,j) = val(4,i,jbeg-1)

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
               val(3,i,j) = val(3,i,nyb+1)
               val(4,i,j) = val(4,i,nyb+1)
200        continue
c
      endif
c
 99   return
      end
