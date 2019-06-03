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
      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,ismp
      dimension val(nrow,ncol,nvar)
      dimension qx(nrow,ncol,nvar), qy(nrow,ncol,nvar)
      dimension irr(nrow,ncol)
      common  /cirr/ poly(10,2,5933),ar(-1:5933),
     *               points(24,2,5933),wt(24,-1:5933),
     *               xcirr(5933),ycirr(5933),ixg(5933),iyg(5933),
     *               nxtirr(5933)
 
 
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
               val(i,j,1) = rho
               val(i,j,2) = u
               val(i,j,3) = v
               val(i,j,4) = pr
400        continue
c
      endif
 
c     top boundary (extrap from interior)
 
      if (ytop .gt. yprob+hymarg) then
 
        nyt = (ytop - yprob + hymarg)/hy
        jbeg = max0(ncol-nyt+1, 1)
 
           do 100 j= jbeg,ncol
           do 100 i    = 1, nrow
             val(i,j,1) = val(i,jbeg-1,1)
             val(i,j,2) = val(i,jbeg-1,2)
             val(i,j,3) = val(i,jbeg-1,3)
             val(i,j,4) = val(i,jbeg-1,4)

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
               val(i,j,1) = val(ibeg-1,j,1) 
               val(i,j,2) = val(ibeg-1,j,2)
               val(i,j,3) = val(ibeg-1,j,3)
               val(i,j,4) = val(ibeg-1,j,4)
300        continue
c
      endif
 
 
c     bottom boundary. 
 
      if (ybot .lt. -hymarg) then
        nyb = (hymarg-ybot)/hy
           do 200 j = 1,nyb
           do 200 i = 1,nrow   ! the first lwidth set in left bc above
               val(i,j,1) = val(i,nyb+1,1)
               val(i,j,2) = val(i,nyb+1,2)
               val(i,j,3) = val(i,nyb+1,3)
               val(i,j,4) = val(i,nyb+1,4)
200        continue
c
      endif
c
 99   return
      end
