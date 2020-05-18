c
!> print solution and aux. variables to output. 
c =======================================================================
      subroutine outval(val,nvar,mitot,mjtot,mptr,outgrd,naux,aux,
     &                  irr,lstgrd)
c =======================================================================
c
      use amr_module
      implicit double precision (a-h,o-z)

      dimension  val(nvar,mitot,mjtot)
      dimension  primval(nvar,mitot,mjtot)
      dimension  aux(naux,mitot,mjtot)
      integer    irr(mitot,mjtot)
      logical    outgrd


c ::::::::::::::::::::::OUTVAL :::::::::::::::::::::::::::::::
c print solution and aux. variables to output. 
c if cell outside domain, don't print soln. value - nothing
c currently in ghost cells.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      if (.not. outgrd) go to 99
      level = node(nestlevel,mptr)
      hx    = hxposs(level)
      hy    = hyposs(level)
      cornx = rnode(cornxlo,mptr) -  nghost*hx
      corny = rnode(cornylo,mptr) -  nghost*hy
      time = rnode(timemult,mptr)

      ! need ghost cells to avoid zero divides
      call bound(time,nvar,nghost,val,mitot,mjtot,mptr,aux,naux)
      call vctoprm(val,primval,mitot,mjtot,nvar)
c
      do 20 i=nghost+1,mitot-nghost
      write(outunit,*)
      do 25 j=nghost+1,mjtot-nghost

          if (irr(i,j) .eq. -1) cycle
          x  = cornx + hx*(dble(i)-.5d0)
          y  = corny + hy*(dble(j)-.5d0)
          if (irr(i,j) .eq. lstgrd) then
          write(outunit,107) x,y,i,j,(primval(ivar,i,j),ivar=1,nvar)
 107      format(1x,2hx=,f7.3,3h y=,f7.3,4h, i=,i4,4h, j=,i4,' a=',
     *           4(e10.4,1x))
          else
             write(outunit,108) x,y,i,j,(primval(ivar,i,j),ivar=1,nvar)
 108         format('*',2hx=,f7.3,3h y=,f7.3,4h, i=,i4,4h, j=,i4,' a=',
     *              4(e10.4,1x))
          endif
          if (naux.gt.0) write(outunit,109) (aux(iaux,i,j),iaux=1,naux)
 109      format(1x,'aux = ',7(e10.3,1x))

 25   continue
 20   continue

 99   return
      end
