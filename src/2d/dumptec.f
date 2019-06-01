c
c -----------------------------------------------------
c
      subroutine dumptec (lst,lend,nvar,naux,nplot,time)
c
      use amr_module
      implicit double precision (a-h,o-z)
      common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold
      logical         flag
      character*23  filename 
c
c dumptec = make tecplot file for finest level soln values over entire domain
c (but for now assume only 1 level  -- output finest level)
c
c nplot == plotnum (0 for initial conditions, 1 (or more) for final conditions
c
      filename = 'disjointCutPlanesxx.dat'
      filename(18:18) = '0'
      filename(19:19) = char(ichar('0')+nplot)
      nplot = nplot+1

      open(14,file=filename,status='unknown',form='formatted')
      write(*,*)" writing graphics file ",filename," at time ",time
      write(14,100) 
 100  format('TITLE = "Extracted cutting Planes through mesh"' )

 10   level = lfine
      mptr = lstart(level)
 20       if (mptr .eq. 0) go to 50
              nx = node(ndihi,mptr)-node(ndilo,mptr)+1
              ny = node(ndjhi,mptr)-node(ndjlo,mptr)+1
              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)
              mitot = nx + 2*nghost
              mjtot = ny + 2*nghost
              locnew = node(store1,mptr)
              locaux = node(storeaux,mptr)
              locirr = node(permstore,mptr)
              lstgrd = node(lstptr,mptr)
              xl   = rnode(cornxlo, mptr)
              yb   = rnode(cornylo, mptr)
              xr   = rnode(cornxhi, mptr)
              yt   = rnode(cornyhi, mptr)
              hx   = hxposs(level)
              hy   = hyposs(level)

             call bound(time,nvar,nghost,alloc(locnew),mitot,mjtot,mptr,
     1                  alloc(locaux),naux)

              ! remember got 3 times size of irr to include other arrays
              locirr = node(permstore,mptr)
              locncount = locirr + mitot*mjtot
              locnumHoods = locncount + mitot*mjtot

              call outtec(alloc(locnew),nvar,mptr,
     1                    alloc(locirr),mitot,mjtot,
     2                    lstgrd,hx,hy,xlow,ylow,time,
     3                    alloc(locncount),alloc(locnumHoods))
c
              mptr = node(levelptr,mptr)
          go to 20
50        continue
c
      close(14)
c
 99   return
      end
