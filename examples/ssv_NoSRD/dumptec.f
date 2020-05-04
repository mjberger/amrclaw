c
c -----------------------------------------------------
c
      subroutine dumptec (lst,lend,nvar,naux,nplot,time)
c
      use amr_module
      implicit double precision (a-h,o-z)
      common /userdt/ cflcart,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                ismp,gradThreshold,pwconst
      common /order2/ ssw, quad, nolimiter
      logical         flag,pwconst,quad,nolimiter
      character*23  filename 
      character*11   filebndry
c
c dumptec = make tecplot file for finest level soln values over entire domain
c (but for now assume only 1 level  -- output finest level)
c
c nplot == plotnum (0 for initial conditions, 1 (or more) for final conditions
c
      filename = 'disjointErrPlanesxx.dat'
      filename(18:18) = '0'
      filename(19:19) = char(ichar('0')+nplot)
      open(13,file=filename,status='unknown',form='formatted')
      write(13,100) 

      filename = 'disjointCutPlanesxx.dat'
      filename(18:18) = '0'
      filename(19:19) = char(ichar('0')+nplot)
      open(14,file=filename,status='unknown',form='formatted')
      write(14,100) 
 100  format('TITLE = "Extracted cutting Planes through mesh"' )

      ! also dump boundary plots at same time
      filebndry = 'bndry0x.dat'
      filebndry(7:7) = char(ichar('0')+nplot)
      ibunit = 24
      open(ibunit,file=filebndry,status='unknown',form='formatted')
      write(ibunit,*)"# writing bndry file at time ",time
      

      nplot = nplot+1

 10   level = lfine

c     initialize for error computation
      volDenErrorL1   = 0.d0
      volExactDenL1   = 0.d0
      volDenErrorMax  = 0.d0
      exactVol        = 0.d0
      bndryDenErrorL1 = 0.d0
      bndryExactDenL1 = 0.d0
      bndryCentExactDenL1 = 0.d0
      bndryReconErrL1 = 0.d0
      exactBndry      = 0.d0
      aftDenErrorL1   = 0.d0

c     fill ghost cells if need them for output or to compute gradients
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

             !!if (ssw .ne. 0 .and. .not. pwconst)
               call bound(time,nvar,nghost,alloc(locnew),mitot,mjtot,
     2                    mptr,alloc(locaux),naux)

              ! remember got 3 times size of irr to include other arrays
              locirr = node(permstore,mptr)
              locncount = locirr + mitot*mjtot
              locnumHoods = locncount + mitot*mjtot

              call outtec(alloc(locnew),nvar,mptr,
     1                    alloc(locirr),mitot,mjtot,
     2                    lstgrd,hx,hy,xlow,ylow,time,
     3                    alloc(locncount),alloc(locnumHoods),ibunit,
     4              volDenErrorL1,volExactDenL1,volDenErrorMax,
     5              exactVol,bndryDenErrorL1,bndryExactDenL1,
     6              bndryCentExactDenL1,bndryReconErrL1,
     7              exactBndry,aftDenErrorL1)
              write(*,*)" after aftDenError  = ",aftDenError
c
              mptr = node(levelptr,mptr)
          go to 20
50        continue
 
c
      close(14)
      close(24)
      write(*,*)"done writing ",filename," at time",time

c output errors
      write(outunit,600) volDenErrorL1,volExactDenL1,
     .  volDenErrorL1/volExactDenL1,exactVol,volDenErrorMax
 600  format("L1 density volume error ", e15.7,/,
     .       "L1 density exact soln   ", e15.7,/,
     .       "L1 Relative density error",e15.7,/,
     .       "Computed volume          ",e15.7,/,
     .       "Max density error        ",e15.7,//)
      write(outunit,601) bndryDenErrorL1, bndryCentExactDenL1,
     .                   bndryDenErrorL1/bndryCentExactDenL1
 601  format("L1 Bndry Density Error   ",e15.7,/,
     .       "L1 Bndry exact soln      ",e15.7,/,
     .       "L1 relative density error",e15.7,//)
      write(outunit,602) bndryReconErrL1, bndryExactDenL1,
     .                   bndryReconErrL1/bndryExactDenL1
 602  format("L1 Recon2Bndry   Error   ",e15.7,/,
     .       "L1 Bndry exact soln      ",e15.7,/,
     .       "L1 reconstructed relative  error",e15.7,//)

      write(outunit,603) exactBndry
 603  format("Length of Bndry segments ", e15.7,//)

      write(outunit,604) aftDenErrorL1
 604  format("Aftosmis relative error  ",e15.7)
c
 99   return
      end
