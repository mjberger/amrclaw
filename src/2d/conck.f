c
c -----------------------------------------------------------
c
!> Conservation check for specified level.
!! This is mostly a debugging tool and assumes grids don't overlap
      subroutine conck(level, nvar, naux, time, rest)
c
      use amr_module
      implicit double precision (a-h,o-z)

      logical  rest

c ## indexing into mcapa assumes cell volume is in mcapa location
      iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)*mitot+i-1)
      iaddaux(i,j) = locaux + mcapa - 1 + naux*(i-1) +
     .                                    naux*mitot*(j-1)
c
c
c ******************************************************************
c conck - conservation check  for specified level
c         mostly a debugging tool
c         this assumes grids don't overlap
c
c ******************************************************************
c
c
c  grid loop for given level
c
      hx      = hxposs(level)
      hy      = hyposs(level)
      dt      = possk(level)
      totmass = 0.d0

      mptr = lstart(level)
 20   if (mptr .eq. 0) go to 85
         loc    = node(store1,mptr)
         locirr = node(permstore,mptr)
         locaux = node(storeaux,mptr)
         nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
         ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
         mitot  = nx + 2*nghost
         mjtot  = ny + 2*nghost
         lstgrd = node(lstptr,mptr)
         ar(lstgrd) = hx * hy ! regular grid cell area
         ar(-1) = 0.d0        ! so solid cells dont count
c
         if (mcapa .eq. 0) then
           do 50 j  = nghost+1, mjtot-nghost
           do 50 i  = nghost+1, mitot-nghost
              k =  kget(mitot,mjtot,i,j,alloc(locirr)) 
              area = ar(k)
              totmass = totmass + area*alloc(iadd(1,i,j)) 
 50           continue
          else
c          # with capa array:
           do 60 j  = nghost+1, mjtot-nghost
           do 60 i  = nghost+1, mitot-nghost
              totmass = totmass + alloc(iadd(1,i,j))*alloc(iaddaux(i,j)) 
 60           continue
          endif
c
       mptr = node(levelptr,mptr)
       go to 20
c
 85    continue
       if (time.eq. t0 .and. (level.eq.1) .and. .not. rest) then
           tmass0 = totmass
           write(6,*) 'Total mass at initial time: ',tmass0
           endif
       write(outunit,777) time, totmass, totmass-tmass0
 777   format('time t = ',e12.5,',  total mass = ',e22.15, '  diff = ',
     &         e11.4)
c
 99   return
      end
c
c ----------------------------------------------
c
      integer function kget(mitot,mjtot,i,j,irr)
      dimension irr(mitot,mjtot)
      kget = irr(i,j)
      return 
      end
