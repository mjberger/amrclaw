c
c ---------------------------------------------------------------
c
       subroutine addneighbors(irr,nlist,nst,nend,newend,mitot,mjtot,
     .            lstgrd,quad,xlow,ylow,hx,hy)

       use amr_module
       implicit double precision (a-h, o-z)

       include "cuserdt.i"
       dimension nlist(25,2), irr(mitot,mjtot)
       logical outofbounds, outside, quad
       logical verbose/.false./
       integer nxb(4) 
       data nxb/-1,1,0,0/
       integer nyb(4) 
       data nyb/0,0,-1,1/
       integer cutoffDist
       outofbounds(i,j) = (i .lt. 1 .or. i .gt. mitot .or.
     &                 j .lt. 1 .or. j .gt. mjtot)
       l1_len(ixn,iyn) = iabs(ixn-nlist(1,1))+iabs(iyn-nlist(1,2))


       newend = nend

       do 10 n = nst, nend
c        ::: add neighbors of this point to nlist
         ix0 = nlist(n,1)
         iy0 = nlist(n,2)
         k0  = irr(ix0,iy0)

         if (k0 .eq. lstgrd) then
            nco = 1 ! reg cell, use 3 by 3 block around cell
         else
            nco = 2 ! use 5 by 5 for cut cell since half are missing
         endif

         do 20 ioff = -nco,nco
         do 20 joff = -nco,nco
           if (ioff .eq. 0 .and. joff .eq. 0) go to 20
c           ::: prerequisite for edge sharing if one of the offsets = 0
c           ## if next line commented out, then diagonal nbors are allwoed
           ixn = ix0 + ioff
           iyn = iy0 + joff
           if (outofbounds(ixn,iyn)) go to 20

           kn  = irr(ixn,iyn)
           if (kn .eq. -1) go to 20
c     
c          ::: make sure new cell not already on list
c25        do 30 ncheck = 1, newend
c             if (ixn .eq. nlist(ncheck,1) .and.
c    .            iyn .eq. nlist(ncheck,2)) go to 20
c30        continue
           if (kn .ne. lstgrd) then
 1            xn = xcirr(kn)
              yn = ycirr(kn)
           else
              xn = xlow + (ixn-.5d0)*hx
              yn = ylow + (iyn-.5d0)*hy
           endif
c          ::: dont use new cell if ghost cell, want 1-sided diffs like cart3d
c          if (xn .gt. xprob .or. xn .lt. 0.) go to 20
c          ### dont use new cell if area way too small (same test for no gradients in this cell)
c          in irreg3hbox
           newend = newend + 1
           nlist(newend,1) = ixn 
           nlist(newend,2) = iyn 
 20     continue

 10   continue

       return
       end
