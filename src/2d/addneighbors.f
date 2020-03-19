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
            nco = 1 ! reg cell has 4 nbors so will have enough 
            !nco = 2 ! do quadratic fit for both cut and one adjacent
            cutoffDistReg = 4 ! for reg cell this means use all nbors
            !cutoffDistReg = 3 ! for reg cell this means use all nbors
            ! cutoffDistReg = 2 ! means use only 2 away in l1 norm
            ! cutoffDistReg = 1 means use edge nbors only
         else
            nco = 2 ! try larger stencil for cut cells?
c           nco = 1 ! will have enough nbors
            ! count how many edge nbors for use below
            nb = 0
            do in = 1, 4
              if (irr(ix0+nxb(in),iy0+nyb(in)) .ne. -1) nb = nb + 1
            end do
            if (nb .ge. 4) then
               cutoffDist = 2
            else
               cutoffDist = 3
            endif
            cutoffDist = 6
         endif
         do 20 ioff = -nco,nco
         do 20 joff = -nco,nco
           if (ioff .eq. 0 .and. joff .eq. 0) go to 20
c           ::: prerequisite for edge sharing if one of the offsets = 0
c           ## if next line commented out, then diagonal nbors are allwoed
           ixn = ix0 + ioff
           iyn = iy0 + joff
           if (outofbounds(ixn,iyn)) go to 20

c          this next line says to use is edge nbors only if not commented out
           dist1 =  l1_len(ixn,iyn)

           ! if reg cell use only edge nbors (dist 1)
           if (k0.eq.lstgrd .and. dist1 .gt. cutoffDistReg) go to 20

           ! if cut cell include vertex nbors (dist 2)
           ! but never allow more than 2 (test in case nco = 2)
           !!if (k0.ne.lstgrd .and. dist1 .gt. 2) go to 20  
           ! new test - use edge nbors only for cut cell if have enough
           if (k0.ne.lstgrd .and. dist1 .gt. cutoffDist) go to 20  

           kn  = irr(ixn,iyn)
           if (kn .eq. -1) go to 20
           if (kn .eq. lstgrd  .or. k0 .eq. lstgrd) go to 25
c           ::: both cells cut. must check that really share
c           ::: a common edge
           go to 25  ! skip next check for larger stencil
           do 23 kside = 1,6
              if (poly(kside+2,1,k0) .eq. -11) go to 20
              x1 = poly(kside,1,k0)
              y1 = poly(kside,2,k0)
              x2 = poly(kside+1,1,k0)
              y2 = poly(kside+1,2,k0)
              if (x1 .ne. x2) then
                 ihoriz = 1
                 ivert  = 0
              else
                 ihoriz = 0
                 ivert  = 1
              endif
c     compute indices of adjacent cell
              ixr = ix0 + ivert*isig(y1-y2)
              iyr = iy0 + ihoriz*isig(x2-x1)
              if ((ixr .eq. ixn) .and. (iyr .eq. iyn)) go to 25
 23        continue
c          ::: should never fall through
c     
c          ::: make sure new cell not already on list
 25        do 30 ncheck = 1, newend
              if (ixn .eq. nlist(ncheck,1) .and.
     .            iyn .eq. nlist(ncheck,2)) go to 20
 30        continue
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
           if (verbose .and. ar(kn)/ar(lstgrd).lt.gradThreshold) then
            go to 20
           endif
           newend = newend + 1
           nlist(newend,1) = ixn 
           nlist(newend,2) = iyn 
 20     continue

 10   continue

       return
       end
