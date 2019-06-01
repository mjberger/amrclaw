c
c ---------------------------------------------------------------
c
       subroutine addneighbors(irr,nlist,nst,nend,newend,mitot,mjtot,
     .            lstgrd,quad,xlow,ylow,hx,hy)

       use amr_module
       implicit double precision (a-h, o-z)

      common   /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                  ismp,gradThreshold
       dimension nlist(25,2), irr(mitot,mjtot)
       logical outside, quad
       logical verbose/.false./
       outside(i,j) = (i .lt. 1 .or. i .gt. mitot .or.
     &                 j .lt. 1 .or. j .gt. mjtot)
       l1_len(ixn,iyn) = iabs(ixn-nlist(1,1))+iabs(iyn-nlist(1,2))


       newend = nend

       do 10 n = nst, nend
c        ::: add neighbors of this point to nlist
         ix0 = nlist(n,1)
         iy0 = nlist(n,2)
         k0  = irr(ix0,iy0)

         do 20 ioff = -1, 1    
         do 20 joff = -1, 1    
           if (ioff .eq. 0 .and. joff .eq. 0) go to 20
c           ::: prerequisite for edge sharing if one of the offsets = 0
c           ## if next line commented out, then diagonal nbors are allwoed
c           if (ioff.ne.0 .and. joff.ne.0) go to 20 
           ixn = ix0 + ioff
           iyn = iy0 + joff
           if (outside(ixn,iyn)) go to 20
c           ::: dont use points 2 away in each coordinate direction
           if (l1_len(ixn,iyn) .ge. 3) go to 20
           kn  = irr(ixn,iyn)
           if (kn .eq. -1) go to 20
           if (kn .eq. lstgrd  .or. k0 .eq. lstgrd) go to 25
c           ::: both cells cut. must check that really share
c           ::: a common edge
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
            write(*,*)"not using ",ixn,iyn," volfrac ",ar(kn)/ar(lstgrd)
            write(*,*)"    in neighbor stencil for cell",ix0,iy0
            go to 20
           endif
           newend = newend + 1
           nlist(newend,1) = ixn 
           nlist(newend,2) = iyn 
 20     continue

 10   continue

       return
       end
