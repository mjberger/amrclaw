c
c-----------------------------------------------------------------------
c
      subroutine channelPtInit(xcen,ycen,rho,u,v,p)
      implicit double precision(a-h,o-z)
c
c     # For the point (x,y), compute the true solution at this point 
c

       rho = ycen**2 
       !rho = ycen**2 - 0.1d0*xcen**2 + 5.d0 + 2.d0*xcen*ycen
       !rho = ycen - 0.1d0*xcen + 0.5d0
       !rho = ycen + 0.5d0
       u = 0.d0
       v = 0.d0
       p = 1.d0

       return
       end
c
c-----------------------------------------------------------------------
c
      subroutine channelAvgInit(poly,xcen,ycen,rho,u,v,p)

      implicit double precision(a-h,o-z)
      dimension poly(10,2)
      dimension q(4,3), qcent(20,4), area(20), celltot(4)
c
c     # For the centroid (x,y), compute the true cell avg.  To
c      compute cell average (to 3rd order). decompose
c      polygon into triangles, with centroid always
c      a vertex.  use approx. that cell avg. is
c      average of the function values on the 3 sides

      do 10 kside = 2, 20
         if (poly(kside,1) .eq. -11) then
            kend = kside - 2
            go to 99
         endif

         xst = poly(kside-1,1)
         yst = poly(kside-1,2)

         xend = poly(kside,1)
         yend = poly(kside,2)

         xavg = (xst+xend)/2.d0
         yavg = (yst+yend)/2.d0
         call channelPtInit(xavg,yavg,rho1,u1,v1,p1)

         xavg = (xend+xcen)/2.d0
         yavg = (yend+ycen)/2.d0
         call channelPtInit(xavg,yavg,rho2,u2,v2,p2)

         xavg = (xst+xcen)/2.d0
         yavg = (yst+ycen)/2.d0
         call channelPtInit(xavg,yavg,rho3,u3,v3,p3)
c     
c     compute area of triangle
         area(kside-1) = .5d0*((yst+ycen)*(xcen-xst) +
     &        (ycen+yend)*(xend-xcen)+(yend+yst)*(xst-xend))
            qcent(kside-1,1) = (rho1+rho2+rho3)/3.d0
            qcent(kside-1,2) = (u1+u2+u3)/3.d0
            qcent(kside-1,3) = (v1+v2+v3)/3.d0
            qcent(kside-1,4) = (p1+p2+p3)/3.d0

 11      continue
 10   continue

 99   continue

         totar   = 0.d0
         celltot = 0.d0
         do 12 k = 1, kend
            celltot(1) = area(k)*qcent(k,1) + celltot(1)
            celltot(2) = area(k)*qcent(k,2) + celltot(2)
            celltot(3) = area(k)*qcent(k,3) + celltot(3)
            celltot(4) = area(k)*qcent(k,4) + celltot(4)
            totar   = totar + area(k)
 12      continue

      rho = celltot(1) / totar
      u   = celltot(2) / totar
      v   = celltot(3) / totar
      p   = celltot(4) / totar

      return
      end

