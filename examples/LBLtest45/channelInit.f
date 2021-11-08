c
c-----------------------------------------------------------------------
c
      subroutine channelInit(xcen,ycen,rho,u,v,p)
      implicit double precision(a-h,o-z)
c
c     # For the point (x,y), compute the true solution at this point 
c

       if (xcen .lt. 2.d0) then
          rho = 1.d0
       else
          rho = 0.0d0
       endif

       ! leave nvar 4 in case can't be changed
       u = 0.d0
       v = 0.d0
       p = 1.d0

       return
       end
