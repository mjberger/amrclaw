c
c-----------------------------------------------------------------------
c
      subroutine channelInit(xcen,ycen,rho,u,v,p)
      implicit double precision(a-h,o-z)
c
c     # For the point (x,y), compute the true solution at this point 
c

       !rho = ycen - 0.1d0*xcen + 0.5d0
       rho = ycen + 0.5d0
       u = 0.d0
       v = 0.d0
       p = 1.d0

       return
       end
