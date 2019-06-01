c
c-----------------------------------------------------------------------
c
      subroutine p20fn(xcen,ycen,q,time)

      implicit double precision(a-h,o-z)
      dimension q(1)
c
c     # For the point (x,y), compute the true solution at this point 
c
       xshift = 1.5
       yshift = 1.5
       pi = 3.14159265357989d0

       yuse = ycen - yshift
       xuse = xcen - xshift
       radius = dsqrt(yuse**2+xuse**2)
       sr     = 2.d0*pi/5.d0   ! goes around one full circle in 5 time units

       theta  = atan2(yuse,xuse)

       thetause = theta - sr*time
       xinit = cos(thetause)
       yinit = sin(thetause)
       
       thetainit = atan2(yinit,xinit)

       w = .5d0*( erf((pi/6.d0-thetainit)/sqrt(4.d0/100.d0)) + 
     &            erf((pi/6.d0+thetainit)/sqrt(4.d0/100.d0)) )

c      theta = atan2(ycen-yshift,xcen-xshift)
c      w = .5d0*( erf((pi/6.d0-theta)/sqrt(4.d0/100.d0)) + 
c    &            erf((pi/6.d0+theta)/sqrt(4.d0/100.d0)) )

       q(1) = w
c      q(1) = 1.0
c   stick in prob 16 into prob 20
c        q(1) = ycen - .1*xcen + .5

       return
       end
