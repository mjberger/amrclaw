c
c ------------------------------------------------------------------------
c
      subroutine checkPhys(q,irr,mitot,mjtot,mptr,istage,lstgrd,str,
     &                     ist,iend,jst,jend)

      implicit real*8 (a-h,o-z)
      dimension q(4,mitot,mjtot)
      integer irr(mitot,mjtot)
      character*14 str

      gamma1 = .4d0
      do j = jst, jend
      do i = ist, iend
         k = irr(i,j)
         if (k .eq. -1) cycle
         rho = q(1,i,j)
         u = q(2,i,j)/rho
         v = q(3,i,j)/rho
         pr = gamma1*(q(4,i,j)-.5d0*rho*(u*u+v*v))
         if (rho < 0 .or. pr < 0) then
           if (k .ne. lstgrd) then
            write(*,901)rho,pr,i,j,mptr,istage,str
 901        format("non-physical den/pr",2e15.7," at i,j*  grid stage ",
     .           4i5,2x,a14)
           else
            write(*,902)rho,pr,i,j,mptr,istage,str
 902        format("non-physical den/pr",2e15.7," at i,j  grid stage ",
     .           4i5,2x,a14)
           endif
         endif
      end do
      end do

      return
      end
