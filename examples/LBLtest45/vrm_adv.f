c
c -------------------------------------------------------------
c  UPWINDING FOR LINEAR ADVECTION
C
      subroutine vrm(qr,ql,rx,ixmin,ixmax,iflip,msize, mptr)

      implicit double precision (a-h,o-z)
      parameter ( m1size = 5000 )

      dimension qr(1,msize),ql(1,msize)
      dimension rx(1,msize)
      dimension dir(2), dn(2)

      if (ixmax .gt. m1size) then
         write(6,*)" need to increase m1size in vrm to ", ixmax
         stop
      endif


      call getdir(dir)

      dn = 0.d0
      if(iflip .eq. 1) then
        dn(1) = 1.d0
      elseif(iflip .eq. 0) then
        dn(2) = 1.d0
      endif

      dot = dir(1)*dn(1) + dir(2)*dn(2)

      rx = 0.d0
      do 10 k = ixmin, ixmax
         rx(1,k) = 0.5d0*dot * (ql(1,k)+qr(1,k))
     .            +0.5d0*abs(dot)*( ql(1,k)-qr(1,k) )
  10     continue

      return
      end


      subroutine getdir(dir)
      implicit double precision (a-h,o-z)
      dimension dir(2)
      dir(1) = 1.d0
      dir(2) = 1.d0
      return
      end
