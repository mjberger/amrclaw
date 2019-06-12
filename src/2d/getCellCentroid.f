c
c -------------------------------------------------------------------------
c
      subroutine getCellCentroid(lstgrd,i,j,xc,yc,xlow,ylow,dx,dy,k)

      use amr_module
      implicit double precision (a-h,o-z)

      if (k .eq. lstgrd .or. k .eq. -1) then
         xc = xlow + (i-0.5d0)*dx
         yc = ylow + (j-0.5d0)*dy
      else
         xc = xcirr(k)
         yc = ycirr(k)
      endif

      return
      end
