c
c
c---------------------------------------------------------------------
c 
      integer function isig(x)
      implicit double precision (a-h,o-z)
c
      if (x.gt.0.d0) then
          isig = 1
        else
          isig = -1
        endif
      return
      end
