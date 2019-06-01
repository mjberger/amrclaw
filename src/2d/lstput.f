c
c -------------------------------------------------------------
c
      subroutine lstput(lstgrd)
c
      use amr_module
      implicit double precision (a-h,o-z)
c
c lstgrd = start of linked list for the grid. is never 0, since
c the "regular cell" info. (at least area) uses at least
c one spot on the info. list.
c
c find the end of the grid's list
c
       lptr = lstgrd
 10    if (nxtirr(lptr) .eq. 0) go to 20
          lptr = iabs(nxtirr(lptr))
          go to 10
c
c now lptr points to end of the grid's list. attach to rest of list
c
 20    nxtirr(lptr) = lhead
       lhead = lstgrd
c
       return
       end
