c
c -------------------------------------------------------------
c
      integer function lstget(dummy)
c
      use amr_module
c
c first free spot in linked list of irregular info. pointed to
c by lhead.
c
      if (lhead .ne. 0) go to 10
         write(6,100)
 100     format(' out of irregular list space ')
         stop
c
 10   lstget = lhead
      lhead = iabs(nxtirr(lhead))
c
      return
      end
