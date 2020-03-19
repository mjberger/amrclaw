c
c ---------------------------------------------------------------------
c
       subroutine qslopes(qp,qx,qy,mitot,mjtot,irr,lstgrd,lwidth,
     &                    hx,hy,xlow,ylow,mptr,nvar) 
      use amr_module
      implicit double precision(a-h,o-z)

      dimension qp(nvar,mitot,mjtot),qx(nvar,mitot,mjtot),
     &          qy(nvar,mitot,mjtot),irr(mitot,mjtot)
      logical nolimiter
      common /order2/ ssw, quad, nolimiter
      include "cuserdt.i"


      dimension a(25,5),at(5,25),c(5,5),b(25,nvar),d(25,nvar)
      dimension rhs(25,nvar), w(5,nvar)
      dimension nlist(25,2)
      logical   prflag, quad, enufNbor, outside, almost_outside
      data      prflag/.false./
      logical ALL_NBORS_EXIST


c
c sandras variables
c
      parameter (max_num_nbor=22)
      logical  LP_success, output_debug
      double precision neighb_x(max_num_nbor)
      double precision neighb_y(max_num_nbor)
      double precision neighb_u(max_num_nbor) 
      double precision Phi_x, Phi_y, LP_Dx,LP_Dy, LS_Dx, LS_Dy
      integer row_A
      logical do_pos_constr ! are we really a cutcell?
      integer my_iters, tot_iters, max_iters   ! maximum no of iterations needed in Simplex
      double precision avg_iters    ! avg number of iterations needed in Simplex
      integer count_cells
      logical added_neighbor 
c
c     outside(x,y) = ((i.lt.1) .or. (j.lt.1.) .or.
c    .                (i.gt.mitot) .or. (j.gt.mjtot))
      outside(i,j) = ((i.lt.1) .or. (j.lt.1.) .or.
     .                (i.gt.mitot) .or. (j.gt.mjtot))

      ALL_NBORS_EXIST(i,j) = (i .gt. 1 .and. i .lt.mitot .and.
     .                        j .gt. 1 .and. j .lt. mjtot)

c
c   ##########
c   #  compute slopes for cut cells using least squares approach
c   ##########
c   
c   ssw = 0 for ctu (no slopes), 1 for muscl (second order).
c   qp contains primitive variables, regular cell slopes already set 
c   all slopes initialized to zero in regular slopes subr.
c
c
c     first save everything into unlim (even though limited already
c     on flow cells. will add unlimited cut gradients

      quad = .false.

      count_cells = 0
      tot_iters = 0
      max_iters = 0

      do 110 ix0 = 2, mitot-1
      do 110 iy0 = 2, mjtot-1
c     do 110 ix0 = 1, mitot
c     do 110 iy0 = 1, mjtot
         k = irr(ix0,iy0)
         if (k .eq. -1) go to 110
         if (k .eq. lstgrd .and. ALL_NBORS_EXIST(ix0,iy0)) then
           if (irr(ix0+1,iy0) .eq. lstgrd .and. 
     .         irr(ix0,iy0+1) .eq. lstgrd .and. 
     .         irr(ix0,iy0-1) .eq. lstgrd .and. 
     .         irr(ix0-1,iy0) .eq. lstgrd) go to 110
         endif
c     
         if (ar(k)/ar(lstgrd) .lt. gradThreshold) then
            qx(:,ix0,iy0) = 0.d0
            qy(:,ix0,iy0) = 0.d0
            go to 110           ! leave 0 gradient:  more stable for teeny cells w/o slopes
         endif
c     
c      # this cell needs a slope
         if (k .ne. lstgrd) then
            x0 = xcirr(k)
            y0 = ycirr(k)
            do 37 kside = 1, 6
               if (poly(kside+2,1,k).eq.-11)then
                  sidex2 = poly(kside,1,k)
                  sidey2 = poly(kside,2,k)
                  sidex1 = poly(kside+1,1,k)
                  sidey1 = poly(kside+1,2,k)
                  go to 39
               endif
 37         continue
 39         rlenb = dsqrt((sidey1-sidey2)**2 + (sidex1-sidex2)**2)
            alf  = (sidey1-sidey2)/rlenb
            beta = (sidex2-sidex1)/rlenb
            bxpt = .5d0*(sidex1+sidex2)
            bypt = .5d0*(sidey1+sidey2)
         else
            x0 = xlow + (ix0-.5d0)*hx
            y0 = ylow + (iy0-.5d0)*hy
         endif
         
c     if (outside(x0,y0)) go to 110   ! to match cart3d, no gradients in ghost cells
         nlist(1,1) = ix0
         nlist(1,2) = iy0
         nst        = 1  
         nend       = 1  
         call addneighbors(irr,nlist,nst,nend,newend,mitot,mjtot,
     .        lstgrd,quad,xlow,ylow,hx,hy)


c        ::: not enough neighbors for quadratic fit
c        ::: add neighbors of neighbors    
c
c         if ((irow .lt. 2)) then
c          write(*,*) "no slope for cut cell ",ix0,iy0
c          go to 110    ! not enough neighbors for gradient recon
c         endif

 16      irow = 0
         if (newend .le. 2) then !could happen since starting at r/c=1
            qx(:,ix0,iy0) = 0.d0
            qy(:,ix0,iy0) = 0.d0
            go to 110    ! leave 0 gradient
         endif

         nToUse = nTerms
         if (nToUse .eq. 5 .and. newend .lt. 7) then
            nToUse = 2
         endif


         do 22 n = 2, newend
            irow = irow + 1
            ixn = nlist(n,1)
            iyn = nlist(n,2)
            kn =  irr(ixn,iyn)
            if (kn .ne. lstgrd) then
               xn = xcirr(kn)
               yn = ycirr(kn)
            else
               xn = xlow + (ixn-.5d0)*hx
               yn = ylow + (iyn-.5d0)*hy
            endif
c     # code treating q as pt.wise values
            delx = xn - x0
            dely = yn - y0
            if (k .eq. lstgrd) then
               a(irow,1) = delx
               a(irow,2) = dely
               if (nToUse .eq. 5) then
c               # code doing quadratic but pointwise not cell averaged right
                 a(irow,3) = (xn-x0)**2
                 a(irow,4) = (xn-x0)*(yn-y0)
                 a(irow,5) = (yn-y0)**2
               endif
            else
               deln = alf*delx + beta*dely
               delt = -beta*delx + alf*dely
               a(irow,1) = deln
               a(irow,2) = delt
               if (nToUse .eq. 5) then
                   a(irow,3) = deln*deln
                   a(irow,4) = deln*delt
                   a(irow,5) = delt*delt
               endif
            endif
c
 21         do m = 1, nvar
               b(irow,m) = qp(m,ixn,iyn) - qp(m,ix0,iy0)
            end do
 22      continue
 
c  if cut cell add ghost cell state  (but not used for limiting?)
           if (.false. .and. ghost_ccg .and. k .ne. lstgrd) then
              distx = x0 - bxpt
              disty = y0 - bypt
              d_dot_n = distx*alf + disty*beta
              dNormx = d_dot_n * alf
              dNormy = d_dot_n * beta
              dx = 2.d0*dNormx
              dy = 2.d0*dNormy
              xghost = x0 - dx
              yghost = y0 - dy

              irow = irow+1
              a(irow,1) = xghost - x0
              a(irow,2) = yghost - y0
              ! check that signs cancel, since my normal is inward pointing
              uvel = qp(2,ix0,iy0)
              vvel = qp(3,ix0,iy0)
              dot = uvel*alf + vvel*beta
              ughost = uvel - 2.d0*dot*alf
              vghost = vvel - 2.d0*dot*beta

              b(irow,1) = 0.d0 ! pw const extrap of density
              b(irow,2) = ughost - uvel
              b(irow,3) = vghost - vvel
              b(irow,4) = 0.d0  ! and pressure
           endif
c

          do 30 it = 1, irow
          do 30 jt = 1, nToUse
             at(jt,it) = a(it,jt)
 30       continue

          do 50 it = 1, nToUse
             do 50 jt = 1, nToUse
                c(it,jt) = 0.d0
                do m = 1, nvar
                   d(it,m)  = 0.d0
                end do
                do 45 kt = 1, irow
                   c(it,jt) = c(it,jt) + at(it,kt)*a(kt,jt)
                   do m = 1, nvar
                      d(it,m) = d(it,m) + at(it,kt)*b(kt,m)
                   end do
 45             continue
 50       continue
c
c       do linear fit

c     # solve C*w = d for least squares slopes. use cholesky
c     # put factors back in a
           if (nToUse .eq. 2) then
             a(1,1) = dsqrt(c(1,1))
             a(1,2) = c(1,2)/a(1,1)
             a(2,2) = dsqrt(c(2,2)-a(1,2)**2)
c     
             do 61 m = 1, nvar
c     
c     # at*a = c. solve at*b = d, aw = b.  reuse b.
                b(1,m) = d(1,m) / a(1,1)
                b(2,m) = (d(2,m) - a(1,2)*b(1,m)) / a(2,2)
                w2 =   b(2,m)/a(2,2)
                w1 = (b(1,m)-a(1,2)*w2)/a(1,1)
                if (k .eq. lstgrd) then
                   qy(m,ix0,iy0) =  w2
                   qx(m,ix0,iy0) = w1
                else
                   sn = w1
                   st = w2
                   ! rotate back for x and y gradients from sn and st
                   qy(m,ix0,iy0) = beta*sn + alf*st
                   qx(m,ix0,iy0) = alf*sn - beta*st  
                endif
 61          continue

           else ! nterms 5, do quadratic
c         # now solve C*w = d for least squares slopes. use cholesky
c         # put factors back in a
             a(1,1) = dsqrt(c(1,1))
             a(1,2) = c(1,2)/a(1,1)
             a(1,3) = c(1,3)/a(1,1)
             a(1,4) = c(1,4)/a(1,1)
             a(1,5) = c(1,5)/a(1,1)

             a(2,2) = dsqrt(c(2,2)-a(1,2)**2)
             a(2,3) = (c(2,3)-a(1,2)*a(1,3))/a(2,2)
             a(2,4) = (c(2,4)-a(1,2)*a(1,4))/a(2,2)
             a(2,5) = (c(2,5)-a(1,2)*a(1,5))/a(2,2)

             a(3,3) = dsqrt(c(3,3)-a(1,3)**2 - a(2,3)**2)
             a(3,4) = (c(3,4)-a(1,3)*a(1,4)-a(2,3)*a(2,4))/a(3,3)
             a(3,5) = (c(3,5)-a(1,3)*a(1,5)-a(2,3)*a(2,5))/a(3,3)

             a(4,4) = dsqrt(c(4,4)-a(1,4)**2-a(2,4)**2-a(3,4)**2)
             a(4,5) = (c(4,5)-a(1,4)*a(1,5)-a(2,4)*a(2,5)-a(3,4)*a(3,5))
     &            /a(4,4)
             a(5,5) = dsqrt(c(5,5)-a(1,5)**2-
     &                      a(2,5)**2-a(3,5)**2-a(4,5)**2)
c     
             do 60 m = 1, nvar
c     
c     # at*a = c. solve at*b = d, aw = b.  reuse b.
                b(1,m) = d(1,m) / a(1,1)
                b(2,m) = (d(2,m) - a(1,2)*b(1,m)) / a(2,2)
                b(3,m) = (d(3,m) - a(1,3)*b(1,m)-a(2,3)*b(2,m))/a(3,3)
                b(4,m) = (d(4,m) - a(1,4)*b(1,m)-a(2,4)*b(2,m)
     &                    -a(3,4)*b(3,m))/a(4,4)
                b(5,m) = (d(5,m) - a(1,5)*b(1,m)-a(2,5)*b(2,m)
     &                    -a(3,5)*b(3,m)-a(4,5)*b(4,m))/a(5,5)


                w(5,m) = b(5,m) / a(5,5)
                w(4,m) = (b(4,m)-a(4,5)*w(5,m))/a(4,4)
                w(3,m) = (b(3,m)-a(3,5)*w(5,m)-a(3,4)*w(4,m))/a(3,3)
                w(2,m) = (b(2,m)-a(2,5)*w(5,m)-a(2,4)*w(4,m)
     &               -a(2,3)*w(3,m))/a(2,2)
                w(1,m) = (b(1,m)-a(1,5)*w(5,m)-a(1,4)*w(4,m)
     &               -a(1,3)*w(3,m)-a(1,2)*w(2,m))/a(1,1)

                if (k .eq. lstgrd) then
                   qy(m,ix0,iy0)  =  w(2,m)
                   qx(m,ix0,iy0)  =  w(1,m)
                else 
                   sn = w(1,m)
                   st = w(2,m)
                   ! rotate back for x and y gradients from sn and st
                   qy(m,ix0,iy0) = beta*sn + alf*st
                   qx(m,ix0,iy0) = alf*sn - beta*st  
                endif
c
 60       continue
              
           endif

c          write(*,401) ix0,iy0,(qx(m,ix0,iy0),qy(m,ix0,iy0),m=1,4)
 401      format("cell ",2i5,  " orig lsq grad ",/,4(2e15.7,/))

c       should we limit?
        if (nolimiter)  go to 110

c
          center_x = x0
          center_y = y0
c         if (added_neighbor) then
c            num_neighb = newend-2  ! don't use diagonal neighbor for limiting
c         else
c            num_neighb = newend-1 ! not counting cell itself
c         end if

c         use diagonal neighbor
          num_neighb = newend-1 ! not counting cell itself
          output_debug = .false.


          do m = 1, 4   ! loop over the 4 variables
            row_A     = 2*num_neighb + 4
            do_pos_constr = .false.

c           if (irr(ix0,iy0) .ne. lstgrd) then
c              if (m .eq. 1 .OR. m .eq. 4) then
c                 do_pos_constr = .true.
c                 row_A = row_A + 1 ! additional pos constraint for rho, p
c              endif
c           endif

c           for NT cut cell limiting, redo so x0,y0 = 0,0

            do inbor = 1, num_neighb
              ixn = nlist(inbor+1,1)
              iyn = nlist(inbor+1,2)
              kn =  irr(ixn,iyn)
              ! in this version subtract centroid, then reset it to 0,0
              if (kn .ne. lstgrd) then
                 xn = xcirr(kn) 
                 yn = ycirr(kn)
              else
                 xn = xlow + (ixn-.5d0)*hx 
                 yn = ylow + (iyn-.5d0)*hy
              endif

              if (k .ne. lstgrd) then ! try limiting in normal and tang
              ! take rotated vals out of at, subtract rotated center
                !dxn = at(1,inbor) - x0
                !dyn = at(2,inbor) - y0
                xn = xn - x0
                yn = yn - y0
                dxn = alf*xn + beta*yn
                dyn = -beta*xn + alf*yn
                xn = dxn
                yn = dyn
              endif

              neighb_x(inbor) = xn ! neighboring centroid, not diffs
              neighb_y(inbor) = yn ! neighboring centroid, not diffs


              neighb_u(inbor) = qp(m,ixn,iyn)  ! neighboring u, need to loop over m
            end do  ! do inbor = 1, num_neighb
              center_u = qp(m,ix0,iy0)   
              LS_Dx = qx(m,ix0,iy0)
              LS_Dy = qy(m,ix0,iy0)
              ! initialize without rotation first
              bxptRot = 0.d0
              byptRot = 0.d0
              if (k .ne. lstgrd) then
                 LS_Dx = w(1,m)
                 LS_Dy = w(2,m)
                 center_x = 0.d0
                 center_y = 0.d0
                 bxptRot =   alf*(bxpt-x0) + beta*(bypt-y0)
                 byptRot = -beta*(bxpt-x0) +  alf*(bypt-y0)
              endif

            LP_success = .false.
            if (k .eq. 201) then
               output_debug = .true.
            endif

             call vector_lim_allineq(num_neighb,neighb_x,neighb_y,
     &               neighb_u,center_x,center_y,center_u,output_debug,
     &               row_A,Phi_x,Phi_y, LP_Dx,LP_Dy,LP_success,
     &           LS_Dx,LS_Dy,bxptRot,byptRot,do_pos_constr,
     &            my_iters,lpChoice)

c           Increase the count of number of iterations
            tot_iters = tot_iters + my_iters
            if (my_iters .GT. max_iters) then
               max_iters = my_iters
            endif
            count_cells = count_cells + 1

c             write(*,402) LS_Dx,LS_Dy,LP_Dx,LP_Dy
 402         format("sandras lsq ",2e15.7,10x," limited ",2e15.7)
         
             if (LP_success) then
             if (k .eq. lstgrd) then
                  qx(m,ix0,iy0) = LP_Dx
                  qy(m,ix0,iy0) = LP_Dy
             else
                 qx(m,ix0,iy0) =  alf*LP_DX - beta*LP_DY
                 qy(m,ix0,iy0) = beta*LP_DX + alf*LP_DY
             endif
c                  write(*,400),m, Phi_x, Phi_y
 400              format("limiter values for m = ",i3," are ",2e15.7)
             else
                  write(*,*)"LP failure for grid cell var ",
     &                      mptr,ix0,iy0,m
                  qx(m,ix0,iy0) = 0.d0
                  qy(m,ix0,iy0) = 0.d0

c                 stop
             endif

          end do   ! do m=1,4
c  
 110  continue

      if (count_cells > 0) then
         avg_iters = dble(tot_iters)/dble(count_cells)
      else
         avg_iters = 0
      endif
      write(21,*) "Avg_no_of_iters", avg_iters
      write(23,*) "Max_no_of_iters", max_iters
c
      if (prflag) then
         write(21,*)' qx '
         do 180 i = 2, mitot-1
         do 180 j = 2, mjtot-1
            if (irr(i,j) .ne. -1) then
               write(21,190)i,j,(qx(m,i,j),m=1,nvar)
 190           format('  i,j  ',2i4,4e14.6)
            endif
 180     continue
         write(21,*)' qy '
         do 181 i = 2, mitot-1
         do 181 j = 2, mjtot-1
            if (irr(i,j) .ne. -1) then
               write(21,190)i,j,(qy(m,i,j),m=1,nvar)
            endif
 181     continue

      endif
c
 99   return
      end
