c
c ---------------------------------------------------------------------
c    limit one cell at a time. When chagne data structure to store
c    nhood for each cell, then put loop around it and wont need nlist
c
       subroutine limitCellLP(qp,qx,qy,mitot,mjtot,irr,lstgrd,
     &                        lwidth,hx,hy,xlow,ylow,mptr,nvar,
     &                        nlist,num_neighb,num_iters,ix0,iy0,x0,y0) 

      use amr_module
      implicit double precision(a-h,o-z)

      dimension qp(nvar,mitot,mjtot),qx(nvar,mitot,mjtot),
     &          qy(nvar,mitot,mjtot),irr(mitot,mjtot)
      dimension nlist(25,2)
      logical nolimiter, quad
      common /order2/ ssw, quad, nolimiter
      include "cuserdt.i"
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
c
c     count_cells = 0
c     tot_iters = 0
c     max_iters = 0
c
          center_x = x0
          center_y = y0
c         if (added_neighbor) then
c            num_neighb = newend-2  ! don't use diagonal neighbor for limiting
c         else
c            num_neighb = newend-1 ! not counting cell itself
c         end if

          num_neighb = newend-1 ! not counting cell itself
          output_debug = .false.


          do m = 1, nvar    ! loop over the 4 variables
            row_A     = 2*num_neighb + 4
            do_pos_constr = .false.

c           if (irr(ix0,iy0) .ne. lstgrd) then
c              if (m .eq. 1 .OR. m .eq. 4) then
c                 do_pos_constr = .true.
c                 row_A = row_A + 1 ! additional pos constraint for rho, p
c              endif
c           endif

            do inbor = 1, num_neighb
              ixn = nlist(inbor+1,1)
              iyn = nlist(inbor+1,2)
              kn =  irr(ixn,iyn)
              if (kn .ne. lstgrd) then
                 xn = xcirr(kn)
                 yn = ycirr(kn)
              else
                 xn = xlow + (ixn-.5d0)*hx
                 yn = ylow + (iyn-.5d0)*hy
              endif

              neighb_x(inbor) = xn ! neighboring centroid, not diffs
              neighb_y(inbor) = yn ! neighboring centroid, not diffs

              neighb_u(inbor) = qp(m,ixn,iyn)  ! neighboring u, need to loop over m
            end do  ! do inbor = 1, num_neighb
              center_u = qp(m,ix0,iy0)   
              LS_Dx = qx(m,ix0,iy0)
              LS_Dy = qy(m,ix0,iy0)

            LP_success = .false.

             call vector_lim_allineq(num_neighb,neighb_x,neighb_y,
     &               neighb_u,center_x,center_y,center_u,output_debug,
     &               row_A,Phi_x,Phi_y, LP_Dx,LP_Dy,LP_success,
     &           LS_Dx,LS_Dy,bxpt,bypt,do_pos_constr,my_iters,lpChoice)

c           Increase the count of number of iterations
c           tot_iters = tot_iters + my_iters
c           if (my_iters .GT. max_iters) then
c              max_iters = my_iters
c           endif
c           count_cells = count_cells + 1

c             write(*,402) LS_Dx,LS_Dy,LP_Dx,LP_Dy
 402         format("sandras lsq ",2e15.7,10x," limited ",2e15.7)
         
             if (LP_success) then
                  qx(m,ix0,iy0) = LP_Dx
                  qy(m,ix0,iy0) = LP_Dy
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
c     if (count_cells > 0) then
c        avg_iters = dble(tot_iters)/dble(count_cells)
c     else
c        avg_iters = 0
c     endif
c     write(21,*) "Avg_no_of_iters", avg_iters
c     write(23,*) "Max_no_of_iters", max_iters
c
      return
      end
