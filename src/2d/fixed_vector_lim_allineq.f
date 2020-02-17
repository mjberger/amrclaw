c     ##########################################################
c     # Program for vector limiting for cutcells               #
c     #                                                        #
c     ##########################################################



      subroutine vector_lim_allineq(num_neighb,neighb_x,neighb_y,
     &    neighb_u,center_x,center_y,center_u,output_debug,row_A,
     & Phi_x,Phi_y,LP_Dx,LP_Dy,LP_success,LS_Dx,LS_Dy,bxpt,bypt,
     & do_pos_constr,my_iters,lpChoice)




c :::::::::::::::::::::::::::: VECTOR_LIM ::::::::::::::::::::::
c
c    This routine does:
c     * set up the 2d LP 
c             - build the matrix A and rhs b
c             - initialize x0 and c appropriately
c             - call the subroutine 'Simplex_allineq' -- that's Simplex in allineq form
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c     uses the restrictive constraints (don't overshoot and correct direction)
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c  INPUT:
c     * num_neighb: number of neighboring cells used for limiting
c     * neighb_x(num_neighb), neighb_y(num_neighb): x- and y-coordinates of neighbors
c     * neighb_u(num_neighb): cell value of neighbors
c     * center_x,center_y,center_u: data for central cell
c     * output_debug: how much output do you want?
c     * row_A,column_A: dimension of constraint matrix A in LP
c     * LS_Dx,LS_Dy: least square gradient
c     * bxpt,bypt: x- & y-coord of the boundary edge midpoint
c     * do_pos_constr: if true, then add pos constraint for boundary edge midpoint
c  
c  OUTPUT:
c     * Phi_x,Phi_y: scalar factors that limit LS_Dx to LP_Dx (same for y)
c     * LP_Dx,LP_Dy: results of LP: limited slopes
c     * LP_success: was the Simplex algorithm successful?
c     * my_iters: number of iterations needed in the Simplex-algorithm
c
c :::::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::



c     ::::::::::::::::::::::::::::
c     Declare variables
c     ::::::::::::::::::::::::::::

      implicit none  ! force yourself to properly declare all var

c     Input
      integer num_neighb
      double precision neighb_x(num_neighb), neighb_y(num_neighb) 
      double precision neighb_u(num_neighb)
      double precision center_x, center_y, center_u
      logical output_debug
      integer row_A ! number of constraints in LP, ie number of rows in matrix
      ! number of rows: 2
      double precision LS_Dx, LS_Dy, bxpt, bypt
      logical do_pos_constr

c     Output
      double precision Phi_x, Phi_y 
      double precision LP_Dx, LP_Dy 
      logical LP_success
      integer my_iters

c     Local variables: check num_neighb not too big
      integer max_num_neighb
      parameter (max_num_neighb = 22)
      logical output_debug_here
      logical temp_output

c     Local variables: for least square problem
c      double precision S(num_neighb,2) ! matrix for LS problem
c      double precision S_rhs(num_neighb),rhs_LS(2) ! rhs's for LS problem
c      double precision LS_Dx, LS_Dy    ! solution of LS problem
c      double precision normal_matrix(2,2) ! S'*S
c      character   TRANSA, TRANSB
c      integer ind_perm(2)  ! permutation array coming back from solving matrix
c      integer INFO ! was subroutine of LAPACK successful?

c     Local variables: Formulating the 2d LP
      double precision A(row_A,2)   ! matrix of constraints for LP
      double precision b(row_A)  ! rhs for LP
      double precision x0(2)  ! number of unknowns
      integer W(2) ! working set
      double precision c(2) ! vector c in objective function
      double precision x(2)     ! solution vector of LP
      double precision     rhsmax, rhsmin
      double precision     pred

c     Local variables: Loop variables
      integer i,j,nn, choice_limiting,lpChoice

c     choice = 1 => standard constraints
c     choice = 2 => relaxed constraints      
c      choice_limiting = 1;
c     choice_limiting = 2;
      choice_limiting = lpChoice


c     ::::::::::::::::::::::::::::
c     Declare finished
c     ::::::::::::::::::::::::::::


 902     format (A,'\t',2(E12.2))
 901     format (A,'\t',F6.3)
 903     format (4e15.7,3x)


c     Check whether num_neighb < max_num_neighb
      if (num_neighb .GT. max_num_neighb) then
         write(*,*) "\n##########\n"
         write(*,901) "num_neighb too big: num_neighb is\t", num_neighb
         write(*,*) "\n##########\n"
         return
      endif


      if (output_debug) then
         write(*,*) "############\n Given data \n ############"
         write(*,*) "center_x       center_y             center_u"
         write(*,903) center_x,center_y,center_u
         write(*,*)"neighb_x          neighb_y      neighb_u   ",
     &             "  predicted"
         do nn = 1, num_neighb
            pred = center_u + (neighb_x(nn)-center_x)*LS_DX +
     &                        (neighb_y(nn)-center_y)*LS_DY 
            write(*,903) neighb_x(nn),neighb_y(nn),neighb_u(nn),pred
         end do
         write(*,*)"bxpt       bypt"
         write(*,903) bxpt, bypt
         write(*,*)"LS_DX      LS_DY"
         write(*,903) LS_DX, LS_DY
      endif


c     #####################################
c     Set up and solve least square problem
c     #####################################


cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     LS solution is passed in
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc


      if (choice_limiting .eq. 1) then
      

c     #######################
c     Formulate 2d LP
c     #######################

         if (output_debug) then
            write(*,*) "############\n Formulate 2d LP \n ############"
         endif

c            ::::::::::::::::::::::
c            initialize A and b
c            ::::::::::::::::::::::

c            Initialize with 0
             do 60 j = 1, 2
                do 70 i = 1, row_A
                   A(i,j) = 0.0d0
  70            continue
  60         continue
             do 80 i = 1, row_A
                b(i) = 0.0d0
  80         continue
           if (LS_Dx .gt. 0.0d0) then
               A(1,1) = 1.0d0
            else
               A(1,1) = -1.0d0
            end if
            if (LS_Dy .gt. 0.0d0) then
               A(2,2) = 1.0d0
            else
               A(2,2) = -1.0d0
            end if            
           if (LS_Dx .gt. 0.0d0) then
              A(3,1) = -1.0d0
              b(3) = -LS_Dx
            else
               A(3,1) = 1.0d0
               b(3) = LS_Dx
            end if
           if (LS_Dy .gt. 0.0d0) then
              A(4,2) = -1.0d0
              b(4) = -LS_Dy
            else
               A(4,2) = 1.0d0
               b(4) = LS_Dy
            end if            

c            Add data constraints
             do 100 i = 1, num_neighb
                if ((neighb_u(i)-center_u).GT.0) then
                   A(2*i+3,1) = -(neighb_x(i)-center_x)
                   A(2*i+3,2) = -(neighb_y(i)-center_y)
                   A(2*i+4,1) = (neighb_x(i)-center_x)
                   A(2*i+4,2) = (neighb_y(i)-center_y)
                   b(2*i+3) = center_u-neighb_u(i) 
                   b(2*i+4) = 0.0d0 
               else
                   A(2*i+3,1) = (neighb_x(i)-center_x)
                   A(2*i+3,2) = (neighb_y(i)-center_y)
                   A(2*i+4,1) = -(neighb_x(i)-center_x)
                   A(2*i+4,2) = -(neighb_y(i)-center_y)
                   b(2*i+3) = neighb_u(i) - center_u
                   b(2*i+4) = 0.0d0
               endif
  100        continue

c          positivity constraint:
            if (do_pos_constr) then
                  A(row_A,1) = (bxpt-center_x)
                  A(row_A,2) = (bypt-center_y)
                  b(row_A) = - center_u
            endif

            x0(1) = 0.0d0
            x0(2) = 0.0d0
            c(1) = -sign(1.0d0,LS_Dx)
            c(2) = -sign(1.0d0,LS_Dy)
            x(1) = 0.0d0
            x(2) = 0.0d0
            W(1) = 1
            W(2) = 2

            temp_output = .FALSE.
            LP_success = .FALSE.
            
c          call 'Simplex_allineq' solver

           if (output_debug) then
            write(*,*) "\n########\n Call Simplex allineq\n#######\n"
           endif

           my_iters = 50  ! should be set to the correct value in Simplex_allineq if that routine was successful

           if (abs(LS_Dx)+abs(LS_Dy) .gt. 1.E-12) then
              call Simplex_allineq(row_A,A, b, c,
     &             x0,W,temp_output,x,LP_success,my_iters)
c     write(*,*) "LP_success",LP_success

              if (LP_success) then
                 LP_Dx = x(1)
                 LP_Dy = x(2)
c     write(*,902) "limited slopes", LP_Dx, LP_Dy
c     write(*,*) "\n ###################### \n"
            else
               LP_Dx = 0.0d0
               LP_Dy = 0.0d0
               my_iters = 0
            end if
            
            if (output_debug) then
               write(*,*)"LP_DX      LP_DY"
               write(*,903) LP_DX, LP_DY
               write(*,*) "center_x       center_y             center_u"
               write(*,903) center_x,center_y,center_u
               write(*,*)"neighb_x          neighb_y       neighb_u   ",
     &             "   predicted"
               do nn = 1, num_neighb
                  pred = center_u + (neighb_x(nn)-center_x)*LP_DX +
     &                              (neighb_y(nn)-center_y)*LP_DY 
                  write(*,903) neighb_x(nn),neighb_y(nn),neighb_u(nn),
     &                         pred
               end do
            endif

           else !nothing done, gradient came in (almost) zero
             LP_success = .TRUE.
             LP_Dx = LS_DX
             LP_Dy = LS_DY
             my_iters = 0
           endif


      end if    ! end of if (choice_limiting .eq. 1) then

            


      if (choice_limiting .eq. 2) then
      

c     #######################
c     Formulate 2d LP
c     #######################

         if (output_debug) then
            write(*,*) "############\n Formulate 2d LP \n ############"
         endif

c            ::::::::::::::::::::::
c            initialize A and b
c            ::::::::::::::::::::::

c            Initialize with 0
             do 62 j = 1, 2
                do 72 i = 1, row_A
                   A(i,j) = 0.0d0
 72             continue
 62          continue
             do 82 i = 1, row_A
                b(i) = 0.0d0
 82          continue
           if (LS_Dx .gt. 0.0d0) then
               A(1,1) = 1.0d0
            else
               A(1,1) = -1.0d0
            end if
            if (LS_Dy .gt. 0.0d0) then
               A(2,2) = 1.0d0
            else
               A(2,2) = -1.0d0
            end if            
           if (LS_Dx .gt. 0.0d0) then
              A(3,1) = -1.0d0
              b(3) = -LS_Dx
            else
               A(3,1) = 1.0d0
               b(3) = LS_Dx
            end if
           if (LS_Dy .gt. 0.0d0) then
              A(4,2) = -1.0d0
              b(4) = -LS_Dy
            else
               A(4,2) = 1.0d0
               b(4) = LS_Dy
            end if

            rhsmax = neighb_u(1) - center_u
            rhsmin = neighb_u(1) - center_u
            do 170 i=2,num_neighb
               rhsmax = dmax1(rhsmax,neighb_u(i) - center_u)
               rhsmin = dmin1(rhsmin,neighb_u(i) - center_u)
 170         continue
            rhsmax = dmax1(rhsmax,0.0d0)
            rhsmin = dmin1(rhsmin,0.0d0)            

c            Add data constraints
             do 102 i = 1, num_neighb
                   A(2*i+3,1) = -(neighb_x(i)-center_x)
                   A(2*i+3,2) = -(neighb_y(i)-center_y)
                   A(2*i+4,1) = (neighb_x(i)-center_x)
                   A(2*i+4,2) = (neighb_y(i)-center_y)
                   b(2*i+3) = - rhsmax
                   b(2*i+4) = rhsmin
 102        continue

c          positivity constraint:
            if (do_pos_constr) then
                  A(row_A,1) = (bxpt-center_x)
                  A(row_A,2) = (bypt-center_y)
                  b(row_A) = - center_u
            endif

            x0(1) = 0.0d0
            x0(2) = 0.0d0
            c(1) = -sign(1.0d0,LS_Dx)
            c(2) = -sign(1.0d0,LS_Dy)
            x(1) = 0.0d0
            x(2) = 0.0d0
            W(1) = 1
            W(2) = 2

            temp_output = .FALSE.
            LP_success = .FALSE.
            
c          call 'Simplex_allineq' solver

           if (output_debug) then
            write(*,*) "\n########\n Call Simplex allineq\n#######\n"
           endif

           my_iters = 50  ! should be set to the correct value in Simplex_allineq if that routine was successful

           if (abs(LS_Dx)+abs(LS_Dy) .gt. 1.E-12) then
              call Simplex_allineq(row_A,A, b, c,
     &             x0,W,temp_output,x,LP_success,my_iters)
c     write(*,*) "LP_success",LP_success

              if (LP_success) then
                 LP_Dx = x(1)
                 LP_Dy = x(2)
c     write(*,902) "limited slopes", LP_Dx, LP_Dy
c     write(*,*) "\n ###################### \n"
            else
               LP_Dx = 0.0d0
               LP_Dy = 0.0d0
               my_iters = 0
            endif
           else !nothing done, gradient came in (almost) zero
             LP_success = .TRUE.
             LP_Dx = LS_DX
             LP_Dy = LS_DY
             my_iters = 0
           endif
            
            if (output_debug) then
               write(*,*)"LP_DX      LP_DY"
               write(*,903) LP_DX, LP_DY
               write(*,*) "center_x       center_y             center_u"
               write(*,903) center_x,center_y,center_u
               write(*,*)"neighb_x          neighb_y       neighb_u   ",
     &             "   predicted"
               do nn = 1, num_neighb
                  pred = center_u + (neighb_x(nn)-center_x)*LP_DX +
     &                              (neighb_y(nn)-center_y)*LP_DY 
                  write(*,903) neighb_x(nn),neighb_y(nn),neighb_u(nn),
     &                         pred
               end do
            endif


            end if    ! end of if (choice_limiting .eq. 2) then

            









      return
      end
    

      

