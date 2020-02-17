c     ##########################################################
c     # 2d LP solver based on Simplex algorithm                #
c     #        for directional limiting of cutcells            #
c     # This solver is based on the allineq form of Simplex    #
c     ##########################################################


      subroutine Simplex_allineq(row_A,A, b, c,
     &     x0,W,output_debug,x,LP_success,my_iters)



c :::::::::::::::::::::: Simplex_allineq :::::::::::::::::::::::::
c
c  INPUT:
c     * row_A: no of rows of matrix A || row_A = number of constraints
c     * A: constraint matrix of LP: A(row_A,2)
c     * b: RHS of LP 
c     * c: vector in objective in LP
c     * x0: admissible starting point for LP
c     * W: working set for Simplex method
c     * output_debug: logical whether we should output or not
c  
c  OUTPUT:
c     * x: solution of LP
c     * LP_success: was the Simplex algorithm successful?
c     * my_iters: number of iterations needed
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     ::::::::::::::::::::::::::::
c     Declare variables
c     ::::::::::::::::::::::::::::

      implicit none  ! force yourself to properly declare all var

c     Input
      integer row_A
      double precision A(row_A,2), b(row_A), c(2) 
      double precision x0(2)  
      integer W(2)
      logical output_debug
      integer my_fort

c     Output
      double precision x(2)
      logical LP_success
      integer my_iters ! number of iterations needed

c     Local variables: Optimization variables
c      double precision matr_B(row_A,row_A), B_trans(row_A,row_A) ! columns of A corresponding to index_B
c      double precision matr_N(row_A,2) ! columns of A corresponding to index_N
c      double precision x_B(row_A), c_B(row_A)  ! components of x,c in set B
c      double precision x_N(2), c_N(2)  ! components of x,c in set N
c      double precision pi_vec(row_A)  ! vector pi of Lagrange multipliers
c      double precision s_N(2)  ! Lagrange multiplier s for set N

c     Local variables: other stuff
      integer i,j,k, shift_int, index_min_alpha
      double precision A_k(2,2)
      double precision det, lambda(2),min_lambda
      double precision e_q(2),p_k(2)
      double precision gamma,alpha,test
      integer min_index_W,size_D_k
c      integer ind_perm(row_A)  ! permutation array coming back from solving matrix
c      integer INFO ! was subroutine of LAPACK successful?
c      double precision product, max_step_dir, min_comp ! help variables
c      integer index_N_to_B ! index in set N chosen to go to set B
c      integer index_B_to_N ! index in set B chosen to go to set N (ie limits step length)
c      double precision rhs_for_dgesv(row_A) ! help vector
c      double precision step_dir(row_A) ! help vector
c      integer min_comp_index_B, min_index_s_N
c      double precision test_admiss(row_A),max_admiss ! test whether B*x_B-b = 0



 902     format (2(E12.2))
 901     format (A,'\t',F6.3)
 908     format (A,'\t',8(F6.3,'\t'))
 910     format (A,'\t',10(F6.3,'\t'))
 912     format (A,'\t',12(F6.3,'\t'))
 918     format (A,'\t',8(E12.3,'\t'))

c     ::::::::::::::::::::::::::::
c     Print given data
c     ::::::::::::::::::::::::::::

      if (output_debug) then
         write(*,*) "###"
         write(*,*) "Subroutine Simplex_allineq"
         write(*,*) "###"
      endif

      if (output_debug) then
         write (*,*) "%%%%%%%"
         write (*,*) "Given data for LP"
         write (*,*) "%%%%%%%"
         write(*,*) "number of constraints", row_A
         write(*,*) "matrix A"
         do 10 i = 1, row_A
            write(*,902) (A(i,j),j=1,2)
 10      continue
c
         if (row_A .EQ. 10) then
             write(*,910) "b",b
         elseif (row_A .EQ. 12) then
             write(*,912) "b",b
         endif
         write(*,*) "x0:", x0(1), x0(2)
         write(*,*) "W:", W
         write(*,*) "c:", c(1), c(2)
      endif   ! end of  if (output_debug) 


c     #########################
c     Beginning of loop
c     #########################


      do 555 k = 1, 20
      
          if (output_debug) then
              write(*,*) "%%%%%%%%"
              write(*,*) "Start iteration ",k
              write(*,*) "%%%%%%%%"
              write(*,*) "x:", x0(1), x0(2)
              write(*,*) "W:", W
              write(*,*) "---"
          endif
            

c      Step 1 of Algorithm -----------------------------
c      solve for Lagrange multiplier lambda: A_k^T lambda = c
          A_k(1,1) = A(W(1),1)
          A_k(1,2) = A(W(1),2)
          A_k(2,1) = A(W(2),1)
          A_k(2,2) = A(W(2),2)
          det = A_k(1,1)*A_k(2,2)-A_k(1,2)*A_k(2,1)
          ! added mjb for robustness
          if (det .eq. 0.0) then
             LP_Success = .false.
             return
          endif          
          lambda(1) = (1.0d0/det)*(A_k(2,2)*c(1)-A_k(2,1)*c(2))
          lambda(2) = (1.0d0/det)*(-A_k(1,2)*c(1)+A_k(1,1)*c(2))
          if (output_debug) then
             write(*,*) "lambda", lambda
          endif



c     Step 2 of Algorithm --------------------------------
c     Test whether lambda >= 0
         min_lambda = min(lambda(1),lambda(2))
         if (min_lambda .GT. -1.E-12) then
c            write(*,*) "\n ###################### \n"
c            write(*,*) "Have found the optimal x\n"
            x(1) = x0(1)
            x(2) = x0(2)
c
c            if ((x(1) > 1.00005) .AND. (x(2) > 1.00005)) then
c               write(*,*) "both slope > 1", x(1), x(2)
c            endif
            LP_success = .TRUE.
c          write(*,*) "LP_success", LP_success
c            write(25,*) "No-steps", k-1
c            write(25,*) "x", x
            my_iters = k-1
            return
         endif




c         Step 3 of Algorithm -----------------------------------
c         select q with (lambda)_q < 0 as index to leave the working set
c         have 2 choices
c         our choice for now: go with the most negative one
          if (lambda(1).GE.lambda(2)) then
            min_index_W = 2
          else
            min_index_W = 1
          endif
          if (output_debug) then
            write(*,*) "Index which leaves W",W(min_index_W)
          endif
c         or go with the first negative index
c         if (lambda(1) < -1.E-12)
c             min_index_W = 1
c         else
c             min_index_W = 2
c         endif


c         Step 4 of Algorithm ---------------------------------
c         Calculate step direction by solving A_k p_k = e_q
          if (min_index_W .EQ. 1) then
             e_q = (/ 1.0d0, 0.0d0 /)
          else
             e_q = (/ 0.0d0, 1.0d0 /)
          endif
          p_k(1) = (1.0d0/det)*(A_k(2,2)*e_q(1)-A_k(1,2)*e_q(2))
          p_k(2) = (1.0d0/det)*(-A_k(2,1)*e_q(1)+A_k(1,1)*e_q(2))
          if (output_debug) then
             write(*,*) "p_k", p_k
          endif



c         Step 5 of Algorithm ----------------------------------
c         Calculate step length
          size_D_k = 0
          alpha = 100000000.0d0
          index_min_alpha = row_A + 10
          do 20 i=1,row_A
             test = A(i,1)*p_k(1) + A(i,2)*p_k(2)
             if (test .LT. -1.E-12) then
                gamma = -(A(i,1)*x0(1) + A(i,2)*x0(2) - b(i))/test
                if (gamma .LT. (alpha-1.E-12)) then
c               Make sure that you don't take an index in W which is already there
c               Contribution for indices in W should be 0 -- but by numerical error it
c               could happen that the error is bigger than the tested 1.E-12
                   if ((i .ne. W(1)) .AND. (i .ne. W(2))) then
                      size_D_k = size_D_k + 1
                      alpha = gamma
                      index_min_alpha = i
                   endif
                endif
             endif
 20       continue
          if (size_D_k .EQ. 0) then
             write(*,*)  "\n LP unbounded !!!!!!!!!!!!\n"
             LP_success = .FALSE.
             return
          endif
          if (index_min_alpha .GT. row_A) then
             write(*,*)  "\n Problem with index_min_alpha\n"
             LP_success = .FALSE.
             return
          endif
          if (output_debug) then
             write(*,*) "Index which enters W",index_min_alpha
          endif
          ! update x0
          x0(1) = x0(1) + alpha*p_k(1)
          x0(2) = x0(2) + alpha*p_k(2)
          if (output_debug) then
             write(*,*) "new iterate x0", x0
          endif


c         Step 6 of Algorithm ----------------------------------
c         Update working set
c         get rid of the index in position "min_index_W" and put
c         index_min_alpha there
          W(min_index_W) = index_min_alpha



c         write(*,*) "After ",k," steps:", "Phi_x, Phi_y =",x0(1),x0(2)


 555  continue    ! end of main loop

c     #########################
c     End of loop
c     #########################


            x(1) = x0(1)
            x(2) = x0(2)

            my_iters = k-1


      return
      end

