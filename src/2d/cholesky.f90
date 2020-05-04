! cholesky_d.f -*-f90-*-
! Using Cholesky decomposition, cholesky_d.f solve a linear equation Ax=b,
! where A is a n by n positive definite real symmetric matrix, x and b are
! real*8 vectors length n.
!
! Time-stamp: <2015-06-25 18:05:47 takeshi>
! Author: Takeshi NISHIMATSU
! Licence: GPLv3
!
! [1] A = G tG, where G is a lower triangular matrix and tG is transpose of G.
! [2] Solve  Gy=b with  forward elimination
! [3] Solve tGx=y with backward elimination
!
! Reference: Taketomo MITSUI: Solvers for linear equations [in Japanese]
!            http://www2.math.human.nagoya-u.ac.jp/~mitsui/syllabi/sis/info_math4_chap2.pdf
!
! Comment:   This Cholesky decomposition is used in src/elastic.F and
!            src/optimize-inho-strain.F of feram http://loto.sourceforge.net/feram/ .
!!
      subroutine cholesky(nn, n, A, G)
      implicit double precision(a-h,o-z)
      dimension A(nn,nn)
      dimension G(nn,nn)

      ! Light check of positive definite


      ! [1]
      G(:,:)=0.0d0
      do j = 1, n
        G(j,j) = sqrt( max(1.d-14, A(j,j) - dot_product(G(j,1:j-1),G(j,1:j-1))))
        do i = j+1, n
          G(i,j) = (A(i,j)-dot_product(G(i,1:j-1),G(j,1:j-1)) ) / G(j,j)
        end do

      end do

      end subroutine
