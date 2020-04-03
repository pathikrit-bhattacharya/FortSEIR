module util_mkl

use lapack95    ! Shipped with Intel MKL, F95 interfaces to the original F77 subroutines
use f95_precision ! Shipped with Intel MKL


IMPLICIT NONE

! If any interfaces to add, here

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Solve a square linear system using dgesv from LAPACK
!!! Only for use when different RHS require different operators, otherwise look at ludcmp_MKL below
!!!
subroutine linsolve(Amat,solvec, info)
 
 implicit none

 real(8), dimension(:,:), intent(inout)   :: Amat ! Input matrix, replaced on output by LU factors
 real(8), dimension(:), intent(inout)     :: solvec ! Input vector b in Ax=b, replaced on output by solution x
 integer, intent(inout)                   :: info
 integer, dimension(:), allocatable       :: ipiv   ! Vector for storing pivot locations
 integer :: i, n
! Similar to the NR version, check if square
	if (size(Amat,2) /= size(Amat,1)) then
       		stop 'size error in sqrtm' ! Check for square matrix
    	end if    
	n = size(Amat,1) ! Amat is square
	allocate(ipiv(n))  ! All of the pivot indices
! Calling LAPACK routine dgesv, the LU factor solver from 
	call dgesv(n, 1, Amat, n, ipiv, solvec, n, info)
	!     Check for the exact singularity.
        if ( info.gt.0 ) then
         	write(*,*)'Warning! Matrix is singular to working precision, not computed'
			solvec = 0.d0
         	! stop
        end if
 end subroutine linsolve
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Solve for the LU factors of a square matrix using dgetrf (or getrf in F90 standard) from LAPACK
!!! Only for use when different RHS require different operators, otherwise look at ludcmp_MKL below
!!!

subroutine ludcmp_mkl(Amat,ipiv, info)
 
 implicit none
 
 real(8), dimension(:,:), intent(inout)   :: Amat ! Input matrix, replaced on output by LU factors
 integer, dimension(:), intent(out)       :: ipiv ! Vector for storing pivots
 integer, intent(inout)                   :: info
 integer :: i, n
! Similar to the NR version, check if square
	if (size(Amat,2) /= size(Amat,1)) then
       		stop 'size error in sqrtm' ! Check for square matrix
    end if    
	n = size(Amat,1) ! Amat is square

! Calling LAPACK routine dgetrf, the LU factors will be stored in Amat
!   call dgetrf( m, n, a, lda, ipiv, info)
!	call dgetrf(n, n, Amat, n, ipiv, info) ! Fortran 77 standard
!   call getrf( a [,ipiv] [,info] )  ! All square brackets imply optional
	call getrf(Amat,ipiv,info) ! Fortran 95, overloaded subroutine
	!     Check for the exact singularity.
        if ( info.gt.0 ) then
         	write(*,*)'Warning! Matrix is singular to working precision, not computed'
         	stop
        end if
 end subroutine ludcmp_mkl
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Solve for a linear system of equations with the LU factors of a square matrix using dgetrf from LAPACK
!!! Uses the LAPACK routine dgetrs (or getrs in Fortran 90 standard)
!!! Only for use when different RHS require different operators, otherwise look at ludcmp_MKL below
!!!

subroutine lubksb_mkl(Amat,ipiv,solvec,info)
 
 implicit none
 
 real(8), dimension(:,:), intent(in)   :: Amat ! Input matrix, contains the permutated LU factor from ludcmp_mkl
 real(8), dimension(:), intent(inout)  :: solvec ! Input vector b in Ax=b, replaced on output by solution x
 integer, dimension(:), intent(in)     :: ipiv ! Stores the locations of the pivot elements as obtained by ludcmp_mkl
 integer, intent(inout)                :: info
 integer :: i, n
! Assumes that will only use output from ludcmp_mkl
n = size(Amat,1) ! Amat is square
! Calling LAPACK routine dgetrs, uses the LU factors in Amat
! dgetrs( trans, n, nrhs, a, lda, ipiv, b, ldb, info ), F77
!	call dgetrs('N', n, 1, Amat, n, ipiv, solvec, n, info)
! getrs( a, ipiv, b [, trans] [,info] ), F95, default value of trans = 'N'
	call getrs(Amat, ipiv, solvec, 'N', info)
!	if ( info.gt.0 ) then
!       write(*,*)'Warning! Matrix is singular to working precision, not computed'
!		solvec = 0.d0
!       stop
!   end if

 end subroutine lubksb_mkl

!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!
!!! Multithreaded matrix multiplication using Intel MKL library function dgemm
!!!
 function matmult_MKL(Amat,Bmat) result(C)

 real(8), dimension(:,:), intent(in)           :: Amat, Bmat ! Input matrices for the product A*B
 real(8), dimension(size(Amat,1),size(Bmat,2)) :: C ! Output matrix, C = A*B
 real(8)                                       :: alpha, beta ! scalars, dgemm allows for the product C = alpha*A*B + beta*C
 integer 			               :: m, k, n

 if (size(Amat,2) /= size(Bmat,1)) then
	write(*,'(A,I4,A,I4)') 'size(Amat,2) = ', size(Amat,2), ', size(Bmat,1) = ', size(Bmat,1)
 	stop 'Matrix multiplication rule violated in C=A*B, ncolA .neq. nrowB' ! Check for valid multiplication
 end if
 m = size(Amat,1); k = size(Amat,2); n = size(Bmat,2) ! Assigning size values
 !----------------------------------------------------------------------------!
 ! DGEMM calculates C = alpha*A*B + beta*C
 !----------------------------------------------------------------------------!
 alpha = 1.d0 ! Unscaled product
 beta  = 0.d0 ! Nothing is being added
 !----------------------------------------------------------------------------!
 CALL DGEMM('N','N',m,n,k,alpha,Amat,m,Bmat,k,beta,C,m)
 
 end function matmult_MKL

!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!
!!! Multithreaded matrix vector multiplication using Intel MKL library function dgemm
!!!
 function matvec_MKL(Amat,v) result(y)

 real(8), dimension(:,:), intent(in) :: Amat ! Input matrix for the matrix vector product A*v
 real(8), dimension(:), intent(in)   :: v
 real(8), dimension(size(Amat,1))    :: y ! Output vector, y = A*v
 real(8)                             :: alpha, beta ! scalars, dgemm allows for the product C = alpha*A*B + beta*C
 integer 			     :: m, k, n

 if (size(Amat,2) /= size(v)) then
	write(*,'(A,I4,A,I4)') 'size(Amat,2) = ', size(Amat,2), ', size(v) = ', size(v)
 	stop 'Matrix-vector multiplication rule violated in C=A*v, ncolA .neq. nrowv' ! Check for valid multiplication
 end if
 m = size(Amat,1); k = size(Amat,2); n = 1 ! Assigning size values
 !----------------------------------------------------------------------------!
 ! DGEMV calculates y = alpha*A*x + beta*y
 ! call dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
 !----------------------------------------------------------------------------!
 alpha = 1.d0 ! Unscaled product
 beta  = 0.d0 ! Nothing is being added
 !----------------------------------------------------------------------------!
 CALL DGEMV('N',m,n,alpha,Amat,m,v,1,beta,y,1)
 
 end function matvec_MKL
 


!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!
end module util_mkl
