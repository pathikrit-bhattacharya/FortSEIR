module SEIR_sub

! Numerical Recipes Typecast

use nrtype ! Numerical Recipes type definitions
use global_data ! Paramter definitions
use omp_lib ! OMP runtime directives

contains

!--------------------------------------------------------------------------------------------!
! Subroutine for defining ageclass model time derivatives
!--------------------------------------------------------------------------------------------!

SUBROUTINE derivs_agc(t,y,dydt)

USE nrtype

implicit none

REAL(8), INTENT(IN)                :: t
REAL(8), DIMENSION(:), INTENT(IN)  :: y
REAL(8), DIMENSION(:), INTENT(OUT) :: dydt
real(8)                            :: locdown_func
real(8), DIMENSION(nac,nac)        :: C_sl, C_wl, C_hl, C_ol
integer                            :: i

C_tot = 0.d0

if (homo_lockdown == 1) then
!--------------------------------------------------------------------------------------------------!
! No age structure in lockdown, all contact matrices are affected equally
!--------------------------------------------------------------------------------------------------!
	locdown_func  = 1.d0 - 0.5d0*(tanh((t-t_lck)/0.1d0) - tanh((t-t_ulck)/0.1d0))
	C_tot = locdown_func*C_s + C_h + locdown_func*C_o + locdown_func*C_w

elseif (homo_lockdown == 2) then
!--------------------------------------------------------------------------------------------------!
! Modify contact matrices for age-dependent effect of lockdown scenarios
!--------------------------------------------------------------------------------------------------!
	if ((t>t_lck).and.(t<t_ulck)) then
		do i = 1,nac
			C_sl(:,i) = w_s(i)*C_s(:,i)   ! Age-group dependent modification at school
			C_wl(:,i) = w_w(i)*C_w(:,i)   ! Age-group dependent modification at work
			C_hl(:,i) = w_h(i)*C_h(:,i)   ! Age-group dependent modification at home
			C_ol(:,i) = w_o(i)*C_o(:,i)   ! Age-group dependent modification at other locations
		enddo
	else
		C_sl(:,:) = C_s(:,:)
		C_wl(:,:) = C_w(:,:)
		C_hl(:,:) = C_h(:,:)
		C_ol(:,:) = C_o(:,:)
	endif
	C_tot = C_sl + C_hl + C_ol + C_wl
endif	
! function for defining derivatives for SEIR model
! y: Independent variables, the time series
! t: Time
! dydt: derivatives
! All proportions with respect to initial population
! yt(1:neq:ncmp) --> Intial proportion of susceptible - S_a's
! yt(2:neq:ncmp) --> Initial # of exposed - E_a's
! yt(3:neq:ncmp) --> Initial # of infected and symptomatic - I_ac's
! yt(4:neq:ncmp) --> Initial # of infected but asymptomatic - I_asc's
! yt(5:neq:ncmp) --> Initial # of recovered - R_a's
! yt(6:neq:ncmp) --> Total proportion - R_a's

dydt(1:neq:ncmp) = lambda*y(6:neq:ncmp) - alpha*beta*y(1:neq:ncmp)*matvec_MKL(C_tot,y(3:neq:ncmp)/Nis) &
  - beta*y(1:neq:ncmp)*matvec_MKL(C_tot,y(4:neq:ncmp)/Nis) - mu_n(:)*y(1:neq:ncmp)
dydt(2:neq:ncmp) = alpha*beta*y(1:neq:ncmp)*matvec_MKL(C_tot,y(3:neq:ncmp)/Nis) + beta*y(1:neq:ncmp)*matvec_MKL(C_tot,y(4:neq:ncmp)/Nis) &
 - kappa(:)*y(2:neq:ncmp) - mu_d(:)*y(2:neq:ncmp)
dydt(3:neq:ncmp) = rho(:)*kappa(:)*y(2:neq:ncmp) - (gamm(:)+mu_d(:))*y(3:neq:ncmp)
dydt(4:neq:ncmp) = (1.d0-rho(:))*kappa(:)*y(2:neq:ncmp) - (gamm(:)+mu_d(:))*y(4:neq:ncmp)
dydt(5:neq:ncmp) = gamm(:)*y(3:neq:ncmp) + gamm(:)*y(4:neq:ncmp) - mu_n(:)*y(5:neq:ncmp)
dydt(6:neq:ncmp) = lambda*y(6:neq:ncmp) - mu_n(:)*(y(1:neq:ncmp)+y(5:neq:ncmp)) - mu_d(:)*(y(2:neq:ncmp)+y(3:neq:ncmp)+y(4:neq:ncmp))

END SUBROUTINE derivs_agc

!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!
!!! Multithreaded matrix vector multiplication using Intel MKL library function dgemm
!!!
 function matvec_MKL(Amat,v) result(y)

	implicit none

 real(8), dimension(:,:), intent(in) :: Amat ! Input matrix for the matrix vector product A*v
 real(8), dimension(:), intent(in)   :: v
 real(8), dimension(size(Amat,1))    :: y ! Output vector, y = A*v
 real(8)                             :: alph, bet ! scalars, dgemm allows for the product C = alpha*A*B + beta*C
 integer 			                 :: m, k, n

 if (size(Amat,2) /= size(v)) then
	write(*,'(A,I4,A,I4)') 'size(Amat,2) = ', size(Amat,2), ', size(v) = ', size(v)
 	stop 'Matrix-vector multiplication rule violated in C=A*v, ncolA .neq. nrowv' ! Check for valid multiplication
 end if
 m = size(Amat,1); n = size(Amat,2) ! Assigning size values
 !----------------------------------------------------------------------------!
 ! DGEMV calculates y = alpha*A*x + beta*y
 ! call dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
 !----------------------------------------------------------------------------!
 alph = 1.d0 ! Unscaled product
 bet  = 0.d0 ! Nothing is being added
 !----------------------------------------------------------------------------!
 CALL DGEMV('N',m,n,alph,Amat,m,v,1,bet,y,1)
 
 end function matvec_MKL


!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!
!!! Printing a matrix sometimes to see whats happening in the code, taken from intel MKL, old style F77
!!!
	  SUBROUTINE print_matrix(DESC, M, N, A)
		
		implicit none

      CHARACTER(len=*)           :: DESC
      INTEGER                    :: M, N, LDA
	  REAL(8)                    :: A(m,n)
	  CHARACTER(len=5)           :: ncols
	  CHARACTER(len=15)          :: FMT
      INTEGER                    :: I, J
      WRITE(*,*) DESC
	  WRITE(ncols,'(I5)') N
	  FMT = ('(')//trim(ncols)//('d13.5)')
      DO I = 1, M
         WRITE(*,adjustl(adjustr(FMT))) (A(I,J), J = 1, N)
      END DO
        
	  END subroutine PRINT_MATRIX

!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!
!!! Output the T (Tm) and Sigma (Sm) matrices for the calculation of R_0
!!!

	  subroutine r0_out(Tm, Sm, a, g, k, r, mu, b, nf, C)

		implicit none
		real(8), dimension(:,:), intent(inout) :: Tm, Sm
		real(8), dimension(:), intent(in)      :: g, k, r, mu, nf
		real(8), dimension(:,:), intent(in)    :: C
		real(8)                                :: a, b
		integer :: i, j

		Tm = 0.d0
		Sm = 0.d0

!-------------------------------------------------------------------------------------------------------!
 ! Block 1 - 1:n, 1:n
!-------------------------------------------------------------------------------------------------------!
		do i = 1,nac
			do j = 1, nac
				Tm(i,j) = 0.d0
				if (i.eq.j) then
					Sm(i,j) = - k(i) - mu(i)
				else
					Sm(i,j) = 0.d0
				endif
			enddo
		enddo

!-------------------------------------------------------------------------------------------------------!
 ! Block 2 - 1:n, n+1:2*n
!-------------------------------------------------------------------------------------------------------!
		do i = 1,nac
			do j = nac+1, 2*nac
				Tm(i,j) = a*b*C(i,j-nac)*nf(i)/nf(j-nac)
				Sm(i,j) = 0.d0
			enddo
		enddo

!-------------------------------------------------------------------------------------------------------!
 ! Block 3 - 1:n, 2*n+1:3*n 
!-------------------------------------------------------------------------------------------------------!
		do i = 1,nac
			do j = 2*nac+1, 3*nac
				Tm(i,j) = b*C(i,j-2*nac)*nf(i)/nf(j-2*nac)
				Sm(i,j) = 0.d0
			enddo
		enddo

!-------------------------------------------------------------------------------------------------------!
 ! Block 4 - n+1:2*n, 1:n
!-------------------------------------------------------------------------------------------------------!
		do i = nac+1,2*nac
			do j = 1, nac
				Tm(i,j) = 0.d0
				if (i-nac.eq.j) then
					Sm(i,j) = r(i-nac)*k(i-nac)
				else
					Sm(i,j) = 0.d0
				endif
			enddo
		enddo

!-------------------------------------------------------------------------------------------------------!
 ! Block 5 - n+1:2*n, n+1:2*n
!-------------------------------------------------------------------------------------------------------!
		do i = nac+1,2*nac
			do j = nac+1,2*nac
				Tm(i,j) = 0.d0
				if (i-nac.eq.j-nac) then
					Sm(i,j) = - g(i-nac) - mu(i-nac)
				else
					Sm(i,j) = 0.d0
				endif
			enddo
		enddo

!-------------------------------------------------------------------------------------------------------!
 ! Block 6 - n+1:2*n, 2*n+1:3*n 
!-------------------------------------------------------------------------------------------------------!
		do i = nac+1,2*nac
			do j = 2*nac+1, 3*nac
				Tm(i,j) = 0.d0
				Sm(i,j) = 0.d0
			enddo
		enddo

!-------------------------------------------------------------------------------------------------------!
 ! Block 7 - 2*n+1:3*n, 1:n
!-------------------------------------------------------------------------------------------------------!
		do i = 2*nac+1,3*nac
			do j = 1, nac
				Tm(i,j) = 0.d0
				if (i-2*nac.eq.j) then
					Sm(i,j) = (1-r(i-2*nac))*k(i-2*nac)
				else
					Sm(i,j) = 0.d0
				endif
			enddo
		enddo

!-------------------------------------------------------------------------------------------------------!
 ! Block 8 - 2*n+1:3*n, n+1:2*n
!-------------------------------------------------------------------------------------------------------!
		do i = 2*nac+1,3*nac
			do j = nac+1,2*nac
				Tm(i,j) = 0.d0
				Sm(i,j) = 0.d0
			enddo
		enddo

!-------------------------------------------------------------------------------------------------------!
 ! Block 9 - n+1:2*n, 2*n+1:3*n 
!-------------------------------------------------------------------------------------------------------!
		do i = 2*nac+1,3*nac
			do j = 2*nac+1, 3*nac
				Tm(i,j) = 0.d0
				if (i-2*nac.eq.j-2*nac) then
					Sm(i,j) = -g(i-2*nac) - mu(i-2*nac)
				else
					Sm(i,j) = 0.d0
				endif
			enddo
		enddo

		! call print_matrix('S_matrix : ',3*nac,3*nac,Tm)


	end subroutine r0_out

end module SEIR_sub
