module SEIR_sub

! Numerical Recipes Typecast

use nrtype ! Numerical Recipes type definitions
use global_data ! Paramter definitions

contains

!-----------------------------------------------------------------------------------------!
! Subroutine for defining time derivatives
!-----------------------------------------------------------------------------------------!


SUBROUTINE derivs(t,y,dydt)

USE nrtype

implicit none

REAL(8), INTENT(IN) :: t
REAL(8), DIMENSION(:), INTENT(IN) :: y
REAL(8), DIMENSION(:), INTENT(OUT) :: dydt

			
! function for defining derivatives for SEIR model
! y: Independent variables, the time series
! t: Time
! dydt: derivatives



dydt(1) = -rho*rti*beta*y(1)*y(3)
dydt(2) = rho*beta*rti*y(1)*y(3) - alpha*y(2)
dydt(3) = alpha*y(2) - gamm*y(3) 
dydt(4) = gamm*y(3)

END SUBROUTINE derivs




end module SEIR_sub
