program seir_main

! Numerical recipes types
USE nrtype
USE nr_bsstep
USE global_data
USE SEIR_sub

implicit none

real(8), dimension(:), allocatable   :: st, et, it, rt, t_loc         ! The SEIR Time Series Predictions and time points
real(8), dimension(4)                :: yt, dydt, yt_scale, ytsv, dydtsv  ! Vectors for RKQS
real(8)                              :: dt_try, dt_next, time, timesv, dt_did
real(8)                              :: dt, t_tot, N_tot, s0, e0, r0
real(8)                              :: i0, rho_in, rho0, rtol, atol, acc
integer                              :: nt
integer                              :: i, j, ioerr1, ioerr2, irecout, sizeout(2), t_count
real                                 :: start, finish, time_tot

call filenames   ! Call the subroutine to assign all filenames
open(unit=10,file=infile,STATUS='OLD',IOSTAT = ioerr1) ! Open the text file and assign it a number


! Now we will read in the model characteristics one-by-one
read(10,*) alpha ! inverse of the incubation period (1/t_incubation)
read(10,*) beta  ! probability of infection from I to S
read(10,*) rti    ! rate of interaction between individuals
read(10,*) gamm  ! inverse of the mean infectious period (1/t_infectious)
read(10,*) rho_in  ! Effects of social distancing
read(10,*) t_tot ! How long do you want to run the code in days
read(10,*) N_tot ! Total population, we will solve ODEs in proportions
read(10,*) e0    ! Intial proportion of exposed
read(10,*) t_lck ! Time for introduction of lockdown
read(10,*) t_ulck ! Time for removing lockdown

! Calculate R0
Rm0 = beta/gamm
rho0 = 1.d0
rho = rho0

! Define initial values for the simulation
s0 = 1.d0 - 2.d0*e0    ! Intial proportion of susceptible
i0 = e0           ! Initial proportion of exposed
r0 = 0.d0         ! Initial proportion of recovered

close(10)         ! Close after reading

dt = 0.01d0   ! You are second order accurate in space

write(*,'(A15, F7.4)') 'Initial time step: ', dt

nt = ceiling(t_tot/dt)  ! Number of time steps
write(*,*) 'Total # of time steps :', nt 

! Now allocate the big variables
allocate(st(nt),et(nt), it(nt), rt(nt), t_loc(nt))
! Define time points
t_loc = 0.d0
do i=2,nt
	t_loc(i) = t_loc(i-1) + dt
enddo

! Now initial condition, so first find the edge of the dike
yt(1) = s0   ! y1 is st --> susceptible
yt(2) = e0   ! y2 is et --> exposed
yt(3) = i0   ! y3 is it --> infected
yt(4) = r0   ! y4 is rt --> recovered

! Open output file
open(unit=20,file=outfile,STATUS='REPLACE',IOSTAT = ioerr2) ! Open the text file and assign it a number
! Header: alpha, beta, gamma, rho, R0 --> First line of values are model parameters
write(20,'(5e16.9)') alpha, beta, gamm, rho_in, rti*beta/gamm
! First output: Format >>> t S(t) E(t) I(t) R(t)
write(20,'(5e16.9)') t_loc(1), N_tot*yt(1:4)

rtol   = 1.d-6 ! When 0.d0, pure absolute error control
atol   = 0.d0 ! When 0.d0, pure relative error control
acc    = 1.d-6
time   = t_loc(1)

do i = 2,nt

	dt_try = dt
	if (t_loc(i)>t_lck) then
		rho = rho_in
	elseif (t_loc(i)>t_ulck) then
		rho = rho0
	endif
!-----------------------------------------------------------------------------------
! Use num. recipes rqks
!-----------------------------------------------------------------------------------
	         do while (abs(time-t_loc(i)) >= acc)
	  	          call derivs(time,yt,dydt)
	!------------------------------------------------------------------------------------------------------------------
	! Define scaling according to Numerical Recipes trick in Chapter 16.2 
	!------------------------------------------------------------------------------------------------------------------
	        	  yt_scale = abs(yt)+abs(dt_try*dydt)
	!---------------------------------------------------------------------------------------------------------------------
	!---------------------------------------------------------------------------------------------------------------------
	  	          ytsv = yt
        	  	  dydtsv = dydt
	          	  timesv = time
    ! Use Runge-Kutta 4,5 integrator to integrate 
		          call rkqs(yt,dydt,time,dt_try,rtol,yt_scale,dt_did,dt_next,derivs)
	  	          if (time.gt.t_loc(i)) then
	  	        	yt   = ytsv
	  	        	dydt = dydtsv
	  	        	time = timesv
	  	         	dt_try = t_loc(i) - time  
	  	         else
	  	        	dt_try = dt_next
	  	         endif
             end do
      ! Time steps output: Format >>> t S(t) E(t) I(t) R(t)
	    write(20,'(5e16.9)') t_loc(i), N_tot*yt(1:4)
	    write(*,'(I8,A8,I8,A6)') i, ' out of ', nt, ' steps'
enddo

close(20)


end program seir_main
