module global_data

! Module to contain global model parameters

implicit none

real(8) :: alpha, beta, gamm, rho, rti, Rm0, t_lck, t_ulck   ! SEIR model parameters, check SEIR documentation
character(len=60), parameter :: outdir = '/out/', indir = '/input/'
character(len=600) :: cwd, infile, outfile

contains

	subroutine filenames
        
    implicit none
        
    CALL get_environment_variable('PWD',cwd) ! Access the current directory
	infile = trim(cwd)//trim(indir)//('inparams_SEIR')  ! The input parameters file
	write(*,*) infile
	outfile = trim(cwd)//trim(outdir)//('pred_out')  ! The output formatted file
    
    end subroutine filenames

end module global_data
