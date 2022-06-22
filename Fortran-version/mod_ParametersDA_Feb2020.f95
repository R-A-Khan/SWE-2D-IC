module ParametersDA_Feb2020

	! All the declared variables for the kappa test 

	implicit none

	integer, parameter :: n_obs_x              = 24          ! Number of observation points
	integer, parameter :: n_obs_y              = 24        ! Number of observation points
	logical :: x_equals_y  			  		   = .FALSE.      ! n_obs_x measurement points chosen along x = y
	logical :: rand_x0     			  		   = .FALSE.      ! Use random observation points
	double precision, parameter :: x0_min      = -2.4D0       ! x-Location of first observation point
	double precision, parameter :: y0_min      = -2.4D0       ! y-Location of first observation point
	integer, parameter :: ntrial      		   = 1          ! Number of trials
	integer, parameter :: Nx                   = 128        ! Number of y grid points
	integer, parameter :: Ny                   = 128        ! Number of y grid points
	double precision, parameter :: tmin        = 0.D0          ! Starting time
	double precision, parameter :: tmax        = 3.D0          ! control time
	double precision, parameter :: xmax        = 3.D0          ! domain size: xmin = -xmax
	double precision, parameter :: ymax        = 3.D0          ! domain size: ymin = -ymax
	double precision :: dmu_x 				   = 0.2D0        ! Spacing between x coord of observation points- can be adjusted
    double precision :: dmu_y    			   = 0.2D0    	 ! Spacing between x coord of observation points- can be adjusted
	integer, parameter :: iter_max    		   = 1000      ! Max number of iterations;    
	double precision, parameter :: cfl         = 1.D0          ! CFL constant
	double precision, parameter :: nu          = 0.D0          ! Kinematic viscosity co-efficient                 
	logical :: smooth_grad 			 		   = .FALSE.
	double precision, parameter :: filt        = 0.9D0
	

end module 