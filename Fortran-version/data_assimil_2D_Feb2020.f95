program data_assimil_2D_Feb2020
	use FWandBWSolvers
	USE ParametersDA_Feb2020
	USE observations
	USE lib_array

implicit none


       
! Declared variables for kappa test
integer :: i,j,k,l, i2, j2, jtrial, iter ! Indices for Loops
integer :: Nt, temp_x, temp_y, index, Ntotal, Ny_temp, n_iter
integer, allocatable, dimension(:):: x0_inds, y0_inds
integer, allocatable, dimension(:,:) ::  obs_inds
double precision :: xmin, ymin, delta_x, delta_y, dt, x0_max, y0_max, amp
double precision, allocatable, dimension(:) :: x, y, T, x0_vals, y0_vals, x0_vals_all, y0_vals_all, x0_pts, y0_pts, T_obs, T_d, Tb
double precision, allocatable, dimension(:) :: cost_x, T_opt, cost_opt
double precision, allocatable, dimension(:,:) :: Xmg, Ymg, obs_pts, obs, eta_x0, eta_adj_t0, grad_J, eta_x0_opt
double precision, allocatable, dimension(:,:,:) :: sol_exct, obs_eta, sol_distorted, eta_all, u_all, v_all, Yb, sol_opt
double precision :: a, tau_n
double precision, dimension(Nx,Ny) :: eta_exct0, u, v, eta_optimum
double precision, dimension(Nx,Ny, iter_max) :: eta0, steepest_direction
double precision, dimension(3*Nx,Ny) :: IC, IC2, IC_bw, IC_iter, IC_iter_opt
double precision, dimension(2) :: tRange
double precision, dimension(iter_max+1, ntrial) :: cost, grad, error



xmin = -xmax 
ymin = -ymax

delta_x = abs(xmax - xmin)/Nx
delta_y = abs(ymax - ymin)/Ny
dt = min(delta_x,delta_y)/3
Nt = 1 + nint(( tmax - tmin)/dt )

allocate(x(Nx), y(Ny), T(Nt) )
x=xmin+delta_x*(/(i,i=0,Nx-1)/)
y=ymin+delta_y*(/(i,i=0,Ny-1)/)
T = tmin+dt*(/(i,i=0,Nt-1)/)

tRange(1) = tmin
tRange(2) = tmax

allocate(Xmg(Ny,Nx), Ymg(Nx,Ny))

Xmg = spread(x, 2, Ny)
Ymg = spread(y, 1, Nx)


!Creating values for x and y observation points
x0_max = x0_min+(n_obs_x-1)*dmu_x 
IF (x0_max > xmax) THEN
  dmu_x = (xmax-x0_min)/n_obs_x
  x0_max = x0_min+(n_obs_x-1)*dmu_x
  WRITE(*,*)  'X-distance between measurement points has been adjusted to fit domain'
ELSE
   x0_max = x0_min+(n_obs_x-1)*dmu_x
END IF

allocate(x0_vals(n_obs_x))
x0_vals = linspace_rk(x0_min, x0_max, n_obs_x)


y0_max = y0_min+(n_obs_y-1)*dmu_y
IF (y0_max > ymax) THEN
  dmu_y = (ymax-y0_min)/n_obs_y
  y0_max = y0_min+(n_obs_y-1)*dmu_y
  WRITE(*,*)  'Y-distance between measurement points has been adjusted to fit domain'
ELSE
   y0_max = y0_min+(n_obs_y-1)*dmu_y
END IF

allocate(y0_vals(n_obs_y))
y0_vals = linspace_rk(y0_min, y0_max, n_obs_y)


temp_x = size(x0_vals)*n_obs_y
temp_y = size(y0_vals)*n_obs_x
allocate(x0_vals_all(temp_x), y0_vals_all(temp_y))
x0_vals_all = reshape(spread(x0_vals, 2, n_obs_y), (/temp_x/))
y0_vals_all = reshape(spread(y0_vals, 1, n_obs_x), (/temp_y/))


IF (x_equals_y) THEN
	call mult_x0_y0(x, y, x0_vals, y0_vals, x0_inds, x0_pts, y0_inds, y0_pts )
	allocate( obs_pts(n_obs_x,2))
	obs_pts(:,1) = x0_pts
	obs_pts(:,2) = y0_pts
ELSE
	call mult_x0_y0(x, y, x0_vals_all, y0_vals_all, x0_inds, x0_pts, y0_inds, y0_pts )
	allocate( obs_pts(temp_x,2))
	obs_pts(:,1) = x0_pts
	obs_pts(:,2) = y0_pts
END IF

allocate( obs_inds(size(x0_inds),2) )
obs_inds(:,1) = x0_inds
obs_inds(:,2) = y0_inds

allocate( sol_exct(3*Nx,Ny,Nt) , sol_distorted(3*Nx,Ny,Nt), sol_opt(3*Nx,Ny,Nt), T_opt(Nt), T_obs(Nt), T_d(Nt) )
allocate( obs_eta(Nx, Ny, Nt), obs(size(x0_inds), size(T_obs)) )
allocate( eta_all(Nx,Ny,Nt), u_all(Nx,Ny,Nt), v_all(Nx,Ny,Nt), eta_x0(size(x0_inds), size(T_obs)) , &
		 eta_x0_opt(size(x0_inds), size(T_obs)) )
allocate( Yb(3*Nx,Ny,Nt), Tb(Nt), eta_adj_t0(Nx,Ny), grad_J(Nx,Ny) )
allocate(cost_x(Nt), cost_opt(Nt))

DO jtrial =1, ntrial
	print*, "Trial:", jtrial

	! STEP I: Run forward solver with exact IC to get observations at measurement points
	amp = 0.05D0
	eta_exct0 = amp*exp( - ( (Xmg)**2 + (Ymg)**2)/ 0.1**2)
	u = 0
	v = 0
	IC(1:Nx,:) = eta_exct0
	IC(Nx+1:2*Nx, :) = u
	IC(2*Nx+1:3*Nx, :) = v
	call FW_Solve_2D(dt, tRange, IC, shape(IC), delta_x, delta_y, sol_exct, T_obs)


	! Exact surface wave at all x,y points for all t
	obs_eta = sol_exct(1:Nx,:,:)
	! Exact wave at measurement points x0,y0 for all t

	DO k= 1, size(x0_inds)
	    obs(k,:) = sol_exct(x0_inds(k),y0_inds(k),:)
	end DO

	! STEP 2a: Distort exact IC to get "guess" IC eta0(:,1)
	! Run forward solver
	!eta0(:,:,1) = 0.9*eta_exct0;
	eta0(:,:,1) = 0 !0.5D0*amp*(sin(Xmg) + sin(Ymg))
	IC2(1:Nx,:) = eta0(:,:,1)
	IC2(Nx+1:2*Nx, :) = u
	IC2(2*Nx+1:3*Nx, :) = v
	call FW_Solve_2D(dt, tRange, IC2, shape(IC2), delta_x, delta_y, sol_distorted, T_d)

	eta_all = sol_distorted(1:Nx,:,:)
	u_all   = sol_distorted(Nx+1:2*Nx,:,:)
	v_all   = sol_distorted(2*Nx+1:3*Nx,:,:)

	DO j= 1, size(x0_inds)
	    eta_x0(j,:) = sol_distorted(x0_inds(j),y0_inds(j),:)
	end DO

	! STEP 2b: Run backward solver to get adjoint height at t0 for all x
	IC_bw = 0.D0
	call BW_Solve_2D(dt, tRange, IC_bw, shape(IC_bw), delta_x, delta_y, obs_eta, x0_inds, y0_inds, u_all, v_all, eta_all, &
	 			Yb, Tb, eta_adj_t0, Ntotal, Ny_temp )

	!Define grad_J
	grad_J = -eta_adj_t0


	    
	!Smooth gradient
	if smooth_grad
	    for k = 1:Ny
	        grad_J(:,k) = grad_smooth(grad_J(:,k),filt);
        end
        for l = 1:Nx
	        grad_J(l,:) = grad_smooth(grad_J(l,:),filt);
	    end
	end

	!Compute Cost Function
	cost_x = 0
	DO l=1, size(x0_inds)
		cost_x = cost_x + ( (obs(l,:) - eta_x0(l,:))**2) 
	END DO



	cost(1,jtrial) = integral(T,0.5*cost_x) ! Get trapz first
	grad(1,jtrial) = norm2(grad_J);

    	error(1,jtrial)  = norm2(eta_exct0-eta0(:,:,1))/norm2(eta_exct0)

    ! Initial estimate for gradient descent stepsize
  	tau_n = 1.D-2
  	iter=0

  	! BEGIN LOOP
    DO WHILE (norm2(grad_J) >= 1.D-20 .AND. iter<iter_max)
    	iter = iter+1
    	!print*, "iter =", iter
		! STEP 3
		! % Run forward solver to get height at x0 for all t given eta0_2
		IC_iter(1:Nx,:) = eta0(:,:,iter)
		IC_iter(Nx+1:2*Nx, :) = u
		IC_iter(2*Nx+1:3*Nx, :) = v

		call FW_Solve_2D(dt, tRange, IC_iter, shape(IC_iter), delta_x, delta_y, sol_distorted, T_d)
		eta_all = sol_distorted(1:Nx,:,:)
		u_all   = sol_distorted(Nx+1:2*Nx,:,:)
		v_all   = sol_distorted(2*Nx+1:3*Nx,:,:)

		DO j= 1, size(x0_inds)
	    	eta_x0(j,:) = sol_distorted(x0_inds(j),y0_inds(j),:)
		end DO

		IC_bw = 0.D0
	
		call BW_Solve_2D(dt, tRange, IC_bw, shape(IC_bw), delta_x, delta_y, obs_eta, x0_inds, y0_inds, u_all, v_all, eta_all, &
	 			Yb, Tb, eta_adj_t0, Ntotal, Ny_temp )

		grad_J = -eta_adj_t0

		!Smooth gradient
		!if smooth_grad
		!    for k = 1:Ny
		!        grad_J(:,k) = grad_smooth(grad_J(:,k),filt);
		!    end
		!    for l = 1:Nx
		!        grad_J(l,:) = grad_smooth(grad_J(l,:),filt);
		!    end
		!end

		steepest_direction(:,:,iter) = - grad_J

		! STEP 7: Steepest Descent Algorithm
        eta0(:,:,iter+1)= eta0(:,:,iter) + tau_n*steepest_direction(:,:,iter)
        eta_optimum = eta0(:,:,iter+1)

        ! STEP 8: Run forward solver to get height at x0 for all t given optimised IC
		IC_iter_opt(1:Nx,:) = eta0(:,:,iter+1)
		IC_iter_opt(Nx+1:2*Nx, :) = u
		IC_iter_opt(2*Nx+1:3*Nx, :) = v

        call FW_Solve_2D(dt, tRange, IC_iter_opt, shape(IC_iter_opt), delta_x, delta_y, sol_opt, T_opt)

        DO i2 = 1, size(x0_inds)
        	eta_x0_opt(i2,:) = sol_opt(x0_inds(i2),y0_inds(i2),:)
        END DO

        ! Compute Cost Function
        cost_opt = 0.D0
        DO j2= 1, size(x0_inds)
        	cost_opt = cost_opt + (obs(j2,:) - eta_x0_opt(j2,:))**2
        END DO


        cost(iter+1,jtrial) = integral(T,0.5*cost_opt) 
        grad(iter+1,jtrial) = norm2(grad_J)
        error(iter+1,jtrial) = norm2(eta_exct0-eta_optimum)/norm2(eta_exct0)


        print*, iter,"	", tau_n, norm2(grad_J),"	", norm2(eta_exct0-eta_optimum)/norm2(eta_exct0)
    
		! DEALLOCATE EVERTHING IN WHILE LOOP ********************
	END DO 

! DEALLOCATE everything in jtrial loop *************************
END do
n_iter = iter + 1

!err_std  = std(err(end,:));
!err      = mean(err,2);

!grad_std = std(grad(end,:));
!grad     = mean(grad,2);	

!cost_std = std(cost(end,:));
!cost     = mean(cost,2);
 
!print*, "N_obs = ", n_obs_x, ",mean error = ", error(ntrial), ",std dev = ", err_std(ntrial)
!print*, "this is the shape of grad   : ", shape(eta_optimum)
!print*, "this is the shape of grad   : ", shape(grad)
!print*, "this is the shape of obs_pts: ", shape(obs_pts)
!print*, "this is the shape of eta0   : ", shape(eta0)



open(unit = 2, file = "ntrial_file.dat", form = 'unformatted', access = 'stream')  
open(unit = 21, file = "res_xyt.dat", form = 'unformatted', access = 'stream')  
open(unit = 3, file = "n_obs_x_file.dat", form = 'unformatted', access = 'stream') 
open(unit = 4, file = "n_obs_y_file.dat", form = 'unformatted', access = 'stream') 
open(unit = 50, file = "rand_x0_file.dat", form = 'unformatted', access = 'stream') 
open(unit = 60, file = "x0_min_file.dat", form = 'unformatted', access = 'stream') 
open(unit = 61, file = "y0_min_file.dat", form = 'unformatted', access = 'stream') 
open(unit = 7, file = "dmu_x_file.dat", form = 'unformatted', access = 'stream') 
open(unit = 8, file = "dmu_y_file.dat", form = 'unformatted', access = 'stream') 
open(unit = 11, file = "iter_max_file.dat", form = 'unformatted', access = 'stream') 
open(unit = 12, file = "obs_vals_file.dat", form = 'unformatted', access = 'stream') 
open(unit = 13, file = "eta_optimum_file.dat", form = 'unformatted', access = 'stream') 
open(unit = 14, file = "eta_exct0_file.dat", form = 'unformatted', access = 'stream') 
open(unit = 15, file = "eta0_file.dat", form = 'unformatted', access = 'stream') 
open(unit = 16, file = "Xmg_file.dat", form = 'unformatted', access = 'stream') 
open(unit = 17, file = "Ymg_file.dat", form = 'unformatted', access = 'stream') 
open(unit = 18, file = "err_file.dat", form = 'unformatted', access = 'stream') 
open(unit = 19, file = "grad_file.dat", form = 'unformatted', access = 'stream') 
open(unit = 20, file = "tau_n_file.dat", form = 'unformatted', access = 'stream') 


write(2)  ntrial
write(21) Nx, Ny, Nt
write(3)  n_obs_x
write(4)  n_obs_y
write(50)  rand_x0
write(60)  x0_min
write(61)  y0_min
write(7)  dmu_x
write(8)  dmu_y
write(11)  iter_max
write(12)  obs_pts
write(13)  eta_optimum
write(14)  eta_exct0
write(15)  eta0
write(16)  Xmg
write(17)  Ymg
write(18)  error
write(19)  grad
write(20)  tau_n














! DON'T FORGET TO DEALLOCATE ALL THE DYNAMIC VARIABLES
deallocate( sol_exct(3*Nx,Ny,Nt) , sol_distorted(3*Nx,Ny,Nt), sol_opt(3*Nx,Ny,Nt), T_opt(Nt), T_obs(Nt), T_d(Nt) )
deallocate( obs_eta(Nx, Ny, Nt), obs(size(x0_inds), size(T_obs)) )
deallocate( eta_all(Nx,Ny,Nt), u_all(Nx,Ny,Nt), v_all(Nx,Ny,Nt), eta_x0(size(x0_inds), size(T_obs)) , &
		 eta_x0_opt(size(x0_inds), size(T_obs)) )
deallocate( Yb(3*Nx,Ny,Nt), Tb(Nt), eta_adj_t0(Nx,Ny), grad_J(Nx,Ny) )
deallocate(cost_x(Nt), cost_opt(Nt))
   
   
end program data_assimil_2D_Feb2020
