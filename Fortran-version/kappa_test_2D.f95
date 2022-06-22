program kappa_test_2D
	use FWandBWSolvers
	USE ParametersDA
	USE observations
	USE lib_array

implicit none


       
! Declared variables for kappa test
integer :: i,j,k,l,i2 ,j2, k2, l2!indexing for loops
integer :: Nt, temp_x, temp_y, index, Ntotal, Ny_temp
integer, allocatable, dimension(:):: x0_inds, y0_inds
integer, allocatable, dimension(:,:) ::  obs_inds
double precision :: xmin, ymin, delta_x, delta_y, dt, x0_max, y0_max, amp, J_eta, temprow2
double precision, allocatable, dimension(:) :: x, y, T, x0_vals, y0_vals, x0_vals_all, y0_vals_all, x0_pts, y0_pts, T_obs, T_d, Tb
double precision, allocatable, dimension(:) :: Tpert, sumJ, sumP, Gateaux_grad, Riesz_grad, temprow, kappa
double precision, allocatable, dimension(:,:) :: Xmg, Ymg, obs_pts, obs, eta_x0, eta_adj_t0, grad_J, etapert_x0_t_mult
double precision, allocatable, dimension(:,:,:) :: sol_exct, obs_eta, sol_distorted, eta_all, u_all, v_all, Yb, solpert
double precision :: a
double precision, dimension(Nx,Ny) :: eta_exct0, u, v, eta0_prime_g, eta_pert, u_p, v_p, riesz_int
double precision, dimension(Nx,Ny, iter_max) :: eta0
double precision, dimension(3*Nx,Ny) :: IC, IC2, IC_bw, IC_pert
double precision, dimension(2) :: tRange
double precision, dimension(10) :: ep


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


! STEP I: Run forward solver with exact IC to get observations at measurement points
amp = 0.05D0
eta_exct0 = amp*exp( - ( (Xmg)**2 + (Ymg)**2)/ 0.1**2)
u = 0
v = 0
IC(1:Nx,:) = eta_exct0
IC(Nx+1:2*Nx, :) = u
IC(2*Nx+1:3*Nx, :) = v

allocate( sol_exct(3*Nx,Ny,Nt) , sol_distorted(3*Nx,Ny,Nt), T_obs(Nt), T_d(Nt) )
call FW_Solve_2D(dt, tRange, IC, shape(IC), delta_x, delta_y, sol_exct, T_obs)

allocate( obs_eta(Nx, Ny, Nt), obs(size(x0_inds), size(T_obs)) )
! Exact surface wave at all x,y points for all t
obs_eta = sol_exct(1:Nx,:,:)
! Exact wave at measurement points x0,y0 for all t

DO k= 1, size(x0_inds)
    obs(k,:) = sol_exct(x0_inds(k),y0_inds(k),:)
end DO

! STEP 2: Distort exact IC to get "guess" IC eta0(:,1)
! Run forward and backwards solvers to get eta_a_t0_x
eta0(:,:,1) = 0.9*eta_exct0;

IC2(1:Nx,:) = eta0(:,:,1)
IC2(Nx+1:2*Nx, :) = u
IC2(2*Nx+1:3*Nx, :) = v

call FW_Solve_2D(dt, tRange, IC2, shape(IC2), delta_x, delta_y, sol_distorted, T_d)


allocate( eta_all(Nx,Ny,Nt), u_all(Nx,Ny,Nt), v_all(Nx,Ny,Nt), eta_x0(size(x0_inds), size(T_obs))  )

eta_all = sol_distorted(1:Nx,:,:)
u_all   = sol_distorted(Nx+1:2*Nx,:,:)
v_all   = sol_distorted(2*Nx+1:3*Nx,:,:)

DO j= 1, size(x0_inds)
    eta_x0(j,:) = sol_distorted(x0_inds(j),y0_inds(j),:)
end DO

IC_bw = 0.D0
allocate( Yb(3*Nx,Ny,Nt), Tb(Nt), eta_adj_t0(Nx,Ny), grad_J(Nx,Ny) )

call BW_Solve_2D(dt, tRange, IC_bw, shape(IC_bw), delta_x, delta_y, obs_eta, x0_inds, y0_inds, u_all, v_all, eta_all, &
 			Yb, Tb, eta_adj_t0, Ntotal, Ny_temp )


grad_J = -eta_adj_t0

!Define grad_J
!grad_J = -eta_adj_t0;
    
!Smooth gradient
!if smooth_grad
!    for k = 1:Ny
!        grad_J(:,k) = grad_smooth(grad_J(:,k),filt);
!    end
!    for l = 1:Nx
!        grad_J(l,:) = grad_smooth(grad_J(l,:),filt);
!    end
!end

 
 !eta0_prime_g = eta0(:,:,1);
 eta0_prime_g = amp*sin(Xmg)+sin(Ymg);

call logspace(10.D-13, 10.D-3, ep)

allocate( Gateaux_grad(size(ep)), Riesz_grad(size(ep)), kappa(size(ep)) )

DO l = 1,size(ep)
	! STEP 1(iv)
	! Compute eta_pert at the obs points given perturbed IC eta + ep*eta_prime
	eta_pert = eta0(:,:,1) + ep(l)*eta0_prime_g
	u_p = 0.D0
	v_p = 0.D0
	! Run forward solver to get eta_t at observation points given "initial guess" IC eta_0
	IC_pert(1:Nx,:) = eta_pert
	IC_pert(Nx+1:2*Nx, :) = u_p
	IC_pert(2*Nx+1:3*Nx, :) = v_p

	allocate(solpert(3*Nx,Ny,Nt), Tpert(Nt), etapert_x0_t_mult(size(x0_inds), size(T_obs)) )

	call FW_Solve_2D(dt, tRange, IC_pert, shape(IC_pert), delta_x, delta_y, solpert, Tpert)

	DO i2= 1, size(x0_inds)
    	etapert_x0_t_mult(i2,:) = solpert(x0_inds(i2),y0_inds(i2),:)
	end DO

	! STEP 2: Compute Gateaux derivative representation of gradient
	allocate(sumJ(Nt), sumP(Nt) )

	sumJ = 0.D0
	DO j2 = 1, size(x0_inds)
		sumJ = sumJ + (obs(j2,:) - eta_x0(j2,:))**2
	END DO
    J_eta = integral(Tpert, 0.5*sumJ);

    sumP = 0.D0

    DO l2 = 1, size(x0_inds)
    	sumP = sumP + (obs(L2,:) - etapert_x0_t_mult(L2,:))**2
    END DO

  
    Gateaux_grad(l) = ( integral(Tpert, 0.5*sumP) - J_eta)/ep(l)

    riesz_int = grad_J*eta0_prime_g

    IF (Ny>1) THEN
    	! trapz(x, riesz_int, 1) integrates every column so you're left with a row with Ny elements
    	! so you're integrating Nx elements with respect to x
    	! so for every column, i.e. DO i = 1,Ny
    	!									temprow(i) = integral(x, riesz_int(:,i))
    	!							END do
    	! Then integrate the remaining dimension(Ny) array with respect to y
    	! Riesz_grad(l) = integral(y, temprow)

    	allocate( temprow(Ny))
    	DO k2 = 1, Ny
    		temprow(k2) = integral(x, riesz_int(:,k2))
    	END DO
    	Riesz_grad(l) = integral(y, temprow)

    ELSE IF (Ny==1) THEN
    	DO k2 = 1, Ny
    		temprow2 = integral(x, riesz_int(:,i))
    	END DO
    	Riesz_grad(l) = temprow2
    ELSE
    	WRITE(*,*)  'You have entered an incorrect value for Ny'
    END IF

    kappa(l) = Gateaux_grad(l)/Riesz_grad(l)

DEALLOCATE(solpert, Tpert, etapert_x0_t_mult, sumJ, sumP,temprow )

END DO




print*, shape(kappa)
print*, kappa





! DON'T FORGET TO DEALLOCATE ALL THE DYNAMIC VARIABLES
   
   
end program kappa_test_2D