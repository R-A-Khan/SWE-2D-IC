module FWandBWSolvers
USE RungeKuttaSolvers

contains

subroutine FW_Solve_2D(dt, tRange, Y0, shapeY0, delta_x, delta_y, sol, tSol)
implicit none 
double precision, intent(in) :: dt, delta_x, delta_y
double precision, dimension(2), intent(in) :: tRange
integer, dimension(2), intent(in) :: shapeY0
double precision, dimension(:,:), intent(in) :: Y0
double precision, allocatable,  dimension(:,:,:), intent(out) :: sol
double precision, allocatable, dimension(:), intent(out) :: tSol


call RK34_FW(dt, tRange, Y0, shapeY0, delta_x, delta_y, sol, tSol)



end subroutine FW_Solve_2D

!______________________________________________________ 


subroutine BW_Solve_2D(dt, tRange, Y0, shapeY0, delta_x, delta_y, &
										   obs_all, x0_inds, y0_inds, u_all, v_all, eta_all, sol, tSol, b, Ntotal, Ny)
implicit none

double precision, intent(in) :: dt, delta_x, delta_y
double precision, dimension(2), intent(in) :: tRange
integer, dimension(2), intent(in) :: shapeY0
double precision, dimension(:,:), intent(in) :: Y0
integer, dimension(:), intent(in) :: x0_inds, y0_inds
double precision,  dimension(:,:,: ), intent(in) :: u_all, v_all, eta_all, obs_all


double precision,  allocatable, dimension(:,:,:), intent(out) :: sol
double precision,  allocatable, dimension(:), intent(out) :: tSol
double precision,  allocatable, dimension(:,:), intent(out) :: b
integer, intent(out) :: Ny, Ntotal

!Local variable declaration

double precision, allocatable, dimension(:,:,:) :: obs_tau, u_all_tau, v_all_tau, eta_all_tau
integer, dimension(3) :: shape_u_all
integer::T, tau, i, Nx

shape_u_all = shape(u_all)
Nx = shape_u_all(1)
Ntotal = shapeY0(1)
Ny = shape_u_all(2)
T  = shape_u_all(3)

allocate(obs_tau(Nx,Ny,T), u_all_tau(Nx,Ny,T), v_all_tau(Nx,Ny,T), eta_all_tau(Nx,Ny,T), b(Nx,Ny) )

tau = T
do i = 1, T
    obs_tau(1:Nx, 1:Ny, i) =  obs_all(1:Nx, 1:Ny, tau); 
    u_all_tau(1:Nx, 1:Ny, i) = u_all(1:Nx, 1:Ny, tau); 
    v_all_tau(1:Nx, 1:Ny, i) = v_all(1:Nx, 1:Ny, tau); 
    eta_all_tau(1:Nx, 1:Ny, i) = eta_all(1:Nx, 1:Ny, tau); 
    tau = tau -1
end do


call RK34_BW(dt, tRange, Y0, shapeY0, delta_x, delta_y , &
						 obs_tau, x0_inds, y0_inds, u_all_tau, v_all_tau, eta_all_tau, sol, tSol, Ntotal, Ny)
b = sol(1:Nx,1:Ny, size(tSol))

!deallocate(obs_tau, u_all_tau, v_all_tau, eta_all_tau)

end subroutine BW_Solve_2D
!______________________________________________________ 

end module


