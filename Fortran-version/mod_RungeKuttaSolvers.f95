Module RungeKuttaSolvers
USE SWEandAdjointSolvers

contains

 
subroutine RK34_FW(dt, tRange, Y0, shapeY0, delta_x, delta_y, sol, tSol)
implicit none

!real, external :: func
double precision, intent(in) :: dt, delta_x, delta_y
double precision, dimension(2), intent(in) :: tRange
integer, dimension(2), intent(in) :: shapeY0
double precision, dimension(:,:), intent(in) :: Y0

double precision, allocatable, intent(out) :: sol(:,:,:)
double precision, allocatable, dimension(:), intent(out) :: tSol



! local variables
integer :: k, numSteps, Nx, Ny
integer, dimension(2):: sizeY0
double precision, dimension(4,4):: alpha, beta
double precision, allocatable, dimension(:,:) :: U0, U1, U2, U3, U4 

alpha = reshape( (/1.D0,   0.D0,    0.D0,  0.D0, &
                   0.D0,   1.D0,    0.D0,  0.D0, &
        (2.D0)/(3.D0),   0.D0, 1.D0/3.D0,  0.D0, &
                   0.D0,   0.D0,    0.D0,  1.D0 /), (/4,4/) )
                   

                   
beta = reshape( (/ 1./2.D0,       0.D0,      0.D0,     0.D0, &
                      0.D0,  1.D0/2.D0,      0.D0,     0.D0, &
                      0.D0,       0.D0, 1.D0/6.D0,     0.D0, &
                      0.D0,       0.D0,    0.D0,  1.D0/2.D0 /), (/4,4/) )     
                    
      

numSteps = floor ( abs( (tRange(2) - tRange(1))/dt ))

Nx = shapeY0(1)
Ny = shapeY0(2)

allocate(tSol(numSteps+1), sol(Nx,Ny, numSteps+1), U0(Nx,Ny), &
U1(Nx,Ny), U2(Nx,Ny), U3(Nx,Ny), U4(Nx,Ny))

sol(1:Nx,1:Ny,1) = Y0
tSol(1) = tRange(1)

do k = 1, numSteps
    tSol(k+1) = tSol(k) + dt

    U0 = sol(1:Nx,1:Ny,k)

    U1 = alpha(1,1)*U0                 + dt*beta(1,1)*solve_swe(U0, shape(U0), delta_x, delta_y)   

    
    U2 = alpha(2,2)*U1                 + dt*beta(2,2)*solve_swe(U1, shape(U1), delta_x, delta_y)

    
    U3 = (2.D0/3.D0)*U0 + alpha(3,3)*U2    + dt*beta(3,3)*solve_swe(U2, shape(U2), delta_x, delta_y)

    U4 = alpha(4,4)*U3                 + dt*beta(4,4)*solve_swe(U3, shape(U3), delta_x, delta_y)
    
		sol(1:Nx, 1:Ny, k+1) = U4

end do

 !deallocate(U0,U1,U2, U3, U4)
end subroutine RK34_FW

!______________________________________________________ 

subroutine RK34_BW(dt, tRange, Y0, shapeY0, delta_x, delta_y ,&
obs, x0_inds, y0_inds, u_all, v_all, eta_all, sol, tSol, Ntotal, Ny)
implicit none

!real,external :: func
double precision, intent(in) :: dt, delta_x, delta_y
double precision, dimension(2), intent(in) :: tRange
integer, dimension(2), intent(in) :: shapeY0
double precision, dimension(:,:), intent(in) :: Y0
integer, dimension(:), intent(in) :: x0_inds, y0_inds
double precision,  dimension(:,:,: ), intent(in) :: u_all, v_all, eta_all, obs


double precision, allocatable, dimension(:,:,:), intent(out) :: sol
double precision, allocatable, dimension(:), intent(out) :: tSol
integer, intent(out) ::  Ntotal, Ny


! local variables
integer :: k, numSteps, Nx
integer, dimension(2):: sizeY0
integer, dimension(3):: shape_u_all 
double precision, dimension(4,4):: alpha, beta
double precision, allocatable, dimension(:,:) :: U0, U1, U2, U3, U4 

shape_u_all = shape(u_all)
Nx = shape_u_all(1)


alpha = reshape( (/1.D0,   0.D0,    0.D0,  0.D0, &
                   0.D0,   1.D0,    0.D0,  0.D0, &
        (2.D0)/(3.D0),   0.D0, 1.D0/3.D0,  0.D0, &
                   0.D0,   0.D0,    0.D0,  1.D0 /), (/4,4/) )
                   

                   
beta = reshape( (/ 1.D0/2.D0,       0.D0,      0.D0,     0.D0, &
                      0.D0,  1.D0/2.D0,      0.D0,     0.D0, &
                      0.D0,       0.D0, 1.D0/6.D0,     0.D0, &
                      0.D0,       0.D0,    0.D0,  1.D0/2.D0 /), (/4,4/) )     
                       
                                    

numSteps = floor ( abs( (tRange(2) - tRange(1))/dt ))


Ntotal = shapeY0(1)
Ny = shapeY0(2)

allocate( tSol(numSteps+1), sol(Ntotal,Ny, numSteps+1), &
U0(Ntotal,Ny), U1(Ntotal,Ny), U2(Ntotal,Ny), U3(Ntotal,Ny), U4(Ntotal,Ny) )


sol(1:Ntotal,1:Ny,1) = Y0
tSol(1) = tRange(1)


do k = 1, numSteps
    tSol(k+1) = tSol(k) + dt
    U0 = sol(1:Ntotal,1:Ny,k)
    
    U1 = 1.D0*U0  & 
    + dt*(1.D0/2.D0)*solve_swe_adj(U0, shape(U0), delta_x, delta_y, x0_inds, y0_inds, & 
    obs(1:Nx,1:Ny,k), u_all(1:Nx,1:Ny,k), v_all(1:Nx,1:Ny,k), eta_all(1:Nx,1:Ny,k))
    !print*, "U1 at iteration", k
    !print*, U1
    
    U2 = 1.D0*U1  & 
    + dt*beta(2,2)*solve_swe_adj(U1, shape(U1), delta_x, delta_y, x0_inds, y0_inds, &
    obs(1:Nx,1:Ny,k), u_all(1:Nx,1:Ny,k), v_all(1:Nx,1:Ny,k), eta_all(1:Nx,1:Ny,k))
    !print*, "U2 at iteration", k
    !print*, U2
       
    
    U3 = (2.D0/3.D0)*U0 + alpha(3,3)*U2 &
    + dt*beta(3,3)*solve_swe_adj(U2, shape(U2), delta_x, delta_y, x0_inds, y0_inds, &
    obs(1:Nx,1:Ny,k), u_all(1:Nx,1:Ny,k), v_all(1:Nx,1:Ny,k), eta_all(1:Nx,1:Ny,k))
    !print*, "U3 at iteration", k
    !print*, U3

    
    U4 = alpha(4,4)*U3   &
    + dt*beta(4,4)*solve_swe_adj(U3, shape(U3), delta_x, delta_y, x0_inds, y0_inds, &
    obs(1:Nx,1:Ny,k), u_all(1:Nx,1:Ny,k), v_all(1:Nx,1:Ny,k), eta_all(1:Nx,1:Ny,k))
		sol(1:Ntotal, 1:Ny, k+1) = U4
		!print*, "solution at iteration", k+1
		!print*, sol(1:Ntotal, 1:Ny, k+1)

end do

 !deallocate(U0,U1,U2, U3, U4)

end subroutine RK34_BW

!______________________________________________________ 

end module
