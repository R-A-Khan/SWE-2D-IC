module SWEandAdjointSolvers

USE procedures

contains

function solve_swe(U0, shape_U0, delta_x, delta_y)
integer, dimension(2) :: shape_U0
double precision, dimension(shape_U0(1),shape_U0(2)), intent(in) :: U0
double precision, intent(in) :: delta_x, delta_y
double precision :: solve_swe(shape_U0(1),shape_U0(2))


!local variables only within scope of function
integer :: Nx, Ny
double precision, dimension((shape_U0(1))/3,shape_U0(2)) :: eta, u, v

Nx = (shape_U0(1))/3
Ny = shape_U0(2)


eta = U0(1:Nx,1:Ny)
u = U0(Nx+1:2*Nx, 1:Ny)
v = U0(2*Nx+1:3*Nx, 1:Ny)

solve_swe(1:Nx,1:Ny) = - (  forw_diff( (m2e_x( 1.D0 +eta,(/Nx,Ny/) )*u) , (/Nx,Ny/), delta_x, 1) &
                          + forw_diff( (m2e_y( 1.D0 +eta,(/Nx,Ny/) )*v) , (/Nx,Ny/), delta_y, 2)  )

solve_swe(Nx+1:2*Nx,1:Ny) = - back_diff(eta,(/Nx,Ny/), delta_x, 1 )  &
                            - back_diff(  5.D-1*( e2m_x(u**2, (/Nx,Ny/)) + e2m_y(v**2,(/Nx,Ny/) ) )  , (/Nx,Ny/), delta_x, 1  ) & 
                            + (1.D0/m2e_x( 1.D0+eta,(/Nx,Ny/) ) )*( back_diff( e2m_y(v,(/Nx,Ny/) ), (/Nx,Ny/), delta_x, 1 ) & 
                            - cent_diff(u, (/Nx,Ny/), delta_y,2) ) &
                            * ( e2m_y( ( m2e_x( ( m2e_y(1.D0+eta,(/Nx,Ny/))*v), (/Nx,Ny/) ) ),  (/Nx,Ny/) ))

solve_swe(2*Nx+1:3*Nx,1:Ny) = - back_diff(eta,(/Nx,Ny/), delta_y, 2 )  &
                              - back_diff(  5.D-1*( e2m_x(u**2, (/Nx,Ny/)) + e2m_y(v**2,(/Nx,Ny/) ) ) , (/Nx,Ny/), delta_y, 2  ) &
                              - (1.D0/m2e_y( 1.D0+eta,(/Nx,Ny/) ) )*( cent_diff(v, (/Nx,Ny/), delta_x,1) & 
                              - back_diff( e2m_x(u,(/Nx,Ny/) ), (/Nx,Ny/), delta_y, 2 ) ) &
                              * ( e2m_x( (m2e_y( (m2e_x(1.D0+eta,(/Nx,Ny/) )*u), (/Nx,Ny/) ) ), (/Nx,Ny/) ) )



end function solve_swe

 !______________________________________________________  

function solve_swe_adj(U0, shape_U0, delta_x, delta_y, x0_inds, y0_inds, obs, u_all, v_all, eta_all)

! Variable declaration
integer, dimension(2) :: shape_U0
integer :: Nx, Ny, i
integer, dimension(:), intent(in) :: x0_inds, y0_inds
double precision, dimension(shape_U0(1),shape_U0(2)), intent(in) :: U0
double precision, dimension(shape_U0(1)/3,shape_U0(2)), intent(in) :: u_all, v_all, eta_all
double precision, dimension(:,:), intent(in) :: obs
double precision, intent(in) :: delta_x, delta_y
double precision :: solve_swe_adj(shape_U0(1),shape_U0(2))
double precision, dimension(shape_U0(1)/3,shape_U0(2)):: eta_adj, u_adj, v_adj


Nx = shape_U0(1)/3
Ny = shape_U0(2)

eta_adj = U0(1:Nx, 1:Ny)
u_adj = U0(Nx+1:2*Nx, 1:Ny)
v_adj = U0(2*Nx+1:3*Nx, 1:Ny)


solve_swe_adj(1:Nx,1:Ny) =  - ( e2m_x(u_all, (/Nx,Ny/)) * cent_diff(eta_adj, (/Nx,Ny/), delta_x, 1) &
                              + e2m_y(v_all, (/Nx,Ny/)) * cent_diff(eta_adj, (/Nx,Ny/), delta_y, 2) &
                              + forw_diff(u_adj, (/Nx,Ny/), delta_x, 1) + forw_diff(v_adj, (/Nx,Ny/), delta_y, 2) )

solve_swe_adj(Nx+1:2*Nx,1:Ny) = - ( u_all * cent_diff(u_adj, (/Nx,Ny/), delta_x, 1) &
                                  + m2e_x( e2m_y( v_all, (/Nx,Ny/) ), (/Nx,Ny/) ) * cent_diff(u_adj, (/Nx,Ny/), delta_y, 2) &
                                  + (1.D0 + m2e_x(eta_all, (/Nx,Ny/) ) ) * back_diff(eta_adj, (/Nx,Ny/), delta_x, 1) )


solve_swe_adj(2*Nx+1:3*Nx,1:Ny) = - ( v_all * cent_diff(v_adj, (/Nx,Ny/), delta_y, 2) &
                                  + e2m_x( m2e_y( u_all, (/Nx,Ny/) ), (/Nx,Ny/) ) * cent_diff(v_adj, (/Nx,Ny/), delta_y, 2) &
                                  + (1.D0 + m2e_y(eta_all, (/Nx,Ny/) ) ) * back_diff(eta_adj, (/Nx,Ny/), delta_y, 2) )



do i = 1, size(x0_inds)
solve_swe_adj(x0_inds(i), y0_inds(i)) = solve_swe_adj(x0_inds(i), y0_inds(i)) &
                     + ( obs(x0_inds(i), y0_inds(i)) - eta_all(x0_inds(i), y0_inds(i)) )/(delta_x*delta_y) 
 end do

end function solve_swe_adj

end module

!----------------------------------------------------------

