module ParametersLineMinimization

implicit none
  
  ! Input variables 
  integer, intent(in) :: Nx, Ny, Nt, N_obs
  double precision, intent(in) :: dt, delta_x, delta_y
  integer, dimension(2) :: tRange
  double precision,  dimension(:, :), intent(in) :: grad_J, eta0
  double precision,  dimension(:,:), intent(in) :: obs
  double precision,  dimension(:) , intent(in) :: T
  integer,  dimension(:), intent(in) :: x0_inds, y0_inds

  ! Output variables:
  double precision, intent(out) :: cost




  ! Local Variables
  double precision,  dimension(Nt)  :: cost_opt, Tsol, Tb
  double precision,  dimension(N_obs, Nt)    :: eta_new_x0
  double precision,  dimension(Nx, Ny) :: eta_new, eta_adj_t0, IC_bw, grad_J_new, grad_J_tau
  double precision,  dimension(Nx*3, Ny, Nt) :: Y, Yb
  double precision,  dimension(Nx, Ny, Nt) :: eta_all, u_all, v_all
  integer :: j,j2, Ntotal, Ny_temp


! INPUT AND OUTPUT VARIABLES
 double precision :: tau


end module
