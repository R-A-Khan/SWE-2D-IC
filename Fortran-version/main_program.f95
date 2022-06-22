program test
    use FWandBWSolvers


    implicit none

     double precision, dimension(6,3) :: a, fw_swe_a, bw_swe_adj_a, IC
     double precision, dimension(2,3) :: obs, eta_all, u_all, v_all, temp
     double precision, dimension(2,3,6) :: obs_full, eta_full, u_full, v_full
     double precision, allocatable, dimension(:,:,:) :: rk_fw_out, rk_bw_out, fw_out, bw_out
     double precision, allocatable, dimension(:,:) :: adjoint_t0
     double precision, allocatable, dimension(:) :: T, Tb
     integer :: M, N, Mb, Nb, j
     double precision:: stepsize_x, stepsize_y, stepsize_t
     integer, dimension(2):: x0_inds, y0_inds
     double precision, dimension(2) :: tRange

     
    IC = 1
    temp = 1
    a = reshape((/ 1.D0,3.D0,3.D0,5.D0,4.D0,7.D0,2.D0,6.D0,5.D0,8.D0,1.D0,4.D0,3.D0,9.D0,2.D0,9.D0,1.D0,11.D0/),(/6,3/))
    eta_all = a(1:2,1:3)
    u_all = a(3:4,1:3)
    v_all = a(5:6,1:3)
    stepsize_x = 3.D0
    stepsize_y = 6.D0
    stepsize_t = 1.D-1
    obs = reshape((/ 2.D0,3.D0,1.D0,3.D0,5.D0,7.D0/), (/2,3/))
    x0_inds = (/1,2/)
    y0_inds = (/1,3/)
    tRange = (/0.D0,5.D-1/)
    
    do j = 1, 6
    eta_full(1:2,1:3,j) = a(1:2,1:3)+ 1.D-1*j/2*temp
    u_full(1:2,1:3,j) = a(3:4,1:3)+ 1.D-1*j/2*temp
    v_full(1:2,1:3,j) = a(5:6,1:3)+ 1.D-1*j/2*temp
    obs_full(1:2,1:3,j) = obs + 1.D-1*j/2*temp
    end do

   ! shape_a = shape(a)
   bw_swe_adj_a =  &
   solve_swe_adj(IC, shape(IC), stepsize_x, stepsize_y, x0_inds, y0_inds, obs, u_all, v_all, eta_all)
   
   !fw_swe_a = solve_SWE(a, shape(a), stepsize_x, stepsize_y)
   
  !call RK34_FW(stepsize_t, tRange,  a, shape(a), stepsize_x, stepsize_y, rk_fw_out, T, M, N)
   
  !call RK34_BW(stepsize_t, tRange,  IC, shape(IC), stepsize_x, stepsize_y, &
        !obs_full, x0_inds, y0_inds, u_full, v_full, eta_full, rk_bw_out, Tb, Mb, Nb)
  
  call FW_Solve_2D(stepsize_t, tRange, a, shape(a), stepsize_x, stepsize_y, fw_out, T)
        
  !call BW_Solve_2D(stepsize_t, tRange,  IC, shape(IC), stepsize_x, stepsize_y, &
             ! obs_full, x0_inds, y0_inds, u_full, v_full, eta_full, bw_out, Tb,  adjoint_t0, Mb, Nb)
        

    
   !print*, "this is the output from function solve_swe"
   !print*,  fw_swe_a
   !print*, "this is the output from function solve_swe_adj"
   !print*,  bw_swe_adj_a
   print*, "this is the output of function bw_solve at the last time step"
   print*,  fw_out(:,:,6)
   print*, shape(fw_out)
   !print*, "this is the time vector T"
   !print*, T

    
end program