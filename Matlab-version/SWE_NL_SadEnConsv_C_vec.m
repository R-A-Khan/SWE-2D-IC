function [trend] = SWE_NL_SadEnConsv_C_vec(t, u0, delta_x, delta_y, nu)
% SWE_NL_2D_FDapp roximates solution for 2D SWE equations using:
% Sadourny energy conserving scheme on a C-grid
% 
% Input Arguments:
% u0 = verticall concatenated h,u,v
% h = height found at cell centres each time step 
% u = velocity in x-direction found at vertical edge midpoints at each time
% step
% v = velocity in y-direction found at horizontal edge midpoints at each time
% step
% nu = kinematic viscosity co-efficient 
%
% Output Arguments:
% dh/dt    = trend(1:Nx,:)
% du/dt = trend(Nx+1:2*N,:) 
% dv/dt = trend(2*Nx+1:3*Nx,:)

[m,n] = size(u0);
Nx = m/3; Ny = n;
trend = zeros(3*Nx,Ny);
eta = u0(1:Nx,:);
u = u0(Nx+1:2*Nx,:);
v = u0(2*Nx+1:3*Nx,:);

% 2D NL SWE equation

trend(1:Nx,:)        = -(  forw_diff_u_2D_vec( mid2edge_2D_x_vec(1+eta).*u, delta_x, 1) + forw_diff_u_2D_vec( mid2edge_2D_y_vec(1+eta).*v, delta_y, 2)  );

trend(Nx+1:2*Nx,:)   = ...
  - back_diff_h_2D_vec(eta, delta_x, 1) ...
  - back_diff_h_2D_vec( 0.5*(edge2mid_2D_x_vec(u.^2)+edge2mid_2D_y_vec(v.^2)), delta_x, 1) ...
  + (1./mid2edge_2D_x_vec(1+eta)).*(back_diff_h_2D_vec(edge2mid_2D_y_vec(v), delta_x, 1) - cent_diff_2D_vec(u, delta_y,2)).*...
  (edge2mid_2D_y_vec(mid2edge_2D_x_vec(mid2edge_2D_y_vec(1+eta).*v)))...
  + nu*( cent_diff_2D_vec( cent_diff_2D_vec(u, delta_x, 1), delta_x, 1) + cent_diff_2D_vec( cent_diff_2D_vec(u, delta_y, 2), delta_y, 2) );
                    
trend(2*Nx+1:3*Nx,:) = ...
  - back_diff_h_2D_vec(eta, delta_y, 2) ...
  - back_diff_h_2D_vec( 0.5*(edge2mid_2D_x_vec(u.^2)+edge2mid_2D_y_vec(v.^2)), delta_y, 2) ...
  - (1./mid2edge_2D_y_vec(1+eta)).*(cent_diff_2D_vec(v, delta_x, 1) - back_diff_h_2D_vec(edge2mid_2D_x_vec(u), delta_y, 2)).*...
  (edge2mid_2D_x_vec(mid2edge_2D_y_vec(mid2edge_2D_x_vec(1+eta).*u))) ...
  + nu*( cent_diff_2D_vec( cent_diff_2D_vec(v, delta_x, 1), delta_x, 1) + cent_diff_2D_vec( cent_diff_2D_vec(v, delta_y, 2), delta_y, 2) );

end