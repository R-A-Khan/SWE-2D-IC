function trend = BW_SWE_NL_2D_vec_bath(t, y, delta_x, delta_y, x0_inds, y0_inds, obs, u_all, v_all, eta_all, beta)
% Usage: trend = bw_SWE_periodic_mult_NL_v2(t, y, h, x0_inds, obs, eta, u_all, eta_all)
%
% Finds trend of backwards (adjoint) nonlinear SWE
% using second order finite difference / finite volume approximation
% grid size is 2*N, with dh/dt evaluated at centres and du/dt at edges
%
% Input:
% t       = time (not used)
% y       = input variables (height and velocity in one vector)
% N       = number of grid points in space
% h       = spacing of grid points
% x0_inds = indices of observation positions
% obs     = height at observation position from t=tmin to t=tmax
% eta     = height from forward solver using approximate initial conditions
% u_all   = velocity at all positions and times
% eta_all = height at all positions and times
%
% Output:
% trend = trend for both height and velocity equations in one vector

[N, Ny] = size(y);
Nx = N/3;
eta_adj = y(1:Nx,:);
u_adj   = y(Nx+1:2*Nx,:);
v_adj   = y(2*Nx+1:3*Nx,:);
u       = u_all;
v       = v_all;
eta     = eta_all;

trend = zeros(3*Nx, Ny);
% eta_size = size(trend(1:N,:) )
% eta_adj_size = size(eta_adj)
% u_size = size(u)
% v_size = size(v)
% u_adj_size = size(u_adj)
% v_adj_size = size(v_adj)

trend(1:Nx,:)        = - ( edge2mid_2D_x_vec(u) .* cent_diff_2D_vec(eta_adj, delta_x, 1) ...
                         + edge2mid_2D_y_vec(v) .* cent_diff_2D_vec(eta_adj, delta_y, 2) ...
                         + forw_diff_u_2D_vec(u_adj, delta_x,1) + forw_diff_u_2D_vec(v_adj, delta_y, 2) );
                     
trend(Nx+1:2*Nx,:)   = - ( u .* cent_diff_2D_vec(u_adj, delta_x, 1) ...
                         + mid2edge_2D_x_vec(edge2mid_2D_y_vec(v)) .* cent_diff_2D_vec(u_adj, delta_y, 2) ...
                         + (1 + mid2edge_2D_x_vec(eta-beta) ) .* back_diff_h_2D_vec(eta_adj, delta_x, 1) );
trend(2*Nx+1:3*Nx,:) = - ( v .* cent_diff_2D_vec(v_adj, delta_y, 2) ...
                         + edge2mid_2D_x_vec(mid2edge_2D_y_vec(u)) .* cent_diff_2D_vec(v_adj, delta_y, 2) ...
                         + (1 + mid2edge_2D_y_vec(eta - beta) ) .* back_diff_h_2D_vec(eta_adj, delta_y, 2) );
% Add observations
for i = 1:length(x0_inds)
        trend(x0_inds(i), y0_inds(i)) = trend(x0_inds(i), y0_inds(i))...
        + (obs(x0_inds(i), y0_inds(i)) - eta(x0_inds(i), y0_inds(i)))/(delta_x*delta_y);
end


% pos = any(trend(:)>0);
% neg = any(trend(:)<0);
% disp([sprintf('%d',pos), '   ', sprintf('%d', neg)]);