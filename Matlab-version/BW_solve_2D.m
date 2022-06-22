function [b, T_adj, Y] = BW_solve_2D(h, hy, dt, u0, tmin, tmax, x0_inds, y0_inds, trend, obs, u_all, v_all, eta_all)
% Usage: 
% [b, T, Y] = BW_solve(h, dt, u0, tmin, tmax, x0_inds, trend, obs, eta, u_all, eta_all, T_obs, T_eta, nl)
%
% Solves backwards (adjoint) equations using observations results from
% forward solver
%
% Input:
% h       = grid spacing
% dt      = time step
% u0      = initial conditions for height and velocity
% tmin    = minimum time
% tmax    = maximum time (control time)
% x0_inds = indices of observation positions
% trend   = function handle for trend (rhs of ode in time)
% obs     = vector of height at observation position from t=tmin to t=tmax
% eta     = height at control time for all positions
% u_all   = velocity at all positions and times
% eta_all = height at all positions and times
% T_obs   = times for all observations
% T_eta   = times for all height data from forward solver
% nl      = true = nonlinear equations, false = linear equations
%
% Output:
% b = height at time tmin for all (x,y) positions
% T = numSteps x 1 vector of discrete time steps between tmin and tmax
% Y = 3N x Ny x numSteps matrix
%     1:N rows = height, N+1 : 2N rows = x velocity, 2N+1 : 3N rows = y velocity
%     Rows   = values of height/velocity at point x_i for all times
%     Colums = values of height/velocity at time t_i for all positions

% Solving backwards with tau = tmax - t
[Nx, Ny, T] = size(u_all);

obs_tau     = zeros(Nx,Ny,T);
u_all_tau   = zeros(Nx,Ny,T);
v_all_tau   = zeros(Nx,Ny,T);
eta_all_tau = zeros(Nx,Ny,T);

tau = T;
for i = 1:T
    obs_tau(:,:,i)   = obs(:,:,tau);
    u_all_tau(:,:,i)   = u_all(:,:,tau);
    v_all_tau(:,:,i)   = v_all(:,:,tau);
    eta_all_tau(:,:,i) = eta_all(:,:,tau);
    tau = tau - 1;
end

% surf( abs(obs_tau(:,:,1)-obs(:,:,end)) )

tRange = [tmin tmax];
[T_adj,Y] = RK34_BW_2D(dt, u0, tRange, trend, h, hy, x0_inds, y0_inds, obs_tau, u_all_tau, v_all_tau, eta_all_tau);
% figure(1);
% for i = 1:round(length(T)/1)
%     Z = Y(1:Nx,:,i);
%     surf(Z, 'EdgeColor','none');
%     zlim([-0.5 0.5])
% %     plot(x,sol(1:Nx,round(Ny/2),i),'linewidth',2)
% %     ylim([-0.05 0.15])
%     xlabel('x'); 
%     ylabel('y');
%     M(i)=getframe(gcf);
% end

b = Y(1:Nx,:,end);
