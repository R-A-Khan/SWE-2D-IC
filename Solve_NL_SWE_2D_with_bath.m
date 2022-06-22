function [Tsol, sol] = Solve_NL_SWE_2D_with_bath( dt, tRange, trend, u0, delta_x, delta_y, nu, beta )
% SOLVE_ADVEC solves the 2D advection equation with velocity
% field v(x,y)=(-y,x) for some transported scalar function u0(x,y,t) 
%   
% Input Arguments
%
% dt = step size for t
% u0 = initial value h(x,y,0)
% tRange = time interval
% trend = approximation function for advection term
% u_x = x-component of velocity field
% u_y = y_component of velocity field
% delta_x = step size for x
% delta_y = step size for y
%
% Output Arguments
%
% h = m x n x T matrix representing solution h(x,y,t)
% Tsol = vector of time steps 


[Tsol, sol] = RK34_FW_SWE_2D_FD_with_bath(dt, tRange, trend, u0, delta_x, delta_y, nu, beta);
end

