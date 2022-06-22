function [dU] = back_diff_h_2D_vec( U, delta, dim)
% BACK_DIFF_H computes the BACKWARD DIFFERENCE approximation for dU0 w.r.t.
% x or y with periodic BC U(:,N+1) = U(:,1), U(N+1,:) = U(1,:).
%
% Input Arguments:
% N     = size of dimension dim for U; x = m, y = n
% U0    = matrix size mxn; x = rows, y = cols
% delta = size of step in x or y
% dim   = dimension U0 differentiated w.r.t.; x=1,y=2
%
% Output Arguments:
% dU    = mxn matrix approximation for dU/dx or dU/dy
%
% Example Usage
%
% delta_x = pi/4; delta_y = pi/8;
% x = 0:delta_x:2*pi;
% y = 0: delta_y:2*pi;
% [X,Y] = meshgrid(x,y);
% f = sin(X)+sin(Y);
% dfdx_exct = cos(X);
% [dfdx_approx] = forw_diff_u_2D(U, delta_x, 1)
% err = abs( dfdx_exct - dfdx_approx );
% figure; surf(X,Y,err);

if dim == 1
    dU = ( U - circshift(U,1,1) )/delta;   
elseif dim == 2 
    dU = ( U - circshift(U,1,2) )/delta;
end
        
    





end

