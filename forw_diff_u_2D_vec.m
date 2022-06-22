function [dU] = forw_diff_u_2D_vec( U, delta, dim)
% FORW_DIFF_U computes the FORWARD DIFFERENCE approximation for dU0 w.r.t.
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


if dim == 1

dU = ( circshift(U,-1,1) - U )/delta;

elseif dim == 2

dU = ( circshift(U,-1,2) -  U )/delta;

end



end