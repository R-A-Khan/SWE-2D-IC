function [mid_Ux] = edge2mid_2D_x_vec(U)
% EDGE2MID_2D shifts all values of U currently found at midpoint of horizontal cell edges of 
% a mxn grid (in a y-direction) to cell midpoints using a forward difference formula mid_U(:,i) = (U(:,i+1) + U(:,i))/2 with
% periodic BC U(:,N+1) = U(:,1). 
%
% Input Arguments:
% U    = matrix size mxn; x = rows, y = cols
%
% Output Arguments:
% mid_Uy    = mxn matrix for cell midpoint values of U after shift in y
% direction
%
% Example Usage
% U = [ 1 2 3; 5 8 9; 4 1 1; 7 4 11];
% [mid_U] = edge2mid_2D_x_vec(U);

% [Nx,Ny] = size(U);
% mid_Ux  = zeros(Nx,Ny);
 

% For each col 1...Ny, shift values at x-edges to midpoints
mid_Ux= 0.5*(U+circshift(U,-1,1));

end
