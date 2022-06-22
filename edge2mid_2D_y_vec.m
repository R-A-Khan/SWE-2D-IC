function [mid_Uy] = edge2mid_2D_y_vec(U)
% EDGE2MID_2D shifts all values of U currently found at midpoints of vertical cell edges of 
% a mxn grid (in the x-direction) to cell midpoints using mid_U(i,:) = (U(i+1,:) + U(i,:))/2 with
% periodic BC UN+1,:) = U(1,:).
%
% Input Arguments:
% U    = matrix size mxn; x = rows, y = cols
%
% Output Arguments:
% mid_Uy    = mxn matrix for cell midpoint values after y-shift
%
% Example Usage
% U = [ 1 2 3; 5 8 9; 4 1 1; 7 4 11];
% [mid_U] = edge2mid_2D_y(U);


mid_Uy = 0.5*(U+circshift(U,-1,2));
