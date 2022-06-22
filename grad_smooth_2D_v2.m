function [smoothed]=grad_smooth_2D_v2(signal,filt)
% Usage:
% smoothed = grad_smooth(signal,filt)
%
% Input:
% signal = signal to be smoothed
% fmin   = set to (1+L^2k^2) all Fourier modes greater than filt*max frequency
% (0<filt<1)
% Output:
% smoothed = smoothed sobolev gradient version of signal
 
[Nx, Ny] = size(signal);
 
% Take fft of signal
fJ=fft2(signal);
 
fs_x = 0:Nx/2; % Frequencies
fs_y = 0:Ny/2; % Frequencies
 
fmin_x = round(filt*max(fs_x)); % Minimum filtering frequency (corresponding to 1/L in the notes)
fmin_y = round(filt*max(fs_y)); 
fJ_H = fJ;
for k = 1:Nx/2+1
    for j =1:Ny/2+1
        fJ_H(k,j) = fJ(k,j)/(1+fs_x(k)^2/fmin_x^2 + fs_y(j)^2/fmin_y^2 );
    end
end
 
% Use symmetry property of fft of real signals to deal with remaining Fourier modes 
% Frequencies associated with the elements 1:N/2+1 of the fft vector are 0:N/2
for k=Nx/2+2:Nx
    for j = Ny/2+2:Ny
        fJ_H(k,j)=conj(fJ_H(Nx-k+2, Ny-j+2));
    end
end
 
% Inverse Fourier transform
smoothed=ifft2(fJ_H);