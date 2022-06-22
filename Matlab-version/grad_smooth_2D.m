function [smoothed]=grad_smooth_2D(signal,filt)
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

% 1. take fftshift of the signal, now the zero freq. component is at centre
fJ = fftshift(fJ);
% figure(1);
% subplot(1,2,1); imagesc(abs(fJ))


% 2. define fs_x and fs_y to be from -N/2 :N/2 -1
fs_x= -Nx/2: Nx/2 - 1;
fs_y= -Ny/2: Ny/2 - 1;
% 3. Apply filter to signal
 
fmin_x = round(filt*max(fs_x)); % Minimum filtering frequency (corresponding to 1/L in the notes)
fmin_y = round(filt*max(fs_y)); 

% fs_x = fftshift(0:Nx-1); % Frequencies
% fs_y = fftshift(0:Ny-1); % Frequencies
fJ_H = fJ;
for k = 1:Nx
    for j =1:Ny
        fJ_H(k,j) = fJ(k,j)/(1+ (fs_x(k)^2/fmin_x^2) + (fs_y(j)^2/fmin_y^2) );
    end
end


% figure(1); subplot(1,2,2);
% imagesc(abs((fJ_H)))

% 4. Apply ifftshift to smoothed signal
fJ_H = ifftshift(fJ_H);
% 5. Apply ifft to smoothed signal



% % Inverse Fourier transform
smoothed=ifft2(fJ_H);