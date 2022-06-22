
fileID = fopen( 'ntrial_file.dat');
a = fread(fileID, inf, 'uint');
ntrial = a;

fileID = fopen( 'res_xyt.dat');
b = fread(fileID, 3, 'uint');
Nx =b(1); Ny = b(2); Nt = b(3);

fileID = fopen( 'n_obs_x_file.dat');
c = fread(fileID, inf, 'uint');
n_obs_x = c;

fileID = fopen( 'n_obs_y_file.dat'); 
d = fread(fileID, inf, 'uint');
n_obs_y = d;

fileID = fopen( 'x0_min_file.dat'); 
f = fread(fileID, inf, 'double');
x0_min = f;

fileID = fopen( 'y0_min_file.dat');
g = fread(fileID, inf, 'double');
y0_min = g;

fileID = fopen( 'dmu_x_file.dat'); 
h = fread(fileID, inf, 'double');
dmu_x = h;

fileID = fopen( 'dmu_y_file.dat');
l = fread(fileID, inf, 'double');
dmu_y = l;

fileID = fopen( 'iter_max_file.dat'); 
o = fread(fileID, inf, 'uint');
iter_max = o;

fileID = fopen( 'obs_vals_file.dat');
p = fread(fileID, inf, 'double');
n = length(p)/2;
obs_vals = reshape(p, [n 2]);

fileID = fopen( 'eta_optimum_file.dat'); 
q = fread(fileID, inf, 'double');
eta_optimum = reshape(q, [Nx Ny]);

fileID = fopen( 'eta_exct0_file.dat'); 
r = fread(fileID, inf, 'double');
eta_exct0 = reshape(r, [Nx Ny]);

fileID = fopen( 'eta0_file.dat'); 
s = fread(fileID, inf, 'double');
eta0 = reshape(s, [Nx Ny iter_max]);

fileID = fopen( 'Xmg_file.dat');
u = fread(fileID, inf, 'double');
X  = reshape(u, [Nx Ny]);

fileID = fopen( 'Ymg_file.dat'); 
v = fread(fileID, inf, 'double');
Y = reshape(v, [Nx Ny]);

fileID = fopen( 'err_file.dat'); 
w = fread(fileID, inf, 'double');
err = reshape(w, [iter_max+1 1]);

fileID = fopen( 'grad_file.dat'); 
z = fread(fileID, inf, 'double');
grad = reshape(z, [iter_max+1 1]);

fileID = fopen( 'tau_n_file.dat'); 
tau_n = fread(fileID, inf, 'double');

save('test_fortran_data2.mat ', 'ntrial', 'Nx', 'Ny', 'Nt', 'n_obs_x','n_obs_y','x0_min', 'dmu_x','dmu_y', 'iter_max', ...
    'obs_vals','eta_optimum', 'eta_exct0', 'eta0', 'X', 'Y', 'err', 'grad')
