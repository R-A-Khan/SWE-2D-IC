% Run data assimilation for various cases and plot results
clear all
% Parameters for the data assimilation calculations
n_obs_x        = [2 3 4 5];          % Number of observation points
n_obs_y        = [2 3 4 5];          % Number of observation points
x_equals_y     = true;       % n_obs_x measurement points chosen along x = y
rand_x0        = false;      % Use random observation points
x0_min         = 0.2;        % x-Location of first observation point
y0_min         = 0.2;        % y-Location of first observation point
dmu_x          = 0.2;        % Spacing between x coord of observation points
dmu_y          = 0.2;        % Spacing between y coord of observation points
ntrial         = 1;          % Number of trials
Nx             = 256;        % Number of y grid points
Ny             = 256;        % Number of y grid points
tmin           = 0;          % Starting time
tmax           = 2;          % control time
xmax           = 3;          % domain size: xmin = -xmax
ymax           = 3;          % domain size: ymin = -ymax
iter_max       = 500;        % Max number of iterations; 
line_min       = true;       % Find optimal step
cfl            = 1;          % CFL constant
nu             = 0;          % Kinematic viscosity co-efficient                 
smooth_grad    = false;
filt           = 0.9;
conj_grad_type = 0;          % 0 = Steepest Descent, 1 = Fletcher Reeves, 2 = Polak-Ribière
iter_chunk     = 10;         % No. of iterations after which descent algorithm resets to steepest descent 


disp(' ');

disp('Nonlinear SWE');
if x_equals_y
    for j_obs_x = n_obs_x
        k = find(n_obs_x==j_obs_x);
        j_obs_y = j_obs_x;
    %         dmu  = abs((xmax-1e-4) - x0_min)/n_obs(k);
        [eta_optimum(:,:,k), eta_exct0(:,:,k), X, Y, err(:,k), err_std(:,k), grad(:,k), grad_std(:,k), cost(:,k), cost_std, n_iter] = ...
    data_assimil_2D(Nx,Ny,j_obs_x,j_obs_y, x_equals_y,rand_x0,x0_min, y0_min,dmu_x,dmu_y,ntrial,tmin,tmax,xmax,ymax,iter_max,line_min,smooth_grad,filt, conj_grad_type, iter_chunk, nu);

    end
end
disp(' ');

% Save results
save('Assimilation', 'ntrial', 'n_obs_x','n_obs_y', 'rand_x0','x0_min', 'dmu_x','dmu_y', 'line_min','smooth_grad', 'iter_max', ...
    'eta_optimum', 'eta_exct0', 'X', 'Y', 'err', 'err_std', 'grad', 'grad_std', 'cost', 'cost_std');
%% Plotting
clear all;
%

load('Assimilation');
clf;
iter = 0:iter_max;
marker1 = {'k-+', 'k-o', 'k-*', 'k-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'k--+','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

skip=iter_max/25;

set(0,'DefaultFigureWindowStyle','docked')

figure(1);
for j_obs = n_obs_x
    k = find(n_obs_x==j_obs);
    h(k)=semilogy(iter(1:skip:end), err(1:skip:end,k), marker1{k},'markersize',10,'linewidth',1.5);hold on;
end
set(gca,'FontSize',16); grid on;
legend(h, {'N_{obs} = 2','N_{obs} = 3','N_{obs} = 4','N_{obs} = 5','N_{obs} = 6'},'location','northeastoutside');
% legend(h, {'N_{obs} = 2','N_{obs} = 3', 'N_{obs} = 4'});
xlabel('Iteration n'); ylabel('$|\eta_0^{(n)}(x,y) - f(x,y)|_2/|f(x,y)|_2$','interpreter','latex');xlim([0 iter_max]);
ylim([1e-3 1e0]);
% print -depsc2 err_assimil.eps
hold off;

figure(2);
for j_obs = n_obs_x
    k = find(n_obs_x==j_obs);
    h(k)=semilogy(iter(1:skip:end), grad(1:skip:end,k)/grad(1,k), marker1{k},'markersize',10,'linewidth',1.5);hold on;
end
set(gca,'FontSize',16); grid on;
legend(h, {'N_{obs} = 2','N_{obs} = 3','N_{obs} = 4','N_{obs} = 5','N_{obs} = 6'},'location','northeastoutside');
xlabel('Iteration n'); ylabel('$|\nabla J^{(n)}(x,y)|_2/|\nabla J^{(0)}(x,y)|_2$','interpreter','latex');xlim([0 iter_max]);
% print -depsc2 err_grad.eps
hold off;

figure(3);
for j_obs = n_obs_x
    k = find(n_obs_x==j_obs);
    h(k) = semilogy(iter(1:skip:end), cost(1:skip:end,k)/cost(1,k), marker1{k},'markersize',10,'linewidth',1.5);hold on;
end
set(gca,'FontSize',16); grid on;
legend(h, {'N_{obs} = 2','N_{obs} = 3','N_{obs} = 4','N_{obs} = 5','N_{obs} = 6'},'location','northeastoutside');
xlabel('Iteration n'); ylabel('$J^{(n)}/J^{(0)}$','interpreter','latex');xlim([0 iter_max]);
% print -depsc2 err_cost.eps
hold off;

figure(4);
surf(X,Y,eta_optimum(:,:,end));grid on; 
set(gca,'FontSize',16);
xlabel('x'); ylabel('y'); zlabel('$\eta_0(x,y)$', 'interpreter','latex');%xlim([xmin xmax]); ylim([ymin ymax]
% print -depsc2 optimum_eta.eps


% figure(5)
% for j_obs = n_obs
%     k = find(n_obs==j_obs);
%     semilogy(j_obs, err_std(end,k,1)/err(end,k,1)*100, marker1{k},'markersize',10,'linewidth',1.5);hold on;
%     semilogy(j_obs, err_std(end,k,2)/err(end,k,2)*100, marker2{k},'markersize',10,'linewidth',1.5);
% end
% xlabel('N_{obs}');ylabel('% error');
% for j_obs = n_obs
%     k = find(n_obs==j_obs);
%     disp([sprintf('%1.f',j_obs),' ',sprintf('%0.3e',err(end,k,1)), ' ', sprintf('%0.3e',err_std(k,1)), ...
%        ' ',sprintf('%0.3e',err(end,k,2)), ' ', sprintf('%0.3e',err_std(k,2))]);
% end

% figure(5); 
% semilogy(X, abs(eta_optimum(:,1)-eta_exct0(:,1)), 'k-');hold on;
% semilogy(X, abs(eta_optimum(:,2)-eta_exct0(:,2)), 'k--');
% set(gca,'FontSize',16); grid on;
% xlabel('x'); ylabel('Error'); xlim([-xmax xmax]);