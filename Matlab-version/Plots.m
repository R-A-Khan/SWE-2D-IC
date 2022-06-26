load('Assimilation_DA2D8iii.mat')
err_8iii = err;
grad_8iii = grad;
cost_8iii = cost;
obs_vals_8iii = obs_vals;
clear('Assimilation_DA2D8iii.mat')

load('Assimilation_DA2D8iv.mat')
err_8iv = err;
grad_8iv = grad;
cost_8iv = cost;
obs_vals_8iv = obs_vals;




iter = 1:iter_max+1;
marker1 = {'k-+', 'b-+', 'r-+', 'g-x', 'm-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'k--+','b--+','r--*','g--x','m--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

skip=iter_max/20;

set(0,'DefaultFigureWindowStyle','docked')

%%

figure(1);
for j_obs = n_obs_x
    k = find(n_obs_x==j_obs);
    h(k)=semilogy(iter(1:skip:end), err_8iii(1:skip:end), marker1{k},'markersize',10,'linewidth',1.5);hold on;
         semilogy(iter(1:skip:end), err_8iv(1:skip:end), marker1{k+1},'markersize',10,'linewidth',1.5);
end
set(gca,'FontSize',20); grid on;
legend( 'Single Sided', 'Double Sided','location','northeast');
% legend(h, {'N_{obs} = 2','N_{obs} = 3', 'N_{obs} = 4'});
xlabel('Iteration n', 'interpreter','latex');
ylabel('$||\phi^{(t)} - \phi^{(b)}||_2/||\phi^{(t)}||_2$','interpreter','latex');
xlim([0 iter_max]);
axis square
str3 = 'Plots/DA2D8/Figures/D8_err_3_4.eps';
print(str3, '-depsc2')
str5 = 'Plots/DA2D8/Figures/D8_err_3_4.eps.fig';
savefig(str5)
hold off;

figure(2);
for j_obs = n_obs_x
    k = find(n_obs_x==j_obs);
    h(k)=semilogy(iter(1:skip:end), cost_8iii(1:skip:end)/cost_8iii(1), marker1{k},'markersize',10,'linewidth',1.5);hold on;
         semilogy(iter(1:skip:end), cost_8iv(1:skip:end)/cost_8iv(1), marker1{k+1},'markersize',10,'linewidth',1.5); 
end
set(gca,'FontSize',20); grid on;
legend( 'Single Sided', 'Double Sided','location','northeast');
xlabel('Iteration n', 'interpreter','latex'); ylabel('$|J^{(n)}|_2/| J^{(0)}|_2$','interpreter','latex');
xlim([0 iter_max]);
% print -depsc2 err_grad.eps
hold off;
axis square
str3 = 'Plots/DA2D8/Figures/D8_cost_3_4.eps';
print(str3, '-depsc2')
str5 = 'Plots/DA2D8/Figures/D8_cost_3_4.fig';
savefig(str5)



figure(3);
% scatter(obs_vals_8iii(:,1),obs_vals_8iii(:,2),'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',[0 .7 .7],...
%               'LineWidth',1.5); hold on
scatter(x0_vals, y0_vals,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5); hold on
xlim([-3,3]);
xlabel('$x$','interpreter','latex')
ylim([-3,3]);
ylabel('$y$','interpreter','latex')
set(gca,'FontSize',20); grid on

% s = surf(X,Y,eta_optimum,'FaceAlpha',0.5)
% s.EdgeColor = 'none'
% view(0,90)
scatter(0,0, 1000,'LineWidth',1.5); 
hold off
axis square
print( 'Plots/DA2D8/Figures/D8_obs_vals_3.eps', '-depsc2')
savefig('Plots/DA2D8/Figures/D8_obs_vals_3.fig')


figure(4);
% scatter(obs_vals_8iv(:,1),obs_vals_8iv(:,2),'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',[0 .7 .7],...
%               'LineWidth',1.5); hold on
scatter(x0_vals, y0_vals,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5); hold on
xlim([-3,3]);
xlabel('$x$','interpreter','latex')
ylim([-3,3]);
ylabel('$y$','interpreter','latex')
set(gca,'FontSize',20); grid on;
scatter(0,0, 1000,'LineWidth',1.5); 
hold off
axis square
print( 'Plots/DA2D8/Figures/D8_obs_vals_4.eps', '-depsc2')
savefig('Plots/DA2D8/Figures/D8_obs_vals_4.fig')


%% 
clear all

load('Plots/DA2D8/Assimilation_DA2D8i_v2.mat')
err_8i = err;
grad_8i = grad;
cost_8i = cost;
obs_vals_8i = obs_vals;
clear('Plots/DA2D8/Assimilation_DA2D8i_v2.mat')

load('Plots/DA2D8/Assimilation_DA2D8ii_v2.mat')
err_8ii = err;
grad_8ii = grad;
cost_8ii = cost;
obs_vals_8ii = obs_vals;



iter = 1:iter_max+1;
marker1 = {'k-+', 'b-+', 'r-+', 'g-x', 'm-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'k--+','b--+','r--*','g--x','m--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

skip=iter_max/20;

set(0,'DefaultFigureWindowStyle','docked')

%%

figure(1);
for j_obs = n_obs_x
    k = find(n_obs_x==j_obs);
    h(k)=semilogy(iter(1:skip:end), err_8i(1:skip:end), marker1{k},'markersize',10,'linewidth',1.5);hold on;
         semilogy(iter(1:skip:end), err_8ii(1:skip:end), marker1{k+1},'markersize',10,'linewidth',1.5);
end
set(gca,'FontSize',20); grid on;
legend( 'Single Sided', 'Double Sided','location','northeast');
% legend(h, {'N_{obs} = 2','N_{obs} = 3', 'N_{obs} = 4'});
xlabel('Iteration n', 'interpreter','latex');
ylabel('$||\phi^{(t)} - \phi^{(b)}||_2/||\phi^{(t)}||_2$','interpreter','latex');
xlim([0 iter_max]);
axis square
str3 = 'Plots/DA2D8/Figures/D8_err_1_2_v2.eps';
print(str3, '-depsc2')
str5 = 'Plots/DA2D8/Figures/D8_err_1_2_v2.fig';
savefig(str5)
hold off;

figure(2);
for j_obs = n_obs_x
    k = find(n_obs_x==j_obs);
    h(k)=semilogy(iter(1:skip:end), cost_8i(1:skip:end)/cost_8i(1), marker1{k},'markersize',10,'linewidth',1.5);hold on;
         semilogy(iter(1:skip:end), cost_8ii(1:skip:end)/cost_8ii(1), marker1{k+1},'markersize',10,'linewidth',1.5); 
end
set(gca,'FontSize',20); grid on;
legend( 'Single Sided', 'Double Sided','location','northeast');
xlabel('Iteration n', 'interpreter','latex'); ylabel('$ |J^{(n)}|_2/|J^{(0)}|_2$','interpreter','latex');
xlim([0 iter_max]);
axis square
% print -depsc2 err_grad.eps
hold off;
str3 = 'Plots/DA2D8/Figures/D8_cost_1_2_v2.eps';
print(str3, '-depsc2')
str5 = 'Plots/DA2D8/Figures/D8_cost_1_2_v2.fig';
savefig(str5)
hold off;



figure(3);
scatter(obs_vals_8i(:,1),obs_vals_8i(:,2),'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5); hold on
xlim([-3,3]);
xlabel('$x$','interpreter','latex')
ylim([-3,3]);
ylabel('$y$','interpreter','latex')
set(gca,'FontSize',20); grid on

% s = surf(X,Y,eta_optimum,'FaceAlpha',0.5)
% s.EdgeColor = 'none'
% view(0,90)
scatter(0,0, 1000,'LineWidth',1.5); 
hold off
axis square
print( 'Plots/DA2D8/Figures/D8_obs_1_v2.eps', '-depsc2')
savefig('Plots/DA2D8/Figures/D8_obs_1_v2.fig')


figure(4);
scatter(obs_vals_8ii(:,1),obs_vals_8ii(:,2),'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5); hold on
xlim([-3,3]);
xlabel('$x$','interpreter','latex')
ylim([-3,3]);
ylabel('$y$','interpreter','latex')
set(gca,'FontSize',20); grid on;
scatter(0,0, 1000,'LineWidth',1.5); 
hold off
axis square
print( 'Plots/DA2D8/Figures/D8_vals_2_v2.eps', '-depsc2')
savefig('Plots/DA2D8/Figures/D8_vals_2_v2.fig')
