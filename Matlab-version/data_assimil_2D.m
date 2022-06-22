function [eta_optimum, eta_exct0, X, Y, err, err_std, grad, grad_std, cost, cost_std, n_iter] = ...
    data_assimil_2D(Nx,Ny,n_obs_x,n_obs_y, x_equals_y,rand_x0,x0_min, y0_min,dmu_x,dmu_y,ntrial,tmin,tmax,xmax,ymax,iter_max,line_min,smooth_grad,filt, conj_grad_type, iter_chunk, nu)

% n_obs_x        = 5;          % Number of observation points
% n_obs_y        = 5;          % Number of observation points
% x_equals_y     = true;       % n_obs_x measurement points chosen along x = y
% rand_x0        = false;      % Use random observation points
% x0_min         = 0.2;        % x-Location of first observation point
% y0_min         = 0.2;        % y-Location of first observation point
% dmu_x          = 0.2;        % Spacing between x coord of observation points
% dmu_y          = 0.2;        % Spacing between x coord of observation points
% ntrial         = 1;          % Number of trials
% Nx             = 128;        % Number of y grid points
% Ny             = 128;        % Number of y grid points
% tmin           = 0;          % Starting time
% tmax           = 2;          % control time
% xmax           = 3;          % domain size: xmin = -xmax
% ymax           = 3;          % domain size: ymin = -ymax
% iter_max       = 1000;        % Max number of iterations; 
% line_min       = true;       % Find optimal step
% cfl            = 1;          % CFL constant
% nu             = 0;          % Kinematic viscosity co-efficient                 
% smooth_grad    = false;
% filt           = 0.9;
% conj_grad_type = 0;          % 0 = Steepest Descent, 1 = Fletcher Reeves, 2 = Polak-Ribière
% iter_chunk     = 10;         % No. of iterations after which descent algorithm resets to steepest descent 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create grid
xmin = -xmax; ymin = -ymax;
delta_x = abs(xmax - xmin)/Nx; delta_y = abs(ymax - ymin)/Ny;
dt = min(delta_x,delta_y)/3;
x = xmin:delta_x:xmax - delta_x;
y = ymin: delta_y:ymax - delta_y;
T = tmin:dt:tmax;

[X,Y] = meshgrid(x,y);
X = X'; Y = Y';

% Observation points
x0_max = x0_min+(n_obs_x-1)*dmu_x;
if x0_max >= xmax
    dmu_x = (xmax-x0_min)/n_obs_x;
    x0_max = x0_min+(n_obs_x-1)*dmu_x;
    disp('X-distance between measurement points has been adjusted to fit domain')
else
    x0_max = x0_min+(n_obs_x-1)*dmu_x;
end

x0_vals = x0_min:dmu_x:x0_max; %  Evenly spaced
x0_vals = reshape(x0_vals, length(x0_vals),1);

y0_max = y0_min+(n_obs_y-1)*dmu_y;
if y0_max >= ymax
    dmu_y = (ymax-y0_min)/n_obs_y;
    y0_max = y0_min+(n_obs_y-1)*dmu_y;
    disp('Y-distance between measurement points has been adjusted to fit domain')
else
    y0_max = y0_min+(n_obs_y-1)*dmu_y;
end

y0_vals = y0_min:dmu_y:y0_max; % Evenly spaced
y0_vals = reshape(y0_vals, length(y0_vals),1);


x0_vals_all = repmat(x0_vals,n_obs_y,1);
y0_vals_all = repelem(y0_vals, n_obs_x);


if x_equals_y
    obs_vals = [ x0_vals y0_vals]; % Matrix of observation coordinates
    [x0_inds, x0_pts, y0_inds, y0_pts] = mult_x0_y0(x, y, x0_vals, y0_vals);
else
    obs_vals = [ x0_vals_all y0_vals_all]; % Matrix of observation coordinates
    [x0_inds, x0_pts, y0_inds, y0_pts] = mult_x0_y0(x, y, x0_vals_all, y0_vals_all);
end

obs_inds = [ x0_inds y0_inds];

fw = @SWE_NL_SadEnConsv_C_vec;
bw = @BW_SWE_NL_2D_vec;

for jtrial=1:ntrial
    disp(['Trial ', sprintf('%i', jtrial)]);

    % STEP I: Run forward solver using exact initial conditions to get observations
    amp =0.05;
    eta_exct0 = amp*exp(- ((X).^2 + (Y).^2)/0.1^2);
    u = zeros(Nx,Ny);
    v = zeros(Nx,Ny);
    H = [eta_exct0 ; u; v];  % initial data is 3*N x Ny size matrix
    [T_obs,sol_exct] = Solve_NL_SWE_2D( dt, [tmin tmax], fw, H, delta_x, delta_y, nu);
    % Exact surface wave at all x,y points for all t
    obs_eta = sol_exct(1:Nx,:,:);
    % Exact wave at measurement points x0,y0 for all t
    obs = zeros(length(x0_inds),length(T_obs));
    for i= 1:length(x0_inds)
        obs(i,:) = sol_exct(x0_inds(i),y0_inds(i),:);
    end



    % STEP 2: Distort exact IC to get "guess" IC eta0(:,1) used in algorithm with analytical grad_J
    eta0(:,:,1) = zeros(Nx, Ny,1);


    % Define eta0(:,2) as initial guess height computed using grad_J
    eta0_2 = zeros(Nx, Ny,1);
    eta0_2(:,:,1) = eta0(:,:,1);

    % Run forward and backwards solvers to get eta_a_t0_

    H2 = [eta0_2(:,:,1); u; v];  % initial data is 3*Nx x Ny x 1 size 3D matrix
    [T,sol] = Solve_NL_SWE_2D( dt, [tmin tmax], fw, H2, delta_x, delta_y, nu);
    eta_all = sol(1:Nx,:,:);
    u_all   = sol(Nx+1:2*Nx,:,:);
    v_all   = sol(2*Nx+1:3*Nx,:,:);
    eta_x0  = zeros(length(x0_inds),length(T));
    for i= 1:length(x0_inds)
        eta_x0(i,:) = sol(x0_inds(i),y0_inds(i),:);
    end

    H3 = zeros(3*Nx, Ny, 1);
    [eta_adj_t0, Tb, Yb] = BW_solve_2D(delta_x, delta_y, dt, H3, tmin, tmax, x0_inds, y0_inds, bw, obs_eta, u_all, v_all, eta_all);

    % Define grad_J
    grad_J = -eta_adj_t0;

    % Smooth gradient
    if smooth_grad
        for k = 1:Ny
            grad_J(:,k) = grad_smooth(grad_J(:,k),filt);
        end
        for l = 1:Nx
            grad_J(l,:) = grad_smooth(grad_J(l,:),filt);
        end
    end

    % Compute Cost Function
    cost_x = zeros(1,length(T));
    for j = 1:length(x0_inds)
        cost_x = cost_x + ( (obs(j,:) - eta_x0(j,:)).^2) ;
    end

    cost(1,jtrial) = trapz(T,0.5*cost_x) ;
    grad(1,jtrial) = norm(grad_J);
    err(1,jtrial)  = norm(eta_exct0-eta0)/norm(eta_exct0);

    % Initial estimate for gradient descent stepsize
    tau_n = 1e-2;
    %1/length(x0_vals);

    % BEGIN LOOP
    iter=0;
        while norm(grad_J) >= 1e-20 && iter<iter_max
            iter = iter+1;
            % STEP 3
            % Run forward solver to get height at x0 for all t given eta0_1
            H2 = [eta0_2(:,:,iter) ; u; v];        
            [T,sol] = Solve_NL_SWE_2D( dt, [tmin tmax], fw, H2, delta_x, delta_y, nu);
            eta_all = sol(1:Nx,:,:);
            u_all   = sol(Nx+1:2*Nx,:,:);
            v_all   = sol(2*Nx+1:3*Nx,:,:);
            eta_x0  = zeros(length(x0_inds),length(T));
            for i= 1:length(x0_inds)
                eta_x0(i,:) = sol(x0_inds(i),y0_inds(i),:);
            end

            % STEP 4
            % Run backward solver to get adjoint height at t0 for all x
            H3 = zeros(3*Nx, Ny, 1);
            [eta_adj_t0, Tb, Yb] = BW_solve_2D(delta_x, delta_y, dt, H3, tmin, tmax, x0_inds, y0_inds, bw, obs_eta, u_all, v_all, eta_all);
            % STEP 5: Define gradient of cost function
            grad_J = -eta_adj_t0;


            % Smooth gradient
            if smooth_grad
                for k = 1:Ny
                    grad_J(:,k) = grad_smooth(grad_J(:,k),filt);
                end
                for l = 1:Nx
                    grad_J(l,:) = grad_smooth(grad_J(l,:),filt);
                end
            end

            steepest_direction(:,:,iter) = - grad_J;
            % for the first iteration, this is steepest descent
            if mod(iter,iter_chunk) == 1 
                    % The conjugate direction is the steepest direction at
                    % the first iteration
                    conj_direction(:,:,iter) =  steepest_direction(:,:,iter); 

                     % STEP 6: Line Minimisation for optimal step size tau_n
                if line_min
                    warning('off','optim:fminunc:SwitchingMethod')
                    options = optimoptions(@fminunc,'Display', 'off', 'Tolfun',1e-6,'TolX', 1e-6);
                    tau0 = tau_n;
                    f = @(tau)line_min_mult_2D(obs, eta_adj_t0, Nx, Ny, delta_x, delta_y, x0_inds, y0_inds, tau, eta0_2(:,:,iter), dt, tmin, tmax, fw, nu);
                    tau_n = fminunc(f,tau0,options);
                    %tau_n = min(tau_n, 3/n_obs);
                end


                % STEP 7: Steepest Descent Algorithm
                eta0_2(:,:,iter+1)= eta0_2(:,:,iter) + tau_n*conj_direction(:,:,iter);
                eta_optimum = eta0_2(:,:,iter+1);

                % STEP 8: Run forward solver to get height at x0 for all t given optimised IC
                H4 = [eta0_2(:,:,iter+1);u; v];
                [T,sol_opt] = Solve_NL_SWE_2D( dt, [tmin tmax], fw, H4, delta_x, delta_y, nu);
                eta_x0_opt  = zeros(length(x0_inds),length(T));
                for i= 1:length(x0_inds)
                    eta_x0_opt(i,:) = sol_opt(x0_inds(i),y0_inds(i),:);
                end

                % Compute Cost Function
                cost_opt = zeros(1,length(T));
                for i = 1:length(x0_inds)
                    cost_opt = cost_opt + (obs(i,:) - eta_x0_opt(i,:)).^2 ;
                end

                cost(iter+1,jtrial) = trapz(T,0.5*cost_opt) ;
                grad(iter+1,jtrial) = norm(grad_J);
                err(iter+1,jtrial) = norm(eta_exct0-eta_optimum)/norm(eta_exct0);

                disp([sprintf('%1.f',iter),'  ',sprintf('%0.3e',tau_n), ' ', sprintf('%0.3e',norm(grad(iter,jtrial))),'  ', sprintf('%0.3e',norm(err(iter,jtrial)))]);


           else
                if conj_grad_type == 0
                   % Steepest Descent
                   conj_direction(:,:, iter) = steepest_direction(:,:,iter);         
                elseif conj_grad_type == 1 
                   % Fletcher-Reeves
                   b = (steepest_direction(:,:,iter).' * steepest_direction(:,:,iter)) / (steepest_direction(:,:,iter-1).' * steepest_direction(:,:,iter-1));
                   conj_direction(:,:,iter) = steepest_direction(:,:,iter) + b*conj_direction(:,:,iter-1);
                elseif conj_grad_type == 2
                   % Polak-Ribière
                   b = (steepest_direction(:,:,iter).' *  ((steepest_direction(:,:,iter)) - steepest_direction(:,:,iter-1) )) / (steepest_direction(:,:,iter-1).' * steepest_direction(:,:,iter-1));
                   conj_direction(:,:,iter) = steepest_direction(:,:,iter) + b*conj_direction(:,:,iter-1);
                end

                % STEP 6: Line Minimisation for optimal step size tau_n
                if line_min
                    warning('off','optim:fminunc:SwitchingMethod')
                    options = optimoptions(@fminunc,'Display', 'off', 'Tolfun',1e-6,'TolX', 1e-6);
                    tau0 = tau_n;
                    f = @(tau)line_min_mult_2D(obs, eta_adj_t0, Nx, Ny, delta_x, delta_y, x0_inds, y0_inds, tau, eta0_2(:,:,iter), dt, tmin, tmax, fw, nu);
                    tau_n = fminunc(f,tau0,options);
                    %tau_n = min(tau_n, 3/n_obs);
                end

                % STEP 7: Steepest Descent Algorithm
                eta0_2(:,:,iter+1)= eta0_2(:,:,iter) + tau_n*conj_direction(:,:,iter);
                eta_optimum = eta0_2(:,:,iter+1);

                % STEP 8: Run forward solver to get height at x0 for all t given optimised IC
                H4 = [eta0_2(:,:,iter+1);u;v];
                [T,sol_opt] = Solve_NL_SWE_2D( dt, [tmin tmax], fw, H4, delta_x, delta_y, nu);
                eta_x0_opt  = zeros(length(x0_inds),length(T));
                for i= 1:length(x0_inds)
                    eta_x0_opt(i,:) = sol_opt(x0_inds(i),y0_inds(i),:);
                end

                % Compute Cost Function
                cost_opt = zeros(1,length(T));
                for i = 1:length(x0_inds)
                    cost_opt = cost_opt + (obs(i,:) - eta_x0_opt(i,:)).^2 ;
                end

                cost(iter+1,jtrial) = trapz(T,0.5*cost_opt) ;
                grad(iter+1,jtrial) = norm(grad_J);
                err(iter+1,jtrial) = norm(eta_exct0-eta_optimum)/norm(eta_exct0);

                disp([sprintf('%1.f',iter),'  ',sprintf('%0.3e',tau_n), ' ', sprintf('%0.3e',norm(grad(iter,jtrial))),'  ', sprintf('%0.3e',norm(err(iter,jtrial)))]);
         end
      end
end

n_iter = iter+1;

err_std  = std(err(end,:));
err      = mean(err,2);

grad_std = std(grad(end,:));
grad     = mean(grad,2);

cost_std = std(cost(end,:));
cost     = mean(cost,2);

disp(['N_obs = ', sprintf('%i', n_obs_x), '  mean error = ', sprintf('%0.4e',err(end)),' std dev = ', sprintf('%0.4e',err_std(end))]);            