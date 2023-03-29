function prb = problem_data(K,scp_iters,wvc,wtr,cost_factor)
    
    prb.K = K;

    % Dimension of double integrator
    prb.n = 2;

    % Number of integrated constraints
    prb.m = 4;

    prb.nx = 2*prb.n + prb.m;
    prb.nu = prb.n+1;
    prb.np = 0;
    
    prb.tau = grid.generate_grid(0,1,K,'uniform'); % Generate grid in [0,1]

    prb.dtau = diff(prb.tau);
    
    prb.h = (1/20)*min(prb.dtau);            % Step size for integration that computes FOH matrices
    prb.Kfine = 1+20*round(1/min(prb.dtau));    % Size of grid on which SCP solution is simulated
    
    % System parameters

    prb.g = [zeros(prb.n-1,1);-1];
    
    prb.c_d = 0.1; 
    
    % Bounds

    prb.rmax = 40;
    prb.vmax = 5;
    prb.umax = 5;
    prb.umin = 2;

    prb.smin = 0.01;
    prb.smax = 10;
    prb.dtmin = 0.1;
    prb.dtmax = 3;
    prb.ToFmax = 20;

    prb.betmin = zeros(prb.m,1);
    prb.betmax = [1;
                  1;
                  1;
                  1];

    % Obstacle avoidance
    prb.nobs = 2;

    prb.robs = [-5 -10;
                 6  20];
    prb.aobs = [6 7];

    % Boundary conditions

    prb.r1 = [0;0];
    prb.v1 = [0;0];
    prb.bet1 = zeros(prb.m,1);
    
    prb.rK = [-15;28];
    prb.vK = -0.1*[1;0];
    prb.betK = zeros(prb.m,1);    

    assert(length(prb.r1) == prb.n && length(prb.v1) == prb.n,"Specified initial position or velocity does match the system dimension.");

    prb.x1 = [prb.r1;prb.v1;prb.bet1];
    prb.xK = [prb.rK;prb.vK;prb.betK];    
    prb.u1 = [ones(prb.n,1);20/K];
    prb.uK = [ones(prb.n,1);20/K];

    % Scaling parameters
    xmin = [-0.5*prb.rmax*ones(prb.n,1); -0.5*prb.vmax*ones(prb.n,1); prb.betmin];
    xmax = [ 0.5*prb.rmax*ones(prb.n,1);  0.5*prb.vmax*ones(prb.n,1); prb.betmax];
    
    umin = [zeros(prb.n,1); 1];
    umax = [prb.umax*ones(prb.n,1); 5];

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[0,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    % Scaled constraint parameters

    prb.cnstr_fun       = @(x,u) [-norm(x(1:prb.n)-prb.robs(:,1)) + prb.aobs(1);
                                  -norm(x(1:prb.n)-prb.robs(:,2)) + prb.aobs(2);
                                   norm(x(prb.n+1:2*prb.n))^2 - prb.vmax^2;
                                  -norm(u(1:prb.n)) + prb.umin];

    prb.cnstr_fun_jac_x = @(x,u) [-(x(1:prb.n)-prb.robs(:,1))'/norm(x(1:prb.n)-prb.robs(:,1)) zeros(1,prb.n);
                                  -(x(1:prb.n)-prb.robs(:,2))'/norm(x(1:prb.n)-prb.robs(:,2)) zeros(1,prb.n);
                                  zeros(1,prb.n) 2*x(prb.n+1:2*prb.n)'
                                  zeros(1,2*prb.n)];

    prb.cnstr_fun_jac_u = @(x,u) [zeros(1,prb.n+1);
                                  zeros(1,prb.n+1);
                                  zeros(1,prb.n+1);
                                  -u(1:prb.n)'/norm(u(1:prb.n)) 0];    

    % SCP parameters

    prb.disc = "FOH";
    prb.foh_type = "v3";
    prb.scp_iters = scp_iters; % Maximum SCP iterations

    prb.solver_settings = sdpsettings('solver','ecos','verbose',false);
    
    prb.tr_norm = 2;
    % prb.tr_norm = inf;
    % prb.tr_norm = 'quad';
    
    prb.wvc = wvc;
    prb.wtr = wtr;
    prb.cost_factor = cost_factor;
    
    prb.epsvc = 1e-8;
    prb.epstr = 1e-7;

    % Takes in unscaled data
    prb.time_of_maneuver = @(z,u) disc.time_of_maneuver(prb.disc,prb.tau,u(prb.n+1,:));    
    prb.time_grid = @(tau,z,u) disc.time_grid(prb.disc,tau,u(prb.n+1,:));    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func = @(t,z,u) evaluate_dyn_func(z,u,prb.n,prb.c_d,prb.g,prb.cnstr_fun);
    prb.dyn_func_linearize = @(tbar,zbar,ubar) evaluate_linearization(zbar,ubar,prb.n,prb.c_d,prb.g,prb.cnstr_fun,...
                                               prb.cnstr_fun_jac_x,prb.cnstr_fun_jac_u);

end

function dz = evaluate_dyn_func(z,u,n,c_d,g,cnstr_fun)
    x = z(1:2*n);
    cnstr_val = cnstr_fun(x,u);
    dx = plant.doubleint.dyn_func(x,u(1:n),u(n+1),n,c_d,g);
    dz = [dx;
          arrayfun(@(y) max(0,y),cnstr_val) .^ 2];
end

function [A,B,w] = evaluate_linearization(z,u,n,c_d,g,cnstr_fun,cnstr_fun_jac_x,cnstr_fun_jac_u)

    x = z(1:2*n);
    [dfdx,dfdu,dfds] = plant.doubleint.compute_linearization(x,u(1:n),u(n+1),n,c_d,g);
    dfdu = [dfdu,dfds];

    cnstr_val = cnstr_fun(x,u);
    cnstr_val_jac_x = cnstr_fun_jac_x(x,u);
    cnstr_val_jac_u = cnstr_fun_jac_u(x,u);

    m = length(cnstr_val);
    absg2 = zeros(m,1);
    absg2_jac_x = zeros(m,2*n);
    absg2_jac_u = zeros(m,n+1);
    for j = 1:m
        absg2(j) = max(0,cnstr_val(j))^2;
        if cnstr_val(j) > 0
            absg2_jac_x(j,:) = 2*cnstr_val(j)*cnstr_val_jac_x(j,:);
            absg2_jac_u(j,:) = 2*cnstr_val(j)*cnstr_val_jac_u(j,:);
        end
    end

    A = [dfdx        zeros(2*n,m);
         absg2_jac_x zeros(m,m)];
    B = [dfdu;
         absg2_jac_u];
    h = evaluate_dyn_func(z,u,n,c_d,g,cnstr_fun);
    w = h - A*z - B*u;
end