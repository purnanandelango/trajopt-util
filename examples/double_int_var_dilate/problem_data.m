function prb = problem_data(K,scp_iters,wvc,wtr,cost_factor)
    
    prb.K = K;

    % Dimension of double integrator
    prb.n = 2;

    prb.nx = 2*prb.n;
    prb.nu = prb.n+1;
    prb.np = 0;
    
    prb.tau = grid.generate_grid(0,1,K,'uniform'); % Generate grid in [0,1]

    prb.dtau = diff(prb.tau);
    
    prb.h = (1/100)*min(prb.dtau);            % Step size for integration that computes FOH matrices
    prb.Kfine = 1+100*round(1/min(prb.dtau));    % Size of grid on which SCP solution is simulated
    
    % System parameters

    prb.g = [zeros(prb.n-1,1);-1];
    
    prb.c_d = 0.1; 
    
    % Bounds

    prb.rmax = 40;
    prb.vmax = 7;
    prb.umax = 8;
    prb.umin = 2;

    prb.smin = 1;
    prb.smax = 15;
    prb.dtmin = 0.1;
    prb.dtmax = 3;
    prb.ToFmax = 20;
    prb.ToFguess = 10;

    % Obstacle avoidance
    prb.nobs = 2;

    prb.robs = [-5 -10;
                 6  20];
    prb.aobs = [6 7];

    % Boundary conditions

    prb.r1 = [0;0];
    prb.v1 = [0;0];
    
    prb.rK = [-15;28];
    prb.vK = -0.1*[1;0];

    assert(length(prb.r1) == prb.n && length(prb.v1) == prb.n,"Specified initial position or velocity does match the system dimension.");

    prb.x1 = [prb.r1;prb.v1];
    prb.xK = [prb.rK;prb.vK];    
    prb.u1 = [ones(prb.n,1);prb.ToFguess];
    prb.uK = [ones(prb.n,1);prb.ToFguess];

    % Scaling parameters
    xmin = [-0.5*prb.rmax*ones(prb.n,1); -0.5*prb.vmax*ones(prb.n,1)];
    xmax = [ 0.5*prb.rmax*ones(prb.n,1);  0.5*prb.vmax*ones(prb.n,1)];
    
    umin = [zeros(prb.n,1); 1];
    umax = [prb.umax*ones(prb.n,1); 5];

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[-1,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    % SCP parameters

    prb.disc = "FOH";
    prb.foh_type = "v1";
    prb.ode_solver = {'ode45',odeset('RelTol',1e-5,'AbsTol',1e-7)};
    prb.scp_iters = scp_iters; % Maximum SCP iterations

    prb.solver_settings = sdpsettings('solver','ecos','verbose',false);
    
    % prb.tr_norm = 2;
    % prb.tr_norm = inf;
    prb.tr_norm = 'quad';
    
    prb.wvc = wvc;
    prb.wtr = wtr;
    prb.cost_factor = cost_factor;
    
    prb.epsvc = 1e-7;
    prb.epstr = 1e-3;

    % Takes in unscaled data
    prb.time_of_maneuver = @(x,u) disc.time_of_maneuver(prb.disc,prb.tau,u(prb.n+1,:));    
    prb.time_grid = @(tau,x,u) disc.time_grid(prb.disc,tau,u(prb.n+1,:));    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func = @(t,x,u) plant.doubleint.dyn_func(x,u(1:prb.n),u(prb.n+1),prb.n,prb.c_d,prb.g);
    prb.dyn_func_linearize = @(tbar,xbar,ubar) evaluate_linearization(xbar,ubar,prb.n,prb.c_d,prb.g);

end

function [A,B,w] = evaluate_linearization(x,u,n,coeff_drag,g)
    [A,B,S,w] = plant.doubleint.compute_linearization(x,u(1:n),u(n+1),n,coeff_drag,g);
    B = [B,S];
end