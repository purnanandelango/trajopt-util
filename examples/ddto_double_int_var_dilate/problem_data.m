function prb = problem_data(K,scp_iters,wvc,wtr,cost_factor)
    
    prb.K = K;

    ntarg = 2;
    prb.ntarg = ntarg;

    prb.nx = 6*ntarg;
    prb.nu = 4*ntarg;
    prb.np = 0;
    
    prb.tau = grid.generate_grid(0,1,K,'uniform'); % Generate grid in [0,1]

    prb.dtau = diff(prb.tau);
    
    prb.h = (1/5)*min(prb.dtau);            % Step size for integration that computes FOH matrices
    prb.Kfine = 1+round(2/min(prb.dtau));    % Size of grid on which SCP solution is simulated
    
    % Deferrability index
    prb.taustr = 5;


    % System parameters

    prb.g = [0;0;0];
    
    prb.c_d = 0.3; 
    
    % Bounds

    prb.rmax = 10;
    prb.vmax = 5;
    prb.umax = 8;

    prb.smin = 0.01;
    prb.smax = 10;
    prb.dtmin = 0.01;
    prb.dtmax = 2;
    prb.ToFmax = 20;
    
    % Boundary conditions

    prb.r1 = [0;0;0];
    prb.v1 = [0.5;0;0];
    
    prb.rK = [5  10;
              5  0;
              2  1];
    prb.vK = -0.1*[1  0;
                   0  1;
                   0  0];

    prb.x1 = [prb.r1;prb.v1;prb.r1;prb.v1];
    prb.xK = [prb.rK(:,1);prb.vK(:,1);prb.rK(:,2);prb.vK(:,2)];    
    prb.u1 = [ones(3,1);20/K;ones(3,1);20/K];
    prb.uK = [ones(3,1);20/K;ones(3,1);20/K];

    % Scaling parameters
    xmin = [-0.5*prb.rmax*ones(3,1); -0.5*prb.vmax*ones(3,1); -0.5*prb.rmax*ones(3,1); -0.5*prb.vmax*ones(3,1)];
    xmax = [ 0.5*prb.rmax*ones(3,1);  0.5*prb.vmax*ones(3,1);  0.5*prb.rmax*ones(3,1);  0.5*prb.vmax*ones(3,1)];
    
    umin = [zeros(3,1); 1; zeros(3,1); 1];
    umax = [prb.umax*ones(3,1); 5; prb.umax*ones(3,1); 5];

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[-1,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

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
    prb.time_of_maneuver = @(x,u) 0.5 * ( disc.time_of_maneuver(prb.disc,prb.tau,u(4,:)) + disc.time_of_maneuver(prb.disc,prb.tau,u(8,:)) );    
    prb.time_grid = @(tau,x,s) disc.time_grid(prb.disc,tau,s);    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func = @(t,x,u) evaluate_dyn_func(x,u,uint64(prb.nx/(2*prb.ntarg)),prb.c_d,prb.g);
    prb.dyn_func_linearize = @(tbar,xbar,ubar) evaluate_linearization(xbar,ubar,uint64(prb.nx/(2*prb.ntarg)),prb.c_d,prb.g);

end

function dx = evaluate_dyn_func(x,u,n,coeff_drag,g)
    x1 = x(1:2*n);
    x2 = x(2*n+1:4*n);
    u1 = u(1:n+1);
    u2 = u(n+2:2*n+2);
    dx1 = plant.doubleint.dyn_func(x1,u1(1:n),u1(n+1),n,coeff_drag,g);
    dx2 = plant.doubleint.dyn_func(x2,u2(1:n),u2(n+1),n,coeff_drag,g);
    dx = [dx1;dx2];
end

function [A,B,w] = evaluate_linearization(x,u,n,coeff_drag,g)
    x1 = x(1:2*n);
    x2 = x(2*n+1:4*n);
    u1 = u(1:n+1);
    u2 = u(n+2:2*n+2);
    
    [A1,B1,S1,w1] = plant.doubleint.compute_linearization(x1,u1(1:n),u1(n+1),n,coeff_drag,g);
    B1 = [B1,S1];
    [A2,B2,S2,w2] = plant.doubleint.compute_linearization(x2,u2(1:n),u2(n+1),n,coeff_drag,g);
    B2 = [B2,S2];

    A = blkdiag(A1,A2);
    B = blkdiag(B1,B2);
    w = [w1;w2];
end