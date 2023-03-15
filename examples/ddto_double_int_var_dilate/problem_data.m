function prb = problem_data(K,scp_iters,wvc,wtr,cost_factor)
    
    prb.K = K;

    prb.ntarg = 3;                                  % Number of targets
    
    prb.n = 3;                                      % Dimension of double integrator

    prb.nx = 2*prb.n*prb.ntarg;
    prb.nu = (prb.n+1)*prb.ntarg;
    prb.np = 0;
    
    prb.tau = grid.generate_grid(0,1,K,'uniform');  % Generate grid in [0,1]

    prb.dtau = diff(prb.tau);
    
    prb.h = (1/5)*min(prb.dtau);                    % Step size for integration that computes FOH matrices
    prb.Kfine = 1+round(2/min(prb.dtau));           % Size of grid on which SCP solution is simulated
    
    % Deferrability index
    prb.Kstr = 5;

    % System parameters

    prb.g = [0;0;0];                                % External acceleration vector
        
    prb.c_d = 0.3;                                  % Drag coefficient
    
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
    prb.v1 = [0.1;0;0];
    
    prb.rK = [5  10   0;
              5   0  10;
              2   1   0];
    prb.vK = -0.1*[1  0  0;
                   0  1  0;
                   0  0  1];

    assert(length(prb.r1) == prb.n && length(prb.v1) == prb.n,"Specified initial position or velocity does match the system dimension.");
    assert(size(prb.rK,2) >= prb.ntarg && size(prb.vK,2) >= prb.ntarg,"Insufficient targets states specified.");

    prb.x1 = repmat([prb.r1;prb.v1],[prb.ntarg,1]);
    prb.xK = reshape([prb.rK;prb.vK],[prb.nx,1]);    
    prb.u1 = repmat([ones(3,1);20/K],[prb.ntarg,1]);
    prb.uK = repmat([ones(3,1);20/K],[prb.ntarg,1]);

    % Scaling parameters
    xmin = repmat([-0.5*prb.rmax*ones(3,1); -0.5*prb.vmax*ones(3,1)],[prb.ntarg,1]);
    xmax = repmat([ 0.5*prb.rmax*ones(3,1);  0.5*prb.vmax*ones(3,1)],[prb.ntarg,1]);
    
    umin = repmat([zeros(3,1); 1],[prb.ntarg,1]);
    umax = repmat([prb.umax*ones(3,1); 5],[prb.ntarg,1]);

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
    prb.time_of_maneuver = @(x,u) disc.time_of_maneuver(prb.disc,prb.tau(1:prb.Kstr),u(prb.n+1,1:prb.Kstr)); % Returns the time available to defer decision
    prb.time_grid = @(tau,x,s) disc.time_grid(prb.disc,tau,s);    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func = @(t,x,u) evaluate_dyn_func(x,u,prb.ntarg,prb.n,prb.c_d,prb.g);
    prb.dyn_func_linearize = @(tbar,xbar,ubar) evaluate_linearization(xbar,ubar,prb.ntarg,prb.n,prb.c_d,prb.g);

end

function dx = evaluate_dyn_func(x,u,ntarg,n,c_d,g)
    x = reshape(x,[2*n,ntarg]);
    u = reshape(u,[n+1,ntarg]);
    dx = zeros(2*n,ntarg);
    for j = 1:ntarg
        dx(:,j) = plant.doubleint.dyn_func(x(:,j),u(1:n,j),u(n+1,j),n,c_d,g);
    end
    dx = reshape(dx,[2*n*ntarg,1]);
end

function [A,B,w] = evaluate_linearization(x,u,ntarg,n,c_d,g)
    x = reshape(x,[2*n,ntarg]);
    u = reshape(u,[n+1,ntarg]);
    A = [];
    B = [];
    w = [];
    for j = 1:ntarg
        [Aj,Bj,Sj,wj] = plant.doubleint.compute_linearization(x(:,j),u(1:n,j),u(n+1,j),n,c_d,g);
        A = blkdiag(A,Aj);
        B = blkdiag(B,[Bj,Sj]);
        w = [w;wj];
    end
end