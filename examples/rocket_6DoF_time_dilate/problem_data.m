function prb = problem_data(K,scp_iters,w_ep,w_px,cost_factor)
    
    prb.K = K;
    prb.Kfine = prb.K*100; 

    prb.nx = 14;
    prb.nu = 3;
    prb.np = 1;
    prb.tau = ((1:prb.K)-1)/(prb.K-1); % Scaled time grid
    
    prb.h = (1/10)*1/(prb.K-1); % Step size for integration that computes FOH matrices
    
    % System parameters

    prb.g0 = 1;
    
    prb.c_ax = 0.5; prb.c_ayz = 1;
    % prb.c_ax = 0; prb.c_ayz = 0;
    
    prb.CA = diag([prb.c_ax,prb.c_ayz,prb.c_ayz]);
    
    prb.JB = 0.168*diag([2e-2,1,1]);
    prb.gI = [-prb.g0;0;0];
    prb.rho = 1;
    prb.SA = 0.5;
    
    prb.rTB = -0.25*[1;0;0];
    prb.rcpB = 0.05*[1;0;0];
    
    prb.Isp = 30;
    prb.alphmdt = 1/(prb.Isp*prb.g0);
    prb.betmdt = 0.01;
    
    % Bounds

    prb.thetmax = 45*pi/180; prb.sinthetmaxby2 = sin(prb.thetmax/2);
    prb.gamgs = 75*pi/180;   prb.cotgamgs = cot(prb.gamgs);
    
    prb.omgmax = 28.6*pi/180;
    prb.delmax = 20*pi/180; prb.cosdelmax = cos(prb.delmax);
    
    prb.Hgam = [0,1,0;0,0,1];
    prb.Hthet = [0,1,0,0;0,0,1,0];
    
    prb.Tmin = 1.5;
    prb.Tmax = 6.5;
    
    prb.Vmax = 3;
    
    prb.mdry = 1;
    prb.mwet = 2;
    
    % Boundary conditions

    prb.rI1 = [5.33;6.5;0];
    prb.vI1 = [-0.5;-2.5;0];
    % prb.vI1 = [-2.5;-0.5;0];
    
    prb.rIK = zeros(3,1);
    prb.vIK = -0.1*[1;0;0];
    prb.omgB1 = zeros(3,1);
    prb.omgBK = prb.omgB1;
    prb.q1 = [0;0;0;1];

    prb.x1 = [prb.mwet;prb.rI1;prb.vI1;prb.q1;prb.omgB1];
    prb.xK = [prb.mdry;prb.rIK;prb.vIK;prb.q1;prb.omgBK];    
    prb.u1 = -prb.mwet*prb.gI;
    prb.uK = -prb.mdry*prb.gI;
    prb.p = 5;

    % Scaling parameters
    xmin = [1;0;-2;-2;-3*ones(3,1);-ones(4,1);-prb.omgmax*ones(3,1)];
    xmax = [4;7; 2; 2; 3*ones(3,1); ones(4,1); prb.omgmax*ones(3,1)];
    umin = ones(prb.nu,1);
    umax = 6*ones(prb.nu,1);
    pmin = 0;
    pmax = 5;
    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax],[pmin,pmax]},[0,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.Sp = Sz{3}; prb.invSp = inv(Sz{3});
    prb.cx = cz{1};
    prb.cu = cz{2};
    prb.cp = cz{3};        

    % Takes in unscaled data
    prb.time_of_maneuver = @(x,u,p) p;

    % SCP parameters

    prb.disc = "ZOH";
    prb.foh_type = "v3";
    prb.ode_solver = {'ode45',odeset('RelTol',1e-5,'AbsTol',1e-7)};
    prb.scp_iters = scp_iters; % Maximum SCP iterations

    prb.solver_settings = sdpsettings('solver','ecos','verbose',false);
    % prb.solver_settings = sdpsettings('solver','ipopt','verbose',false);
    
    % prb.px_norm = 2;
    % prb.px_norm = inf;
    prb.px_norm = 'quad';
    
    prb.w_ep = w_ep;
    prb.w_px = w_px;
    prb.cost_factor = cost_factor;
    
    prb.eps_ep = 1e-7;
    prb.eps_px = 5e-4;
    
    % convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func = @(t,x,u,s) plant.rocket6DoF.dyn_func(x,u,s,prb.c_ax,prb.c_ayz,diag(prb.JB),prb.gI,prb.rho,prb.SA,prb.rTB,prb.rcpB,prb.alphmdt,prb.betmdt);
    prb.dyn_func_linearize = @(tbar,xbar,ubar,sbar) plant.rocket6DoF.compute_linearization(xbar,ubar,sbar,prb.c_ax,prb.c_ayz,diag(prb.JB),prb.gI,prb.rho,prb.SA,prb.rTB,prb.rcpB,prb.alphmdt,prb.betmdt);

end