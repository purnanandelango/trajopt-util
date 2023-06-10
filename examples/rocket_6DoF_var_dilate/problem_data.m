function prb = problem_data(K,scp_iters,wvc,wtr,cost_factor)
    
    prb.K = K;

    prb.nx = 14;
    prb.nu = 4;
    prb.np = 0;
    
    prb.tau = grid.generate_grid(0,1,K,'uniform'); % Generate grid in [0,1]

    prb.dtau = diff(prb.tau);
    
    prb.h = (1/10)*min(prb.dtau);            % Step size for integration that computes FOH matrices
    prb.Kfine = 1+round(2/min(prb.dtau));   % Size of grid on which SCP solution is simulated
    
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
    
    prb.Tmin = 1;
    prb.Tmax = 6.5;
    
    prb.Vmax = 3;
    
    prb.mdry = 1;
    prb.mwet = 2;

    prb.smin = 0.01;
    prb.smax = 10;
    prb.dtmin = 0.01;
    prb.dtmax = 2;
    prb.ToFmax = 20;
    
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
    prb.u1 = [-prb.mwet*prb.gI;5/K];
    prb.uK = [-prb.mdry*prb.gI;5/K];

    % Scaling parameters
    xmin1 = [1;   0;  0; -1; -2*ones(3,1); -ones(4,1); -prb.omgmax*ones(3,1)];
    xmax1 = [2;   7;  7;  1;  2*ones(3,1);  ones(4,1);  prb.omgmax*ones(3,1)];

    xminK = [1;   0;  0; -1; -2*ones(3,1); -ones(4,1); -prb.omgmax*ones(3,1)];
    xmaxK = [2;   7;  7;  1;  2*ones(3,1);  ones(4,1);  prb.omgmax*ones(3,1)];
    
    umin1 = [1*ones(3,1); 1];
    umax1 = [6*ones(3,1); 15];

    uminK = [1*ones(3,1); 5];
    umaxK = [6*ones(3,1); 10];    

    xmin  = grid.ends2interp(xmin1,xminK,prb.tau,'poly',1);
    xmax  = grid.ends2interp(xmax1,xmaxK,prb.tau,'poly',1);
    xbnd  = linalg.matcat(xmin,xmax,3);

    umin  = grid.ends2interp(umin1,uminK,prb.tau,'poly',1);
    umax  = grid.ends2interp(umax1,umaxK,prb.tau,'poly',1);
    ubnd  = linalg.matcat(umin,umax,3);

    [Sz,cz] = misc.generate_varscaling({xbnd,ubnd},[-1,1]);

    prb.Sx = Sz(:,1);
    prb.Su = Sz(:,2);
    prb.cx = cz(:,1);
    prb.cu = cz(:,2);
    prb.invSx = cellfun(@inv,prb.Sx,'UniformOutput',false);
    prb.invSu = cellfun(@inv,prb.Su,'UniformOutput',false);

    % SCP parameters

    prb.disc = "FOH";
    prb.foh_type = "v3";
    prb.scp_iters = scp_iters; % Maximum SCP iterations

    % prb.solver_settings = sdpsettings('solver','ecos','verbose',false,'ecos.AbsTol',1e-8,'ecos.RelTol',1e-8,'ecos.FeasTol',1e-9);
    prb.solver_settings = sdpsettings('solver','ipopt','verbose',false);

    prb.tr_norm = 2;
    % prb.tr_norm = inf;
    % prb.tr_norm = 'quad';
    
    prb.wvc = wvc;
    prb.wtr = wtr;
    prb.cost_factor = cost_factor;
    
    prb.epsvc = 1e-8;
    prb.epstr = 1e-7;

    % Takes in unscaled data
    prb.time_of_maneuver = @(x,u) disc.time_of_maneuver(prb.disc,prb.tau,u(4,:));    
    prb.time_grid = @(tau,x,u) disc.time_grid(prb.disc,tau,u(4,:));    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func = @(t,x,u) plant.rocket6DoF.dyn_func(x,u(1:3),u(4),prb.c_ax,prb.c_ayz,diag(prb.JB),prb.gI,prb.rho,prb.SA,prb.rTB,prb.rcpB,prb.alphmdt,prb.betmdt);
    prb.dyn_func_linearize = @(tbar,xbar,ubar) evaluate_linearization(xbar,ubar,prb.c_ax,prb.c_ayz,diag(prb.JB),prb.gI,prb.rho,prb.SA,prb.rTB,prb.rcpB,prb.alphmdt,prb.betmdt);

end

function [A,B,w] = evaluate_linearization(x,u,c_ax,c_ayz,JB,gI,rho,SA,rTB,rcpB,alphmdt,betmdt)
    [A,B,S,w] = plant.rocket6DoF.compute_linearization(x,u(1:3),u(4),c_ax,c_ayz,JB,gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
    B = [B,S];
end