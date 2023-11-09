function prb = problem_data(K,scp_iters,wvc,wtr,cost_factor)
    
    prb.K = K;

    prb.nx = 14;
    prb.nu = 4;
    prb.np = 0;
    
    prb.tau = grid.generate_grid(0,1,K,'uniform'); % Generate grid in [0,1]

    prb.dtau = diff(prb.tau);
    
    prb.h = (1/10)*prb.dtau;                % Step size for integration that computes FOH matrices
    prb.Kfine = 1+round(10/min(prb.dtau));   % Size of grid on which SCP solution is simulated
    
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

    prb.thetmax = 90*pi/180;    prb.sinthetmaxby2 = sin(prb.thetmax/2);     % Vehicle tilt
    prb.gamgs   = 75*pi/180;    prb.cotgamgs = cot(prb.gamgs);              % Glide-slope
    
    prb.omgmax  = 28.6*pi/180;                                              % Angular velocity
    prb.delmax  = 20*pi/180;    prb.cosdelmax = cos(prb.delmax);            % Gimbal
    
    prb.Hgam    = [0,1,0;0,0,1];
    prb.Hthet   = [0,1,0,0;0,0,1,0];
    
    prb.Tmin    = 1.5;
    prb.Tmax    = 6.5;
    
    prb.Vmax     = 3;
    
    prb.Vmax_STC    =  2;
    prb.cosaoamax   = cosd( 5 ); 
    prb.STC_flag    = "v1";

    prb.mdry    = 1;
    prb.mwet    = 2;

    prb.smin    = 01;
    prb.smax    = 20;
    prb.dtmin   = 0.01;
    prb.dtmax   = 3;
    prb.ToFmax  = 15;
   
    prb.snom    = [1,20];
    prb.ToFguess= 05;    
    
    % Boundary conditions

    prb.rI1     = [5.5;4.5;0];
    prb.vI1     = [-0.5;-2.5;0];
    % prb.vI1    = [-2.5;-0.5;0];
    
    prb.rIK     = zeros(3,1);
    prb.vIK     = -0.1*[1;0;0];
    prb.omgB1   = zeros(3,1);
    prb.omgBK   = prb.omgB1;
    prb.q1      = [0;0;0;1];

    % Generate straight-line initialization (unscaled)
    prb.x1      = [prb.mwet;prb.rI1;prb.vI1;prb.q1;prb.omgB1];
    prb.xK      = [prb.mdry;prb.rIK;prb.vIK;prb.q1;prb.omgBK];    
    prb.u1      = [-prb.mwet*prb.gI;prb.ToFguess];
    prb.uK      = [-prb.mdry*prb.gI;prb.ToFguess];

    % Scaling parameters
    xmin = 0*[prb.mdry; -10*ones(3,1); -2*ones(3,1); -ones(4,1); -prb.omgmax*ones(3,1)];
    xmax =   [prb.mwet;  10*ones(3,1);  2*ones(3,1);  ones(4,1);  prb.omgmax*ones(3,1)];
    
    umin = 0*[prb.Tmin*ones(3,1); prb.snom(1)];
    umax =   [prb.Tmax*ones(3,1); prb.snom(2)];

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[0,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    % SCP parameters

    prb.disc = "FOH";
    prb.foh_type = "v3";
    prb.ode_solver = {'ode45',odeset('RelTol',1e-5,'AbsTol',1e-7)};
    prb.scp_iters = scp_iters; % Maximum SCP iterations
    

    prb.solver_settings = sdpsettings('solver','gurobi','verbose',false);
    % prb.solver_settings = sdpsettings('solver','ecos','verbose',false,'ecos.AbsTol',1e-8,'ecos.RelTol',1e-8,'ecos.FeasTol',1e-8);
    % prb.solver_settings = sdpsettings('solver','ipopt','verbose',false);

    % prb.tr_norm = 2;
    % prb.tr_norm = inf;
    prb.tr_norm = 'quad';
    
    prb.wvc = wvc;
    prb.wtr = wtr;
    prb.cost_factor = cost_factor;
    
    prb.epsvc = 1e-8;
    prb.epstr = 1e-4;

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