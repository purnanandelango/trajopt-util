function [cnstr,cost_fun,vc_cnstr] = sys_cnstr_cost(x,u,prb,...
                                                    ~,ubar)
% m     = x(1)
% rI    = x(2:4)
% vI    = x(5:7)
% qBI   = x(8:11) 
% omgB  = x(12:14)

    K = prb.K;

    % Unscaled variables
    m    = sdpvar(1,K);
    rI   = sdpvar(3,K);
    vI   = sdpvar(3,K);
    qBI  = sdpvar(4,K);
    omgB = sdpvar(3,K);
    TB   = sdpvar(3,K);
    s    = sdpvar(1,K);

    for k = 1:K
        m(k)      = prb.Sx{k}(1,1)          *x(1,k)     + prb.cx{k}(1);
        rI(:,k)   = prb.Sx{k}(2:4,2:4)      *x(2:4,k)   + prb.cx{k}(2:4);
        vI(:,k)   = prb.Sx{k}(5:7,5:7)      *x(5:7,k)   + prb.cx{k}(5:7);
        qBI(:,k)  = prb.Sx{k}(8:11,8:11)    *x(8:11,k)  + prb.cx{k}(8:11);
        omgB(:,k) = prb.Sx{k}(12:14,12:14)  *x(12:14,k) + prb.cx{k}(12:14);

        TB(:,k)   = prb.Su{k}(1:3,1:3)      *u(1:3,k)   + prb.cu{k}(1:3);        
        s(k)      = prb.Su{k}(4,4)          *u(4,k)     + prb.cu{k}(4);        
    end
    
    % Boundary conditions
    cnstr = [m(1)      == prb.mwet;
             rI(:,1)   == prb.rI1;
             vI(:,1)   == prb.vI1;
             rI(:,K)   == prb.rIK;
             vI(:,K)   == prb.vIK;
             qBI(:,K)  == prb.q1;
             omgB(:,1) == prb.omgB1;
             omgB(:,K) == prb.omgBK];

    cost_fun = 0;

    for k = 1:K    
        
        cnstr = [cnstr;
                 m(k) >= prb.mdry                                                                               % Vehicle mass lower bound
                 norm(prb.Hgam*rI(:,k)) <= rI(1,k)/prb.cotgamgs;                                                % Glide slope constraint
                 norm(prb.Hthet*qBI(:,k)) <= prb.sinthetmaxby2;                                                % Vehicle tilt angle constraint
                 norm(omgB(:,k),inf) <= prb.omgmax;                                                             % Angular velocity magnitude upper bound
                 norm(TB(:,k)) <= prb.Tmax;                                                                     % Thrust magnitude upper bound
                 prb.cosdelmax*norm(TB(:,k)) <= TB(1,k);                                                        % Thrust pointing constraint
                 norm(vI(:,k)) <= prb.Vmax;                                                                     % Velocity magnitude upper bound
                 prb.smin <= s(k) <= prb.smax;                                                                  % Lower and upper bounds on dilation factor
                 (TB(:,k))'*(-ubar(1:3,k)/norm(ubar(1:3,k))) + prb.Tmin <= 0];                                  % Linearized thrust magnitude lower bound
        
        cost_fun = cost_fun + prb.cost_factor*(norm(u(1:3,k))+2*abs(u(4,k)));
    
    end  

    % Compute time of maneuver and constrain time step
    ToF = 0;
    switch prb.disc
        case "ZOH"
            for k = 1:prb.K-1
                ToF = ToF + prb.dtau(k)*s(k);
                cnstr = [cnstr; prb.dtmin <= prb.dtau(k)*s(k) <= prb.dtmax]; 
            end
        case "FOH"
            for k = 1:prb.K-1
                ToF = ToF + 0.5*prb.dtau(k)*(s(k+1)+s(k));
                cnstr = [cnstr; prb.dtmin <= 0.5*prb.dtau(k)*(s(k+1)+s(k)) <= prb.dtmax];
            end
    end    

    % Time of maneuver upper bound
    cnstr = [cnstr;ToF <= prb.ToFmax];

    vc_cnstr = 0;

end