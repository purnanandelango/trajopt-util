function [cnstr,cost_fun] = sys_cnstr_cost(x,u,~,prb,...
                                           ~,ubar,~)
% m     = x(1)
% rI    = x(2:4)
% vI    = x(5:7)
% qBI   = x(8:11) 
% omgB  = x(12:14)
% TB    = u

    K = prb.K;

    % Unscaled variables
    m    = prb.Sx(1,1)*x(1,:)               + prb.cx(1);
    rI   = prb.Sx(2:4,2:4)*x(2:4,:)         + repmat(prb.cx(2:4),[1,K]);
    vI   = prb.Sx(5:7,5:7)*x(5:7,:)         + repmat(prb.cx(5:7),[1,K]);
    qBI  = prb.Sx(8:11,8:11)*x(8:11,:)      + repmat(prb.cx(8:11),[1,K]);
    omgB = prb.Sx(12:14,12:14)*x(12:14,:)   + repmat(prb.cx(12:14),[1,K]);
    TB   = prb.Su*u                         + repmat(prb.cu,[1,K]);

    % Boundary conditions
    cnstr = [m(1) == prb.mwet;
             rI(:,1) == prb.rI1;
             vI(:,1) == prb.vI1;
             rI(:,K) == prb.rIK;
             vI(:,K) == prb.vIK;
             qBI(:,1) == prb.q1;
             omgB(:,1) == prb.omgB1;
             omgB(:,K) == prb.omgBK];

    cost_fun = 0;

    for k = 1:K 
        cnstr = [cnstr;
                 m(k) >= prb.mdry;                                          % Vehicle mass lower bound
                 norm(prb.Hgam*rI(:,k)) <= rI(1,k)/prb.cotgamgs;            % Glide slope constraint
                 2*norm(prb.Hthet*qBI(:,k)) <= 1-prb.costhetmax;            % Vehicle tilt angle constraint
                 norm(omgB(:,k)) <= prb.omgmax;                             % Angular velocity magnitude upper bound  
                 norm(TB(:,k)) <= prb.Tmax;                                 % Thrust magnitude upper bound      
                 prb.cosdelmax*norm(TB(:,k)) <= TB(1,k);                    % Thrust pointing constraint 
                 norm(vI(:,k)) <= prb.Vmax;                                 % Speed upper bound
                 TB(:,k)'*(-ubar(:,k)/norm(ubar(:,k))) + prb.Tmin <= 0;     % Linearized thrust magnitude lower bound
                 ];
        
        cost_fun = cost_fun + prb.cost_factor*norm(u(:,k));
    
    end  

end