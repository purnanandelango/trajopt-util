function [cnstr,cost_fun,vc_cnstr] = sys_cnstr_cost(x,u,prb,...
                                                    ~,~)
% r    = x(1:3)
% v    = x(4:6)
% T    = u(1:3)
% s    = u(4)

    K = prb.K;

    % Unscaled variables
    r   = sdpvar(3,K);
    v   = sdpvar(3,K);
    T   = sdpvar(3,K);
    s   = sdpvar(1,K);

    for k = 1:K
        r(:,k)   = prb.Sx(1:3,1:3)      *x(1:3,k)   + prb.cx(1:3);
        v(:,k)   = prb.Sx(4:6,4:6)      *x(4:6,k)   + prb.cx(4:6);

        T(:,k)   = prb.Su(1:3,1:3)      *u(1:3,k)   + prb.cu(1:3);        
        s(k)     = prb.Su(4,4)         *u(4,k)     + prb.cu(4);        
    end
    
    % Boundary conditions
    cnstr = [r(:,1)   == prb.r1;
             v(:,1)   == prb.v1;
             r(:,K)   == prb.rK;
             v(:,K)   == prb.vK];

    cost_fun = 0;

    for k = 1:K    
        
        cnstr = [cnstr;
                 norm(T(:,k)) <= prb.umax;                                                                     % Thrust magnitude upper bound
                 norm(v(:,k)) <= prb.vmax;                                                                     % Velocity magnitude upper bound
                 norm(r(:,k),'inf') <= prb.rmax; 
                 prb.smin <= s(k) <= prb.smax];                                                                % Lower and upper bounds on dilation factor
        
        cost_fun = cost_fun + prb.cost_factor*(norm(u(1:3,k)) + 2*(u(4,k)));
    
    end  

    % cost_fun = cost_fun + prb.cost_factor*norm(u(:));

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