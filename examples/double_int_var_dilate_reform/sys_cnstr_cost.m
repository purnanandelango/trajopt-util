function [cnstr,cost_fun,vc_cnstr] = sys_cnstr_cost(z,u,prb,...
                                                    ~,~)
% r    = z(1:n)
% v    = z(n+1:2*n)
% bet  = z(2*n+1:2*n+m)
% T    = u(1:n)
% s    = u(n+1)

    K = prb.K;

    % Unscaled variables
    r   = sdpvar(prb.n,K);
    v   = sdpvar(prb.n,K);
    T   = sdpvar(prb.n,K);
    s   = sdpvar(1,K);

    for k = 1:K
        r(:,k)   = prb.Sx(1:prb.n,1:prb.n)                  *z(1:prb.n,k)          + prb.cx(1:prb.n);
        v(:,k)   = prb.Sx(prb.n+1:2*prb.n,prb.n+1:2*prb.n)  *z(prb.n+1:2*prb.n,k)  + prb.cx(prb.n+1:2*prb.n);

        T(:,k)   = prb.Su(1:prb.n,1:prb.n)                  *u(1:prb.n,k)          + prb.cu(1:prb.n);        
        s(k)     = prb.Su(prb.n+1,prb.n+1)                  *u(prb.n+1,k)          + prb.cu(prb.n+1); 
    end
    bet = z(2*prb.n+1:2*prb.n+prb.m,:);
    
    % Boundary conditions
    cnstr = [r(:,1)   == prb.r1;
             v(:,1)   == prb.v1;
             r(:,K)   == prb.rK;
             v(:,K)   == prb.vK;
             bet(:,1) == zeros(prb.m,1)];

    cost_fun = 0;

    for k = 1:K    
        
        cnstr = [cnstr;
                 norm(T(:,k)) <= prb.umax;                                                                     % Thrust magnitude upper bound
                 % norm(v(:,k)) <= prb.vmax;                                                                     % Velocity magnitude upper bound
                 norm(r(:,k),'inf') <= prb.rmax; 
                 prb.smin <= s(k) <= prb.smax];                                                                % Lower and upper bounds on dilation factor
        
        % cost_fun = cost_fun + prb.cost_factor*(norm(u(1:prb.n,k)) + 2*(u(prb.n+1,k)));

        if k < K
            cnstr = [cnstr;
                     bet(:,k+1) == bet(:,k)];
        end

    end  

    cost_fun = cost_fun + prb.cost_factor*norm(u(:));

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