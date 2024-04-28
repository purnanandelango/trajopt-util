function [cnstr,cost_fun,ep_cnstr] = sys_cnstr_cost(x,u,prb,...
                                                    xbar,ubar)

    K = prb.K;

    r    = x(1:2,:);
    v    = x(3:4,:);
    T    = u(1:2,:);
    s    = u(3,:);

    % Obstacle avoidance buffer
    nu_ncvx = sdpvar(prb.nobs+1,K);
    
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
        
        % cost_fun = cost_fun + prb.cost_factor*(norm(u(1:prb.n,k)) + 2*(u(prb.n+1,k)));

        for j = 1:prb.nobs
            cnstr = [cnstr;
                     norm(xbar(1:prb.n,k)-prb.robs(:,j)) - prb.aobs(j) + dot(xbar(1:prb.n,k)-prb.robs(:,j),r(:,k)-xbar(1:prb.n,k))/norm(xbar(1:prb.n,k)-prb.robs(:,j)) + nu_ncvx(j,k) >= 0;
                     nu_ncvx(j,k) >= 0];
        end
        cnstr = [cnstr;
                norm(ubar(1:prb.n,k)) - prb.umin + dot(ubar(1:prb.n,k),T(:,k)-ubar(1:prb.n,k))/norm(ubar(1:prb.n,k)) + nu_ncvx(end,k) >= 0;
                nu_ncvx(end,k) >= 0];

    end  

    ep_cnstr = sum(nu_ncvx(:));    

    cost_fun = cost_fun + prb.cost_factor*norm(u(:)) + prb.w_ep*ep_cnstr;

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

end