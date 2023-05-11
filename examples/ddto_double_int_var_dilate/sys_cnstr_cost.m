function [cnstr,cost_fun,vb_cnstr] = sys_cnstr_cost(x,u,prb,...
                                                    xbar,~)

    K = prb.K;
    n = prb.n;
    ntarg = prb.ntarg;

    % Containers for unscaled variables
    r   = sdpvar(n,K,ntarg);
    v   = sdpvar(n,K,ntarg);
    T   = sdpvar(n,K,ntarg);
    s   = sdpvar(K,ntarg);

    % Obstacle avoidance buffer
    nu_ncvx = sdpvar(prb.nobs,ntarg,K);    

    rbar = zeros(n,K);
    vbar = zeros(n,K);

    % Define unscaled states and control inputs
    for j = 1:ntarg

        for k = 1:K
            r(:,k,j)   = prb.Sx(prb.idx_r(:,j),prb.idx_r(:,j)) *x(prb.idx_r(:,j),k)  + prb.cx(prb.idx_r(:,j));
            v(:,k,j)   = prb.Sx(prb.idx_v(:,j),prb.idx_v(:,j)) *x(prb.idx_v(:,j),k)  + prb.cx(prb.idx_v(:,j));
    
            T(:,k,j)   = prb.Su(prb.idx_T(:,j),prb.idx_T(:,j)) *u(prb.idx_T(:,j),k)  + prb.cu(prb.idx_T(:,j));        
            s(k,j)     = prb.Su(prb.idx_s(j),prb.idx_s(j))     *u(prb.idx_s(j),k)    + prb.cu(prb.idx_s(j));

            rbar(:,k,j) = xbar(prb.idx_r(:,j),k);
            vbar(:,k,j) = xbar(prb.idx_v(:,j),k);            
        end

    end

    cnstr = [];

    for j = 1:ntarg

        % Boundary conditions
        cnstr = [cnstr;
                 r(:,1,j)   == prb.r1;
                 v(:,1,j)   == prb.v1;
                 r(:,K,j)   == prb.rK(:,j);
                 v(:,K,j)   == prb.vK(:,j)];       

        % Constraints
        for k = 1:K   
            
            cnstr = [cnstr;
                     norm(T(:,k,j)) <= prb.umax;                                                                    % Thrust magnitude upper bound
                     norm(v(:,k,j)) <= prb.vmax;                                                                    % Velocity magnitude upper bound
                     -prb.rmax <= r(:,k,j) <= prb.rmax;                                                             % Bounds on position 
                     prb.smin <= s(k,j) <= prb.smax];                                                               % Lower and upper bounds on dilation factor

            % Deferrability
            if j > 1 && k <= prb.Kstr
                cnstr = [cnstr;
                         r(:,k,1) == r(:,k,j);
                         v(:,k,1) == v(:,k,j);
                         T(:,k,1) == T(:,k,j)];
            end

            for i = 1:prb.nobs
                cnstr = [cnstr;
                         norm(rbar(:,k,j)-prb.robs(:,i)) - prb.aobs(i) + dot(rbar(:,k,j)-prb.robs(:,i),r(:,k,j)-rbar(:,k,j))/norm(rbar(:,k,j)-prb.robs(:,i)) + nu_ncvx(i,j,k) >= 0;
                         nu_ncvx(i,j,k) >= 0];
            end            
        
        end
    
    end  

    % Suboptimality constraint    
    switch prb.subopt_type
        case 'sum_stage_cost'
            stage_cost = sdpvar(K,ntarg);            
            for j = 1:ntarg
                for k = 1:K
                    % stage_cost(k,j) = prb.cost_term(u(prb.idx_T(:,j),k));       % Scaled
                    stage_cost(k,j) = prb.cost_term(T(:,k,j));                  % Unscaled
                end
                cnstr = [cnstr; sum(stage_cost(:,j)) <= prb.cost_bound(j)];                
            end
        case 'sum_quad_u'
            for j = 1:ntarg
                % Tj = u(prb.idx_T(:,j),:);                                       % Scaled
                Tj = T(:,:,j);                                                  % Unscaled
                cnstr = [cnstr; norm(Tj(:)) <= prb.cost_bound(j)];
            end
    end

    % Time elapsed in each interval
    dt = sdpvar(K-1,ntarg);

    % Compute time of maneuver, constrain time step and rate of change of thrust magnitude
    for j = 1:ntarg
        for k = 1:prb.K-1
            switch prb.disc
                case "ZOH"
                    dt(k,j) = prb.dtau(k)*s(k,j); 
                case "FOH"
                    dt(k,j) = 0.5*prb.dtau(k)*(s(k+1,j)+s(k,j));
            end
            % Upper and lower bound on time elapsed in each interval
            cnstr = [cnstr; prb.dtmin <= dt(k,j) <= prb.dtmax];             
            % Upper bound on magnitude of rate of change of thrust
            cnstr = [cnstr; norm(T(:,k+1,j)-T(:,k,j)) <= dt(k,j)*prb.dTmax];            
        end
        % Time of maneuver upper bound
        cnstr = [cnstr; sum(dt(:,j)) <= prb.ToFmax];                        
    end

    % Time available to defer decision
    defer_time = sum(dt(1:prb.Kstr,1));

    vb_cnstr = sum(nu_ncvx(:));    

    cost_fun = prb.cost_factor*defer_time + prb.wvb*vb_cnstr;    

end