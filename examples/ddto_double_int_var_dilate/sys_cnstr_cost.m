function [cnstr,cost_fun,vc_cnstr] = sys_cnstr_cost(x,u,prb,...
                                                    ~,~)

    K = prb.K;
    n = prb.n;
    ntarg = prb.ntarg;

    % Containers for unscaled variables
    r   = sdpvar(n,K,ntarg);
    v   = sdpvar(n,K,ntarg);
    T   = sdpvar(n,K,ntarg);
    s   = sdpvar(K,ntarg);

    for j = 1:ntarg

        % Define unscaled states and control inputs
        for k = 1:K
            r(:,k,j)   = prb.Sx(prb.idx_r(:,j),prb.idx_r(:,j)) *x(prb.idx_r(:,j),k)  + prb.cx(prb.idx_r(:,j));
            v(:,k,j)   = prb.Sx(prb.idx_v(:,j),prb.idx_v(:,j)) *x(prb.idx_v(:,j),k)  + prb.cx(prb.idx_v(:,j));
    
            T(:,k,j)   = prb.Su(prb.idx_T(:,j),prb.idx_T(:,j)) *u(prb.idx_T(:,j),k)      + prb.cu(prb.idx_T(:,j));        
            s(k,j)     = prb.Su(prb.idx_s(j),prb.idx_s(j)) *u(prb.idx_s(j),k)      + prb.cu(prb.idx_s(j));  
        end

    end

    cnstr = [];
    cost_fun = 0;
    defer_time = 0;

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
            
            cost_fun = cost_fun + prb.cost_factor*prb.cost_term(u(prb.idx_T(:,j),k));

            % Deferrability
            if j > 1 && k <= prb.Kstr
                defer_time = defer_time + u(prb.idx_s(j),k);
                cnstr = [cnstr;
                         r(:,k,1) == r(:,k,j);
                         v(:,k,1) == v(:,k,j)];
            end            
        
        end        
    
    end  

    cost_fun = cost_fun - 0.2*prb.cost_factor*defer_time; 

    % Compute time of maneuver and constrain time step
    for j = 1:ntarg
        ToFj = 0;
        switch prb.disc
            case "ZOH"
                for k = 1:prb.K-1
                    ToFj = ToFj + prb.dtau(k)*s(k,j);
                    cnstr = [cnstr; prb.dtmin <= prb.dtau(k)*s(k,j) <= prb.dtmax]; 
                end
            case "FOH"
                for k = 1:prb.K-1
                    ToFj = ToFj + 0.5*prb.dtau(k)*(s(k+1,j)+s(k,j));
                    cnstr = [cnstr; prb.dtmin <= 0.5*prb.dtau(k)*(s(k+1,j)+s(k,j)) <= prb.dtmax];
                end
        end    
        % Time of maneuver upper bound
        cnstr = [cnstr; ToFj <= prb.ToFmax];                        
    end

    vc_cnstr = 0;

end