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
    % defer_time = 0;
    stage_cost = sdpvar(K,ntarg);

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
            
            % cost_fun = cost_fun + prb.cost_factor*prb.cost_term(u(prb.idx_T(:,j),k));
            stage_cost(k,j) = prb.cost_term(u(prb.idx_T(:,j),k));

            % Deferrability
            if j > 1 && k <= prb.Kstr
                cnstr = [cnstr;
                         r(:,k,1) == r(:,k,j);
                         v(:,k,1) == v(:,k,j)];
            end            
        
        end    

        cnstr = [cnstr; sum(stage_cost(:,j)) <= prb.cost_bound(j)];
    
    end  

    dt = sdpvar(K-1,ntarg);

    % Compute time of maneuver and constrain time step
    for j = 1:ntarg
        switch prb.disc
            case "ZOH"
                for k = 1:prb.K-1
                    dt(k,j) = prb.dtau(k)*s(k,j);
                    cnstr = [cnstr; prb.dtmin <= prb.dtau(k)*s(k,j) <= prb.dtmax]; 
                end
            case "FOH"
                for k = 1:prb.K-1
                    dt(k,j) = 0.5*prb.dtau(k)*(s(k+1,j)+s(k,j));
                    cnstr = [cnstr; prb.dtmin <= 0.5*prb.dtau(k)*(s(k+1,j)+s(k,j)) <= prb.dtmax];
                end
        end        
        % Time of maneuver upper bound
        cnstr = [cnstr; sum(dt(:,j)) <= prb.ToFmax];                        
    end
    defer_time = sum(dt(1:prb.Kstr,1));

    % cost_fun = cost_fun - prb.cost_factor*defer_time;     
    cost_fun = cost_fun + prb.cost_factor*(sum(dt(:,2))-sum(dt(:,3)));     

    vc_cnstr = 0;

end