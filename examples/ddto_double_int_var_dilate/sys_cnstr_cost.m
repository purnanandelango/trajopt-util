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

    cnstr = [];
    cost_fun = 0;

    for j = 1:ntarg

        idx_r = (j-1)*(2*n)+1:(j-1)*(2*n)+n;
        idx_v = (j-1)*(2*n)+n+1:j*(2*n);
        idx_T = (j-1)*(n+1)+1:(j-1)*(n+1)+n;
        idx_s = j*(n+1);

        % Define unscaled states and control inputs
        for k = 1:K
            r(:,k,j)   = prb.Sx(idx_r,idx_r)          *x(idx_r,k)      + prb.cx(idx_r);
            v(:,k,j)   = prb.Sx(idx_v,idx_v)          *x(idx_v,k)      + prb.cx(idx_v);
    
            T(:,k,j)   = prb.Su(idx_T,idx_T)          *u(idx_T,k)      + prb.cu(idx_T);        
            s(k,j)     = prb.Su(idx_s,idx_s)          *u(idx_s,k)      + prb.cu(idx_s);  
        end

        % Boundary conditions
        cnstr = [cnstr;
                 r(:,1,j)   == prb.r1;
                 v(:,1,j)   == prb.v1;
                 r(:,K,j)   == prb.rK(:,j);
                 v(:,K,j)   == prb.vK(:,j)];       

        for k = 1:K    
            
            cnstr = [cnstr;
                     norm(T(:,k,j)) <= prb.umax;                                                                    % Thrust magnitude upper bound
                     norm(v(:,k,j)) <= prb.vmax;                                                                    % Velocity magnitude upper bound
                     -prb.rmax <= r(:,k,j) <= prb.rmax;                                                             % Bounds on position 
                     prb.smin <= s(k,j) <= prb.smax];                                                               % Lower and upper bounds on dilation factor
            
            cost_fun = cost_fun + prb.cost_factor*(norm(u(idx_T,k)) + 2*u(idx_s,k));

            % Deferrability
            if j > 1 && k <= prb.Kstr
                cnstr = [cnstr;
                         r(:,k,1) == r(:,k,j);
                         v(:,k,1) == v(:,k,j)];
            end            
        
        end        
    
    end  

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