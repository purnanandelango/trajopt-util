function [cnstr,cost_fun,vc_cnstr] = sys_cnstr_cost(x,u,prb,...
                                                    ~,~)
% r1    = x(1:3)
% v1    = x(4:6)
% r2    = x(7:9)
% v2    = x(10:12)
% T1    = u(1:3)
% s1    = u(4)
% T1    = u(5:7)
% s1    = u(8)

    K = prb.K;

    % Unscaled variables
    r1   = sdpvar(3,K);
    v1   = sdpvar(3,K);
    r2   = sdpvar(3,K);
    v2   = sdpvar(3,K);    
    T1   = sdpvar(3,K);
    s1   = sdpvar(1,K);
    T2   = sdpvar(3,K);
    s2   = sdpvar(1,K);    

    for k = 1:K
        r1(:,k)   = prb.Sx(1:3,1:3)          *x(1:3,k)      + prb.cx(1:3);
        v1(:,k)   = prb.Sx(4:6,4:6)          *x(4:6,k)      + prb.cx(4:6);
        r2(:,k)   = prb.Sx(7:9,7:9)          *x(7:9,k)      + prb.cx(7:9);
        v2(:,k)   = prb.Sx(10:12,10:12)      *x(10:12,k)    + prb.cx(10:12);        

        T1(:,k)   = prb.Su(1:3,1:3)      *u(1:3,k)   + prb.cu(1:3);        
        s1(k)     = prb.Su(4,4)          *u(4,k)     + prb.cu(4);  
        T2(:,k)   = prb.Su(5:7,5:7)      *u(5:7,k)   + prb.cu(5:7);        
        s2(k)     = prb.Su(8,8)          *u(8,k)     + prb.cu(8);          
    end
    
    % Boundary conditions
    cnstr = [r1(:,1)   == prb.r1;
             v1(:,1)   == prb.v1;
             r2(:,1)   == prb.r1;
             v2(:,1)   == prb.v1;
             r1(:,K)   == prb.rK(:,1);
             v1(:,K)   == prb.vK(:,1);
             r2(:,K)   == prb.rK(:,2);
             v2(:,K)   == prb.vK(:,2)];

    cost_fun = 0;

    for k = 1:K    
        
        cnstr = [cnstr;
                 norm(T1(:,k)) <= prb.umax;                                                                     % Thrust magnitude upper bound
                 norm(v1(:,k)) <= prb.vmax;                                                                     % Velocity magnitude upper bound
                 norm(r1(:,k),'inf') <= prb.rmax; 
                 prb.smin <= s1(k) <= prb.smax;                                                                 % Lower and upper bounds on dilation factor
                 norm(T2(:,k)) <= prb.umax;                                                                     % Thrust magnitude upper bound
                 norm(v2(:,k)) <= prb.vmax;                                                                     % Velocity magnitude upper bound
                 norm(r2(:,k),'inf') <= prb.rmax; 
                 prb.smin <= s2(k) <= prb.smax];

        if k <= prb.taustr
            cnstr = [cnstr;
                     r1(:,k) == r2(:,k);
                     v1(:,k) == v2(:,k)];
        end
        
        cost_fun = cost_fun + prb.cost_factor*(norm(u([1:3,5:7],k)) + 2*(u(4,k)+u(8,k)));
    
    end  

    % cost_fun = cost_fun + prb.cost_factor*norm(u(:));

    % Compute time of maneuver and constrain time step
    ToF1 = 0;
    ToF2 = 0;
    switch prb.disc
        case "ZOH"
            for k = 1:prb.K-1
                ToF1 = ToF1 + prb.dtau(k)*s1(k);
                ToF2 = ToF2 + prb.dtau(k)*s2(k);
                cnstr = [cnstr; prb.dtmin <= prb.dtau(k)*s1(k) <= prb.dtmax;
                                prb.dtmin <= prb.dtau(k)*s2(k) <= prb.dtmax]; 
            end
        case "FOH"
            for k = 1:prb.K-1
                ToF1 = ToF1 + 0.5*prb.dtau(k)*(s1(k+1)+s1(k));
                ToF2 = ToF2 + 0.5*prb.dtau(k)*(s2(k+1)+s2(k));
                cnstr = [cnstr; prb.dtmin <= 0.5*prb.dtau(k)*(s1(k+1)+s1(k)) <= prb.dtmax;
                                prb.dtmin <= 0.5*prb.dtau(k)*(s2(k+1)+s2(k)) <= prb.dtmax];
            end
    end    

    % Time of maneuver upper bound
    cnstr = [cnstr; ToF1 <= prb.ToFmax;
                    ToF2 <= prb.ToFmax];

    vc_cnstr = 0;

end