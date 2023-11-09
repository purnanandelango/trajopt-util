function [cnstr,cost_fun,vcvb_cnstr] = sys_cnstr_cost(x,u,prb,...
                                                      xbar,ubar)

    K = prb.K;

    m    = x(1,:);
    rI   = x(2:4,:);
    vI   = x(5:7,:);
    qBI  = x(8:11,:);
    omgB = x(12:14,:);
    TB   = u(1:3,:);
    s    = u(4,:);
    
    % Boundary conditions
    cnstr = [m(1)      == prb.mwet;
             rI(:,1)   == prb.rI1;
             vI(:,1)   == prb.vI1;
             rI(:,K)   == prb.rIK;
             vI(:,K)   == prb.vIK;
             qBI(:,K)  == prb.q1;
             omgB(:,1) == prb.omgB1;
             omgB(:,K) == prb.omgBK];

    vb_Tmin = sdpvar(1,K);       

    cost_fun = 0;

    for k = 1:K    
        
        cnstr = [cnstr;
                 m(k) >= prb.mdry                                                                           % Vehicle mass lower bound
                 norm(prb.Hgam*rI(:,k)) <= rI(1,k)/prb.cotgamgs;                                            % Glide slope constraint
                 norm(prb.Hthet*qBI(:,k)) <= prb.sinthetmaxby2;                                             % Vehicle tilt angle constraint
                 norm(omgB(:,k),inf) <= prb.omgmax;                                                         % Angular velocity magnitude upper bound
                 norm(TB(:,k)) <= prb.Tmax;                                                                 % Thrust magnitude upper bound
                 prb.cosdelmax*norm(TB(:,k)) <= TB(1,k);                                                    % Thrust pointing constraint
                 norm(vI(:,k)) <= prb.Vmax;                                                                 % Velocity magnitude upper bound
                 prb.smin <= s(k) <= prb.smax;                                                              % Lower and upper bounds on dilation factor
                 (TB(:,k))'*(-ubar(1:3,k)/norm(ubar(1:3,k))) + prb.Tmin <= vb_Tmin(k);                      % Linearized thrust magnitude lower bound
                 vb_Tmin(k) >= 0];
    
    end  
    cost_fun = cost_fun + prb.cost_factor*x(1,K)*prb.invSx(1,1);    

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

    vcvb_cnstr = 0;

    vcvb_cnstr = vcvb_cnstr + prb.wvc*sum(vb_Tmin(:));

    % Airspeed-triggered angle-of-attack STC
    if prb.STC_flag == "v1" || prb.STC_flag == "v2"
        vb_STC  = sdpvar(1,K);
        vcvb_cnstr = vcvb_cnstr + prb.wvc*sum(vb_STC(:));
        for k = 1:K
            [h,dh] = plant.rocket6DoF.q_aoa_cnstr(xbar(5:7,k),xbar(8:11,k),prb.Vmax_STC,prb.cosaoamax,prb.STC_flag);
            cnstr = [cnstr;
                     h + dh*(x(:,k)-xbar(:,k)) <= vb_STC(k);
                     vb_STC(k) >= 0];    
        end
    end

    cost_fun = cost_fun + vcvb_cnstr;    

end