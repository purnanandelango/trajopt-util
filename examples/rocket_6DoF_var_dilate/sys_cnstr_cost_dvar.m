function [cnstr,cost_fun,vcvb_cnstr] = sys_cnstr_cost_dvar(x,u,prb,...
                                                           xbar,ubar)
% m     = x(1)
% rI    = x(2:4)
% vI    = x(5:7)
% qBI   = x(8:11) 
% omgB  = x(12:14)

    K = prb.K;
    
    % Boundary conditions
    cnstr = [x(1,1)     == prb.invSx(1,1)*prb.mwet;
             x(2:4,1)   == prb.invSx(2:4,2:4)*prb.rI1;
             x(5:7,1)   == prb.invSx(5:7,5:7)*prb.vI1;
             x(2:4,K)   == prb.invSx(2:4,2:4)*prb.rIK;
             x(5:7,K)   == prb.invSx(5:7,5:7)*prb.vIK;
             x(8:11,K)  == prb.invSx(8:11,8:11)*prb.q1;
             x(12:14,1) == prb.invSx(12:14,12:14)*prb.omgB1;
             x(12:14,K) == prb.invSx(12:14,12:14)*prb.omgBK];

    vb_Tmin = sdpvar(1,K);       
    % vb_STC  = sdpvar(1,K);

    cost_fun = 0;

    for k = 1:K    
        
        cnstr = [cnstr;
                 x(1,k) >= prb.invSx(1,1)*prb.mdry                                                                           % Vehicle mass lower bound
                 % norm(prb.Hgam*prb.Sx(2:4,2:4)*x(2:4,k)) <= prb.Sx(2,2)*x(2,k)/prb.cotgamgs;                                            % Glide slope constraint
                 % norm(prb.Hthet*prb.Sx(8:11,8:11)*x(8:11,k)) <= prb.sinthetmaxby2;                                             % Vehicle tilt angle constraint
                 % norm(prb.Sx(12:14,12:14)*x(12:14,k),inf) <= prb.omgmax;                                                         % Angular velocity magnitude upper bound
                 norm(u(1:3,k)) <= prb.invSu(1,1)*prb.Tmax;                                                                 % Thrust magnitude upper bound
                 prb.cosdelmax*norm(u(1:3,k)) <= u(1,k);                                                    % Thrust pointing constraint
                 norm(x(5:7,k)) <= prb.invSx(5,5)*prb.Vmax;                                                                 % Velocity magnitude upper bound
                 prb.invSu(4,4)*prb.smin <= u(4,k) <= prb.invSu(4,4)*prb.smax;                                                              % Lower and upper bounds on dilation factor
                 (u(1:3,k))'*(-ubar(1:3,k)/norm(ubar(1:3,k))) + prb.invSu(1,1)*prb.Tmin <= vb_Tmin(k);                      % Linearized thrust magnitude lower bound
                 vb_Tmin(k) >= 0];

        % Airspeed-triggered angle-of-attack STC
        % [h,dh] = plant.rocket6DoF.q_aoa_cnstr(xbar(5:7,k),xbar(8:11,k),prb.Vmax_STC,prb.cosaoamax,prb.STC_flag);
        % cnstr = [cnstr;
        %          h + dh*(prb.Sx*x(:,k)+prb.cx-xbar(:,k)) <= vb_STC(k);
        %          vb_STC(k) >= 0];
        
        % cost_fun = cost_fun + prb.cost_factor*(norm(u(1:3,k)));
        cost_fun = cost_fun + prb.cost_factor*x(1,K);
    
    end  

    % Compute time of maneuver and constrain time step
    ToF = 0;
    switch prb.disc
        case "ZOH"
            for k = 1:prb.K-1
                ToF = ToF + prb.dtau(k)*prb.Su(4,4)*u(4,k);
                cnstr = [cnstr; prb.dtmin <= prb.dtau(k)*prb.Su(4,4)*u(4,k) <= prb.dtmax]; 
            end
        case "FOH"
            for k = 1:prb.K-1
                ToF = ToF + 0.5*prb.dtau(k)*prb.Su(4,4)*(u(4,k+1)+u(4,k));
                cnstr = [cnstr; prb.dtmin <= 0.5*prb.dtau(k)*prb.Su(4,4)*(u(4,k+1)+u(4,k)) <= prb.dtmax];
            end
    end    

    % Time of maneuver upper bound
    cnstr = [cnstr;ToF <= prb.ToFmax];

    vcvb_cnstr = 0;

    vcvb_cnstr = vcvb_cnstr + prb.wvc*sum(vb_Tmin(:));

    % vcvb_cnstr = vcvb_cnstr + prb.wvc*sum(vb_STC(:));

    cost_fun = cost_fun + vcvb_cnstr;

end