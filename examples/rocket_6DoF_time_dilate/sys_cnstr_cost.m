function [cnstr,cost_fun,ep_cnstr] = sys_cnstr_cost(x,u,~,prb,...
                                                      ~,ubar,~)

    K = prb.K;

    m    = x(1,:);
    rI   = x(2:4,:);
    vI   = x(5:7,:);
    qBI  = x(8:11,:);
    omgB = x(12:14,:);
    TB   = u(1:3,:);

    % Boundary conditions
    cnstr = [m(1) == prb.mwet;
             rI(:,1) == prb.rI1;
             vI(:,1) == prb.vI1;
             rI(:,K) == prb.rIK;
             vI(:,K) == prb.vIK;
             qBI(:,1) == prb.q1;
             qBI(:,K) == prb.q1;
             omgB(:,1) == prb.omgB1;
             omgB(:,K) == prb.omgBK];

    cost_fun = 0;

    ep_Tmin = sdpvar(1,K);

    for k = 1:K 
        cnstr = [cnstr;
                 m(k) >= prb.mdry;                                                      % Vehicle mass lower bound
                 norm(prb.Hgam*rI(:,k)) <= rI(1,k)/prb.cotgamgs;                        % Glide slope constraint
                 norm(prb.Hthet*qBI(:,k)) <= prb.sinthetmaxby2;                         % Vehicle tilt angle constraint
                 norm(omgB(:,k)) <= prb.omgmax;                                         % Angular velocity magnitude upper bound  
                 norm(TB(:,k)) <= prb.Tmax;                                             % Thrust magnitude upper bound      
                 prb.cosdelmax*norm(TB(:,k)) <= TB(1,k);                                % Thrust pointing constraint 
                 norm(vI(:,k)) <= prb.Vmax;                                             % Speed upper bound
                 TB(:,k)'*(-ubar(:,k)/norm(ubar(:,k))) + prb.Tmin <= ep_Tmin(k);        % Linearized thrust magnitude lower bound
                 ep_Tmin(k) >= 0];
    
    end  

    cost_fun = cost_fun + prb.cost_factor*x(1,K)*prb.invSx(1,1);

    ep_cnstr = prb.w_ep*sum(ep_Tmin(:)); % Slack for exactly penalized constraint

    cost_fun = cost_fun + ep_cnstr;

end