clearvars
clc

prb = problem_data_2D(15,15,10,0.01,0.01);

load('recent_solution','x','u','tau');
[xbar,ubar] = misc.create_initialization(prb,1,x,u,tau);

[xbar,ubar] = scp.run_ptr_noparam(xbar,ubar,prb,@sys_cnstr_cost);

% Simulate solution on fine grid in [0,1]
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,1],prb.Kfine,prb.disc);

r     = zeros(prb.n,prb.Kfine,prb.ntarg);
v     = zeros(prb.n,prb.Kfine,prb.ntarg);
T     = zeros(prb.n,prb.Kfine,prb.ntarg);
nrm_v = zeros(prb.Kfine,prb.ntarg);
nrm_T = zeros(prb.Kfine,prb.ntarg);
s     = zeros(prb.Kfine,prb.ntarg);
tvec  = zeros(prb.Kfine,prb.ntarg);

rbar     = zeros(prb.n,prb.K,prb.ntarg);
vbar     = zeros(prb.n,prb.K,prb.ntarg);
Tbar     = zeros(prb.n,prb.K,prb.ntarg);
nrm_vbar = zeros(prb.K,prb.ntarg);
nrm_Tbar = zeros(prb.K,prb.ntarg);
sbar     = zeros(prb.K,prb.ntarg);
tvecbar  = zeros(prb.K,prb.ntarg);
traj_cost = zeros(1,prb.ntarg);

for j = 1:prb.ntarg
    for k = 1:prb.Kfine
        r(:,k,j)   = x(prb.idx_r(:,j),k);
        v(:,k,j)   = x(prb.idx_v(:,j),k);
        T(:,k,j)   = u(prb.idx_T(:,j),k);
        nrm_v(k,j) = norm(v(:,k,j));        
        nrm_T(k,j) = norm(T(:,k,j));
        s(k,j)     = u(prb.idx_s(j),k);
    end
    tvec(:,j) = prb.time_grid(tau,x,s(:,j));

    for k = 1:prb.K
        rbar(:,k,j)   = xbar(prb.idx_r(:,j),k);
        vbar(:,k,j)   = xbar(prb.idx_v(:,j),k);
        Tbar(:,k,j)   = ubar(prb.idx_T(:,j),k);
        nrm_vbar(k,j) = norm(vbar(:,k,j));        
        nrm_Tbar(k,j) = norm(Tbar(:,k,j));
        sbar(k,j)     = ubar(prb.idx_s(j),k);
    end
    tvecbar(:,j) = prb.time_grid(prb.tau,xbar,sbar(:,j));

    % Trajectory cost    
    switch prb.subopt_type
        case 'sum_stage_cost'
            for k = 1:prb.K
                traj_cost(j) = traj_cost(j) + prb.cost_term(Tbar(:,k,j));
            end
        case 'sum_quad_u'
            Tj = Tbar(:,:,j);                                              
            traj_cost(j) = norm(Tj(:));
    end    

    fprintf('\nTarget %d - Final position error: %.3f\n         - Final velocity error: %.3f\n         -      Trajectory cost: %.3f\n',j,norm(r(:,end,j)-prb.rK(:,j)),norm(v(:,end,j)-prb.vK(:,j)),traj_cost(j));    

end

save('recent_solution','r','v','nrm_v','T','nrm_T','s','tvec','x','u','prb','tau',...
                       'rbar','vbar','nrm_vbar','Tbar','nrm_Tbar','sbar','tvecbar','xbar','ubar','traj_cost');

plot_solution;