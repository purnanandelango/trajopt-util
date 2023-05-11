clearvars
close all
clc

scp_iters       = 10;
num_targ        = 4; 
targ_pref       = [1,2,4,3];
Nhorz           = 20;

rtarg           = [20   50   0   50; 
                   70   0    50  30];    
vtarg           = [1  0  0  0.5;
                   0  1  0  1];

cost_targ       = 40*[1 1 1 1];

Kbranch         = @(M) round(M/4);

save('scenario','num_targ','targ_pref','Nhorz','rtarg','vtarg');

% Initialization
Nbranch_prev    = 0; 
rbranch_prev    = [0;0];           
vbranch_prev    = [0;0];
traj_cost_str   = 0;

for stage = 1:num_targ-1

    N           = Nhorz - Nbranch_prev;
    ntarg       = num_targ - stage + 1;
    r1          = rbranch_prev;           
    v1          = vbranch_prev;
    rK          = rtarg(:,targ_pref(stage:end));    
    vK          = vtarg(:,targ_pref(stage:end));
    cost_bound  = cost_targ(targ_pref(stage:end))-traj_cost_str;
    Kstr        = Kbranch(N);
    
    prb = problem_data_2D(N,scp_iters,25,10,0.01,0.01,...
                          ntarg,r1,v1,rK,vK,cost_bound,Kstr);
    

    % load("solve_"+num2str(stage),'x','u','tau');
    [xbar,ubar] = misc.create_initialization(prb,1);%,x,u,tau);
    
    [xbar,ubar] = scp.run_ptr_noparam(xbar,ubar,prb,@sys_cnstr_cost);
    
    [tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,1],prb.Kfine,prb.disc);
    
    [r,   v,   T,   s,   tvec,   nrm_v,   nrm_T,...
     rbar,vbar,Tbar,sbar,tvecbar,nrm_vbar,nrm_Tbar,...
     traj_cost,traj_cost_str]                               = disassemble(x,u,tau,xbar,ubar,prb);

    rbranch_prev = rbar(:,Kstr,1);
    vbranch_prev = vbar(:,Kstr,1);
    Nbranch_prev = Nbranch_prev + Kstr - 1;
    
    save("solve_"+num2str(stage),'r','v','nrm_v','T','nrm_T','s','tvec','x','u','prb','tau',...
                           'rbar','vbar','nrm_vbar','Tbar','nrm_Tbar','sbar','tvecbar','xbar','ubar','traj_cost','traj_cost_str');
    
    if stage == 1
        figure
    end
    plot_ddto_stage_traj(rbar,r,prb);    

end

if isfield(prb,'robs')
    if prb.n == 2
        th = linspace(0,2*pi);
        for j = 1:prb.nobs
            pobs = prb.robs(:,j) + prb.aobs(j)*[cos(th);sin(th)];
            plot(pobs(1,:),pobs(2,:),'-r','LineWidth',2);
        end
    end
end
