clearvars
clc

N = 15;

prb = problem_data_2D(N,10,25,10,0.01,0.01);

load('recent_solution','x','u','tau');
[xbar,ubar] = misc.create_initialization(prb,1,x,u,tau);

[xbar,ubar] = scp.run_ptr_noparam(xbar,ubar,prb,@sys_cnstr_cost);

% Simulate solution on fine grid in [0,1]
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,1],prb.Kfine,prb.disc);

[r,   v,   T,   s,   tvec,   nrm_v,   nrm_T,...
 rbar,vbar,Tbar,sbar,tvecbar,nrm_vbar,nrm_Tbar,...
 traj_cost] = disassemble(x,u,tau,xbar,ubar,prb);    

save('recent_solution','r','v','nrm_v','T','nrm_T','s','tvec','x','u','prb','tau',...
                       'rbar','vbar','nrm_vbar','Tbar','nrm_Tbar','sbar','tvecbar','xbar','ubar','traj_cost');

plot_solution;