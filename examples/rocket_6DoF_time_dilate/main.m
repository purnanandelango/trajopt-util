clearvars
clc

prb = problem_data(30,10,2e2,0.1,0.1);

load('recent_solution','x','u','tvec');
[xbar,ubar,sbar] = misc.create_initialization(prb,1,x,u,tvec(end));

[xbar,ubar,pbar] = scp.run_ptr(xbar,ubar,sbar,prb,@sys_cnstr_cost);

% Simulate solution on fine grid
tvec = prb.tau*pbar;
[t,x,u] = disc.simulate_dyn(xbar(:,1),{tvec,ubar},@(t,x,u) prb.dyn_func(t,x,u,1.0),[tvec(1),tvec(end)],prb.Kfine,prb.disc);

m = x(1,:);
rI = x(2:4,:);
vI = x(5:7,:);
qBI = x(8:11,:);
omgB = x(12:14,:);
tvec = t;

fprintf('\nFinal position error: %.3f\nFinal velocity error: %.3f\n',norm(rI(:,end)-prb.rIK),norm(vI(:,end)-prb.vIK));

save('recent_solution','m','rI','vI','qBI','omgB','tvec','u','x','prb');

plot_solution;