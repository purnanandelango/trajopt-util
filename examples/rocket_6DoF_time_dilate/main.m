clearvars
clc

prb = problem_data(20, ...              % K
                   10, ...              % scp_iters
                   2e2, ...             % wvc
                   0.1, ...             % wtr
                   0.1 ...              % cost_factor
                   );

load('recent_solution','xbar','ubar','taubar','tvecbar');
[xbar,ubar,sbar] = misc.create_initialization(prb,1,xbar,ubar,taubar,tvecbar(end));

[xbar,ubar,pbar] = scp.run_ptr(xbar,ubar,sbar,prb,@sys_cnstr_cost);
taubar = prb.tau;
tvecbar = taubar*pbar;

% Simulate solution on fine grid
[tvec,x,u] = disc.simulate_dyn(xbar(:,1),{tvecbar,ubar},@(t,x,u) prb.dyn_func(t,x,u,1.0),[tvecbar(1),tvecbar(end)],prb.Kfine,prb.disc,prb.ode_solver);

m = x(1,:);
rI = x(2:4,:);
vI = x(5:7,:);
qBI = x(8:11,:);
omgB = x(12:14,:);

fprintf('\nFinal position error: %.3f\nFinal velocity error: %.3f\n',norm(rI(:,end)-prb.rIK),norm(vI(:,end)-prb.vIK));

save('recent_solution','tvec','m','rI','vI','qBI','omgB','u','x', ...
                       'xbar','ubar','taubar','tvecbar', ...
                       'prb');

% plot_solution;