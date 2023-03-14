clearvars
clc

prb = problem_data(20,10,5e1,0.4,0.3);

load('recent_solution','x','u','tau');
[xbar,ubar] = misc.create_initialization(prb,1,x,u,tau);

[xbar,ubar] = scp.run_ptr_noparam(xbar,ubar,prb,@sys_cnstr_cost);
 

% Simulate solution on fine grid

% Simulate on [0,1] grid
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,1],prb.Kfine,prb.disc);

r1 = x(1:3,:);
v1 = x(4:6,:);
r2 = x(7:9,:);
v2 = x(10:12,:);
T1 = u(1:3,:);
s1 = u(4,:);
T2 = u(5:7,:);
s2 = u(8,:);

tvec1 = prb.time_grid(tau,x,s1);
tvec2 = prb.time_grid(tau,x,s2);

fprintf('\n1 - Final position error: %.3f\n1 - Final velocity error: %.3f\n',norm(r1(:,end)-prb.rK(:,1)),norm(v1(:,end)-prb.vK(:,1)));
fprintf('\n2 - Final position error: %.3f\n2 - Final velocity error: %.3f\n',norm(r2(:,end)-prb.rK(:,2)),norm(v2(:,end)-prb.vK(:,2)));

save('recent_solution','r1','v1','tvec1','r2','v2','tvec2','tau','T1','T2','s1','s2','u','x','prb',...
                       'xbar','ubar','tvecbar');

% plot_solution;