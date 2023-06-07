clearvars 
clc
close all


n = 3;
t = randn()^2;
x = randn(2*n,1); 
% x(n+1:2*n) = 0; % Speed is zero
u = randn(n,1);
T = t + randn()^2;

fun = @(z) plant.doubleint_ltv.dyn_func(t,z(1:2*n),z(2*n+1:3*n),T,n);

z = [x;u];

J1 = misc.num_jacobian(fun,z)
[A2,B2,w2] = plant.doubleint_ltv.compute_linearization(t,x,u,T,n);
J2 = [A2,B2]
 
norm(J1-J2)