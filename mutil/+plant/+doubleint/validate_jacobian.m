clearvars 
clc
close all

n = 3;
x = randn(2*n,1); 
x(n+1:2*n) = 0; % Speed is zero
u = randn(n,1);
s = randn()^2;
accl = randn(n,1);
c_d = randn()^2;
fun = @(z) plant.doubleint.dyn_func(z(1:2*n),z(2*n+1:3*n),z(3*n+1),n,c_d,accl);

z = [x;u;s];

J1 = misc.num_jacobian(fun,z)
[A2,B2,S2,w2] = plant.doubleint.compute_linearization(x,u,s,n,c_d,accl);
J2 = [A2,B2,S2]
 
norm(J1-J2)