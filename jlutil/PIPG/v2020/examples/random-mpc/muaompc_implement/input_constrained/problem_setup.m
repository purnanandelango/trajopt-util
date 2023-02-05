%% Random mpc 
% - System module definition
% - Compilation


clear all
clc
rng default

% constansts
nx = 4;
nu = 8;
dt = 0.1;
N = 10;
T = dt*(N-1);

Nmpc = 30;
noise_std = 0.001;

% discrete time system
Ad_temp = randn(nx,nx); singval_Ad_temp = svd(Ad_temp);
Ad = Ad_temp/singval_Ad_temp(1);
Bd = randn(nx,nu);

% input constraints
u_lb = -10*ones(nu,1);
u_ub = 10*ones(nu,1);

% reference
x_ref = ones(nx*(Nmpc+N),1);
u_ref = randn(nu*(Nmpc+N-1),1);

% weighting matrices
Q = eye(nx);
R = eye(nu);
P = eye(nx);

% system module
save problem_data.mat

rmdir('cmpc','s');

% status = system('python3 problem_solve.py &'); % run command on terminal
% if status
%     error('Something went wrong.');
% end
% !python3 problem_solve.py &

%% compile it and make it available
cd cmpc/matlab
make
cd ../..
addpath ./cmpc/matlab