%% Oscillating masses
% - System module definition
% - Compilation

clear all
clc
rng default

% constansts
nx = 20;
nu = 10;
dt = 0.1;
N = 10;
T = dt*(N-1);

% mpc
noise_std = 0.001;
Nmpc = 30;

% discrete time system
Ad = [zeros(nx/2) eye(nx/2);
      (diag(ones(nx/2-1,1),-1)+diag(-2*ones(nx/2,1))+diag(ones(nx/2-1,1),1)) zeros(nx/2)];
Bd = [zeros(nx/2);
     eye(nx/2)];

% state and input constraints
% u_lb = -2*ones(nu*N,1);
% u_ub = 2*ones(nu*N,1);
% e_ub = repmat([2*ones(nx/2,1);0.8*ones(nx/2,1)],[nx*N,1]);
% e_lb = -e_ub; 
% Kx = eye(nx);

u_lb = -2*ones(nu,1);
u_ub = 2*ones(nu,1);
e_ub = [2*ones(nx/2,1);0.8*ones(nx/2,1)];
e_lb = -e_ub;
Kx = eye(nx);

% terminal state constraints
f_lb = e_lb(1:nx);
f_ub = e_ub(1:nx);
F = eye(nx);

% reference
x_ref = repmat([ones(nx/2,1);zeros(nx/2,1)],[Nmpc+N,1]);
% x_ref = zeros(nu*(Nmpc+N),1);
u_ref = zeros(nu*(Nmpc+N-1),1);

% weighting matrices
Q = eye(nx);
R = eye(nu);
P = eye(nx);

% system module
save problem_data.mat

rmdir('cmpc','s')

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