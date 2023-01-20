clearvars
clc

%% Parameters

% Mass-depletion parameters
syms alphmdt betmdt

% Density and reference surface area
syms rho S_A

% Aerodynamics matrix
syms c_ax c_ayz
CA = diag([c_ax,c_ayz,c_ayz]);

% Inertial gravity vector
gI = sym('g_I%d',[3,1]);

% Inertia matrix
JBvec = sym('J_B%d',[3,1]);
JB = diag(JBvec);
JBinv = diag(1./JBvec);

% Location of engine gimbal pivot point
rTB = sym('r_TB%d',[3,1]);

% Location of center of pressure
rcpB = sym('r_cpB%d',[3,1]);

%% Variables

% Time dilation and mass
syms s m

% Intertial position
rI = sym('r_I%d',[3,1]);

% Inertial velocity
vI = sym('v_I%d',[3,1]);

% Quaternion (scalar-last)
q = sym('q%d',[4,1]);

% Angular velocity
omegaB = sym('omega_B%d',[3,1]);

% Control inputs
TB = sym('T_B%d',[3,1]);

%% Objects constructed from variables

% Direction cosine matrix I to B (assumes scalar-first quaternion)
% CBI = diag([q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2, q(1)^2 + q(3)^2 - q(2)^2 - q(4)^2, ...
%             q(1)^2 + q(4)^2 - q(2)^2 - q(3)^2]);
% CBI(1,2) = 2*( q(2)*q(3) + q(4)*q(1) );
% CBI(1,3) = 2*( q(2)*q(4) - q(3)*q(1) );
% CBI(2,1) = 2*( q(3)*q(2) - q(4)*q(1) );
% CBI(2,3) = 2*( q(3)*q(4) + q(2)*q(1) );
% CBI(3,1) = 2*( q(4)*q(2) + q(3)*q(1) );
% CBI(3,2) = 2*( q(4)*q(3) - q(2)*q(1) );
%
% display(CBI)

CBI = qlib.q_dcm(q);

display(CBI)

save('sym_data/sym_var','s','m','rI','vI','q','omegaB','TB');
save('sym_data/sym_var_obj','CBI');
save('sym_data/sym_param','c_ax','c_ayz','CA','JBvec','JB','JBinv','gI','rho','S_A','rTB','rcpB','alphmdt','betmdt');
