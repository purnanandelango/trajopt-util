% Testing the symbolically obtained linearization

clearvars
clc

% Random operating point

%%% IMPORTANT!
% Suppress quaternion normalization in sysdyn_func and sysdyn_func_naive before testing

q_ = randn(4,1); q_ = q_/norm(q_);
z_mean = [1.5;4*rand(3,1);2*rand(3,1);q_;0.2*rand(3,1)];
% z_ = mvnrnd(z_mean,eye(14))';
z_ = z_mean;
TB = [1;2;3];
s = 4;
c_ax = 0.5; c_ayz = 1; CA = diag([c_ax,c_ayz,c_ayz]);
JB = 0.168*diag([2e-2,1,0.5]);
gI = [-1;0;0];
rho = 1;
SA = 0.5;
rTB = -0.25*[1;0;0];
rcpB = 0.05*[1;0;0];
alphmdt = 1/30;
betmdt = 0;

tic
zprime_1 = plant.rocket6DoF.dyn_func(z_,TB,s,c_ax,c_ayz,diag(JB),gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
toc

tic
zprime_2 = plant.rocket6DoF.dyn_func_naive(z_,TB,s,CA,JB,gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
toc

fprintf('Error in output from sysdyn_func and sysdyn_func_naive:\n')
display(norm(zprime_1-zprime_2,inf))

% testing linearization
tic
[A,B,S,w] = plant.rocket6DoF.compute_linearization(z_,TB,s,c_ax,c_ayz,diag(JB),gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
toc

tic
[A2,B2,S2,w2] = numJacobian(@plant.rocket6DoF.dyn_func,z_,TB,s,c_ax,c_ayz,diag(JB),gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
toc

fprintf('\nError in A:\n')
display(max(max(abs(A-A2))))

fprintf('\nError in B:\n')
display(max(max(abs(B-B2))))

fprintf('\nError in S:\n')
display(max(max(abs(S-S2))))

%% function

function [A,B,S,w] = numJacobian(fun,z_,u_,s_,c_ax,c_ayz,JBvec,gI,rho,SA,rTB,rcpB,alphmdt,betmdt)

f = fun(z_,u_,s_,c_ax,c_ayz,JBvec,gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
pert = 0.001*min(abs(f));
if pert<1e-10
    error("Required perturbation is too small")
end


% compute A
pIz = pert*eye(14);
A = zeros(14,14);
for i = 1:14
    fp  = fun(z_+pIz(:,i),u_,s_,c_ax,c_ayz,JBvec,gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
    fm  = fun(z_-pIz(:,i),u_,s_,c_ax,c_ayz,JBvec,gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
    A(:,i) = (fp-fm)/(2*pert);
end

% compute B
pIu = pert*eye(3);
B = zeros(14,3);
for i = 1:3
    fp  = fun(z_,u_+pIu(:,i),s_,c_ax,c_ayz,JBvec,gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
    fm  = fun(z_,u_-pIu(:,i),s_,c_ax,c_ayz,JBvec,gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
    B(:,i) = (fp-fm)/(2*pert);
end

% compute S
pIs = pert;
fp  = fun(z_,u_,s_+pIs,c_ax,c_ayz,JBvec,gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
fm  = fun(z_,u_,s_-pIs,c_ax,c_ayz,JBvec,gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
S = (fp-fm)/(2*pert);

w = -A*z_-B*u_;

end