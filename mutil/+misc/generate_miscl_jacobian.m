clearvars
close all
clc

% Verify alternate expressions for the RHS of the quaternion dynamics

omg_rand = randn(3,1);
q_rand = qlib.q_rand_unit();

a = qlib.q_skew_star([omg_rand;0])*q_rand;
b = qlib.q_mul(q_rand,[omg_rand;0]);

norm(a-b)

% Jacobian of q_skew_star([w;0]) with respect to w
% This terms appears in the quaternion dynamics

omg = sym('omg_%d',[3,1]);
assume(omg,'real');
Omg = qlib.q_skew_star([omg;0]);
dOmg_omg_vec = double(jacobian(Omg(:),omg));
dOmg_omg = zeros(4,4,3);
for j = 1:3
    dOmg_omg(:,:,j) = reshape(dOmg_omg_vec(:,j),[4,4]);
end

norm(misc.q_skew_star_jacobian_times_q(q_rand) - linalg.sltimes(dOmg_omg,q_rand))

% Jacobian of skew(omg)*A*omg with respect to w
% where A is a constant matrix
% This term appears in the rotational dynamics

A = sym('A_%d%d',[3,3]);
assume(A,'real');
cross_term_jacobian = simplify(jacobian(cross(omg,A*omg),omg),'Steps',100);

matlabFunction(cross_term_jacobian,'File','omega_cross_term_jacobian.m','Vars',{omg,A});