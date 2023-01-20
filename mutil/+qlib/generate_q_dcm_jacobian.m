% Symbolically generate function for evaluating Jacobian of DCM wrt quaternion

clearvars 
close all 
clc

% Quaternion with scalar-last convention
qBI = sym('qBI_%d',[4,1]);
assume(qBI,'real');

% BULCG -> NED
CIB = transpose(qlib.q_dcm(qBI)); 

dCIB_dqBI_vec = simplify(jacobian(CIB(:),qBI),'Steps',100);
dCIB_dqBI = sym('dcIB_dqBI_%d%d%d',[3,3,4]);
for j = 1:4
    dCIB_dqBI(:,:,j) = reshape(dCIB_dqBI_vec(:,j),[3,3]);
end

matlabFunction(dCIB_dqBI,'File','q_dcm_transpose_jacobian.m','Vars',{qBI});

% BULCG <- NED
CBI = qlib.q_dcm(qBI); 

dCBI_dqBI_vec = simplify(jacobian(CBI(:),qBI),'Steps',100);
dCBI_dqBI = sym('dcBI_dqBI_%d%d%d',[3,3,4]);
for j = 1:4
    dCBI_dqBI(:,:,j) = reshape(dCBI_dqBI_vec(:,j),[3,3]);
end

matlabFunction(dCBI_dqBI,'File','q_dcm_jacobian.m','Vars',{qBI});

% Numerical validation:

q = qlib.q_rand_unit();
pert = 1e-8;
I = eye(4);
dCIB_dqBI_num = zeros(3,3,4);
for j = 1:4
    dC_p = transpose(qlib.q_dcm(q+pert*I(:,j)));
    dC_m = transpose(qlib.q_dcm(q-pert*I(:,j)));
    dCIB_dqBI_num(:,:,j) = (dC_p-dC_m)/(2*pert);
end

dCIB_dqBI_exact = qlib.q_dcm_transpose_jacobian(q);

for j = 1:4
    fprintf("Slice %d error: %.3e\n",j,norm(dCIB_dqBI_exact(:,:,j) - dCIB_dqBI_num(:,:,j)))
end
