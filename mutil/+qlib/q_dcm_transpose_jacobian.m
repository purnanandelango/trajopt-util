function dCIB_dqBI = q_dcm_transpose_jacobian(in1)
% Compute Jacobian of transpose of DCM wrt quaternion
% Q_DCM_TRANSPOSE_JACOBIAN
%    dCIB_dqBI = Q_DCM_TRANSPOSE_JACOBIAN(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    19-Jan-2023 13:58:45

qBI_1 = in1(1,:);
qBI_2 = in1(2,:);
qBI_3 = in1(3,:);
qBI_4 = in1(4,:);
t2 = qBI_1.*2.0;
t3 = qBI_2.*2.0;
t4 = qBI_3.*2.0;
t5 = qBI_4.*2.0;
t6 = -t2;
t7 = -t3;
t8 = -t4;
t9 = -t5;
dCIB_dqBI = reshape([t2,t3,t4,t3,t6,t5,t4,t9,t6,t7,t2,t9,t2,t3,t4,t5,t4,t7,t8,t5,t2,t9,t8,t3,t2,t3,t4,t5,t4,t7,t8,t5,t2,t3,t6,t5],[3,3,4]);
