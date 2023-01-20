clearvars
clc

load sym_data/sym_param.mat
load sym_data/sym_var.mat
load sym_data/sym_var_obj.mat

plant.rocket6DoF.set_assumptions;

% Mass depletion
mprime = simplify(s*(-alphmdt*sqrt(sum(TB.*TB)) - betmdt),'Steps',100);

% Kinematics
rIprime = s*vI;

% Translational dynamics
AB = simplify(-0.5*rho*S_A*sqrt(sum(vI.*vI))*(CA*CBI)*vI,'Steps',100);
FI = (CBI')*(TB + AB);
vIprime = simplify(s*FI/m + s*gI,'Steps',100);

% Quaternion dynamics (scalar-last)
qprime = simplify(0.5*s*qlib.q_mul(q,[omegaB;0]),'Steps',100);                    

% Rotational dynamics
omegaBprime = simplify(s*JBinv*(cross(rTB,TB) + cross(rcpB,AB) - cross(omegaB,JB*omegaB)),'Steps',100);

save('sym_data/sym_dyn','mprime','rIprime','vIprime','qprime','omegaBprime');
