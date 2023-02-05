clearvars
close all
clc

rsc = sym('rsc_%d',[3,1]);
vsc = sym('vsc_%d',[3,1]);
rE = sym('rE_%d',[3,1]);
rS = sym('rS_%d',[3,1]);
vE = sym('vE_%d',[3,1]);
vS = sym('vS_%d',[3,1]);
syms GME GMM GMS SRP

drscdt         = vsc;
dvscdt         = -GMM*rsc/norm(rsc)^3 + GME*( (rE-rsc)/norm(rE-rsc)^3 - rE/norm(rE)^3 ) + GMS*( (rS-rsc)/norm(rS-rsc)^3 - rS/norm(rS)^3 ) ...
                 -SRP*GMS*(rS-rsc)/norm(rS-rsc)^3 ...
                 +;

jac_x_dxdt     = simplify(jacobian([drscdt;dvscdt],[rsc;vsc]),Steps=100);
jac_t_dxdt     = [zeros(3,1);
                  simplify(jacobian(dvscdt,rE)*vE,Steps=100) + simplify(jacobian(dvscdt,rS)*vS,Steps=100)];

matlabFunction(jac_t_dxdt,jac_x_dxdt,'File','dyn_fun_jacobian_inert.m',)
