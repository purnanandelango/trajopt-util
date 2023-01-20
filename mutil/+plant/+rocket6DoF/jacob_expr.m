clearvars
clc

load sym_data/sym_param.mat
load sym_data/sym_var.mat
load sym_data/sym_var_obj.mat
load sym_data/sym_dyn.mat

plant.rocket6DoF.set_assumptions;

z_ = [m;rI;vI;q;omegaB];
u_ = TB;
s_ = s;
zprime_ = [mprime;rIprime;vIprime;qprime;omegaBprime];

A = simplify(jacobian(zprime_,z_),'Steps',100);

%%

B = simplify(jacobian(zprime_,u_),'Steps',100);

%%

S = simplify(jacobian(zprime_,s_),'Steps',100);

%%

matlabFunction(zprime_,'File','dyn_func.m','Vars',{[m;rI;vI;q;omegaB]...
    ,TB,'s','c_ax','c_ayz',JBvec,gI,'rho','S_A',rTB...
    ,rcpB,'alphmdt','betmdt'});

% z,u,s,c_ax,c_ayz,JBvec,rho,SA,rcpB
matlabFunction(A,'File','compute_A.m','Vars',{[m;rI;vI;q;omegaB]...
    ,TB,'s','c_ax','c_ayz',JBvec,'rho','S_A',rcpB});

% z,u,s,JBvec,rTB,alphmdt
matlabFunction(B,'File','compute_B.m','Vars',{[m;rI;vI;q;omegaB],TB,'s',...
    JBvec,rTB,'alphmdt'});

% z,u,c_ax,c_ayz,JBvec,gI,rho,SA,rTB,rcpB,alphmdt,betmdt
matlabFunction(S,'File','compute_S.m','Vars',{[m;rI;vI;q;omegaB]...
    ,TB,'c_ax','c_ayz',JBvec,gI,'rho','S_A',rTB...
    ,rcpB,'alphmdt','betmdt'});