clearvars
close all
clc

srp_flg = true;
mJ2_flg = false;

rsc = sym('rsc_%d',[3,1]);  assume(rsc,'real');
vsc = sym('vsc_%d',[3,1]);  assume(vsc,'real');
rE = sym('rE_%d',[3,1]);    assume(rE,'real');
rS = sym('rS_%d',[3,1]);    assume(rS,'real');
vE = sym('vE_%d',[3,1]);    assume(vE,'real');
vS = sym('vS_%d',[3,1]);    assume(vS,'real');
syms GME GMM GMS SRP MJ2 RM Meqincl

assume(GME,'positive');
assume(GMM,'positive');
assume(GMS,'positive');
assume(SRP,'positive');
assume(MJ2,'positive');
assume(RM,'positive');
assume(Meqincl,'positive');

x = [rsc;vsc];

vEprp          = -cross(rE,cross(rE,vE));
rscEM          = rsc - dot(rsc,vEprp)*vEprp/norm(vEprp)^2;
lam_sc         = acos( dot(rE,rscEM)/(norm(rE)*norm(rscEM)) ) + Meqincl;

drscdt         = vsc;
dvscdt         = -GMM*rsc/norm(rsc)^3 ...
                 + GME*( (rE-rsc)/norm(rE-rsc)^3 - rE/norm(rE)^3 ) ...
                 + GMS*( (rS-rsc)/norm(rS-rsc)^3 - rS/norm(rS)^3 );

a_SRP          = -SRP*(rS-rsc)/norm(rS-rsc)^3;
a_J2           = 1.5*GMM*MJ2*(RM^2)*( 3*sin(lam_sc)^2 - 1 )*rsc/norm(rsc)^5;

if srp_flg; dvscdt = dvscdt + a_SRP; end
if mJ2_flg; dvscdt = dvscdt + a_J2; end

jac_x_dxdt     = simplify(jacobian([drscdt;dvscdt],[rsc;vsc]),Steps=100);
jac_t_dxdt     = [zeros(3,1);
                  simplify(jacobian(dvscdt,rE)*vE,Steps=100) + simplify(jacobian(dvscdt,rS)*vS,Steps=100)];

file_name = "dyn_func_inert";
arg_list = {x,rE,vE,rS,vS,'GME','GMM','GMS'};
if srp_flg
    file_name = file_name + "_SRP";
    arg_list{end+1} = 'SRP';
    if mJ2_flg
        arg_list(end+1:end+3) = {'MJ2','RM','Meqincl'};
        file_name = file_name + "_J2";
    end
end

matlabFunction(jac_t_dxdt,jac_x_dxdt,'File',file_name+"_jac.m",'Vars',arg_list);
