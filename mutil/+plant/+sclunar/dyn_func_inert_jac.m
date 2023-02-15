% Compute Jacobian of the spacraft dynamics in N body problem in the J2000 MEME inertial frame (centered at Moon)
function [jac_t_dxdt,jac_x_dxdt] = dyn_func_inert_jac(t,x,astro)
%dyn_func_inert_jac
%    [JAC_T_DXDT,JAC_X_DXDT] = dyn_func_inert_jac(t,x,astro)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    10-Feb-2023 13:21:09

t_ephm      = astro.start_JD + t*astro.t_star/(24*3600);     % [day]
ephm_state  = plant.sclunar.ephemeris_state(0,t_ephm,astro.Earth,astro.Moon);
rE          = ephm_state(1:3)'/astro.l_star;            % Normalized Earth position
vE          = ephm_state(4:6)'/astro.v_star;            % Normalized Earth velocity
ephm_state  = plant.sclunar.ephemeris_state(0,t_ephm,astro.Sun,astro.Moon);
rS          = ephm_state(1:3)'/astro.l_star;            % Normalized Sun position
vS          = ephm_state(4:6)'/astro.l_star;            % Normalized Sun velocity
GME         = astro.ndGM_Earth;
GMM         = astro.ndGM_Moon;
GMS         = astro.ndGM_Sun;

rE_1 = rE(1,:);
rE_2 = rE(2,:);
rE_3 = rE(3,:);
rS_1 = rS(1,:);
rS_2 = rS(2,:);
rS_3 = rS(3,:);
rsc_1 = x(1,:);
rsc_2 = x(2,:);
rsc_3 = x(3,:);
vE_1 = vE(1,:);
vE_2 = vE(2,:);
vE_3 = vE(3,:);
vS_1 = vS(1,:);
vS_2 = vS(2,:);
vS_3 = vS(3,:);
t2 = rE_1.^2;
t3 = rE_2.^2;
t4 = rE_3.^2;
t5 = rS_1.^2;
t6 = rS_2.^2;
t7 = rS_3.^2;
t8 = rsc_1.^2;
t9 = rsc_2.^2;
t10 = rsc_3.^2;
t11 = rE_1.*rsc_1.*2.0;
t12 = rE_2.*rsc_2.*2.0;
t13 = rE_3.*rsc_3.*2.0;
t14 = rS_1.*rsc_1.*2.0;
t15 = rS_2.*rsc_2.*2.0;
t16 = rS_3.*rsc_3.*2.0;
t17 = -rsc_1;
t18 = -rsc_2;
t19 = -rsc_3;
t20 = t8.*2.0;
t21 = t9.*2.0;
t22 = t10.*2.0;
t23 = -t11;
t24 = -t12;
t25 = -t13;
t26 = -t14;
t27 = -t15;
t28 = -t16;
t32 = rE_1+t17;
t33 = rE_2+t18;
t34 = rE_3+t19;
t35 = rS_1+t17;
t36 = rS_2+t18;
t37 = rS_3+t19;
t44 = t2+t3+t4;
t45 = t5+t6+t7;
t46 = t8+t9+t10;
t29 = -t20;
t30 = -t21;
t31 = -t22;
t38 = t32.^2;
t39 = t33.^2;
t40 = t34.^2;
t41 = t35.^2;
t42 = t36.^2;
t43 = t37.^2;
t47 = 1.0./t44.^(3.0./2.0);
t48 = 1.0./t44.^(5.0./2.0);
t49 = 1.0./t45.^(3.0./2.0);
t50 = 1.0./t45.^(5.0./2.0);
t51 = 1.0./t46.^(3.0./2.0);
t52 = 1.0./t46.^(5.0./2.0);
t53 = GMM.*t51;
t54 = rE_1.*rE_2.*t48.*3.0;
t55 = rE_1.*rE_3.*t48.*3.0;
t56 = rE_2.*rE_3.*t48.*3.0;
t57 = rS_1.*rS_2.*t50.*3.0;
t58 = rS_1.*rS_3.*t50.*3.0;
t59 = rS_2.*rS_3.*t50.*3.0;
t67 = GMM.*rsc_1.*rsc_2.*t52.*3.0;
t68 = GMM.*rsc_1.*rsc_3.*t52.*3.0;
t69 = GMM.*rsc_2.*rsc_3.*t52.*3.0;
t70 = t38+t39+t40;
t71 = t41+t42+t43;
t60 = -t53;
t61 = -t54;
t62 = -t55;
t63 = -t56;
t64 = -t57;
t65 = -t58;
t66 = -t59;
t72 = 1.0./t70.^(3.0./2.0);
t73 = 1.0./t70.^(5.0./2.0);
t74 = 1.0./t71.^(3.0./2.0);
t75 = 1.0./t71.^(5.0./2.0);
t76 = -t72;
t77 = -t74;
t78 = t32.*t33.*t73.*3.0;
t79 = t32.*t34.*t73.*3.0;
t80 = t33.*t34.*t73.*3.0;
t81 = t35.*t36.*t75.*3.0;
t82 = t35.*t37.*t75.*3.0;
t83 = t36.*t37.*t75.*3.0;
mt1 = [0.0;0.0;0.0;GME.*vE_2.*(t54-t78)+GME.*vE_3.*(t55-t79)+GMS.*vS_2.*(t57-t81)+GMS.*vS_3.*(t58-t82)-GME.*vE_1.*(t47+t76-t2.*t48.*3.0+t38.*t73.*3.0)-GMS.*vS_1.*(t49+t77-t5.*t50.*3.0+t41.*t75.*3.0);GME.*vE_1.*(t54-t78)+GME.*vE_3.*(t56-t80)+GMS.*vS_1.*(t57-t81)+GMS.*vS_3.*(t59-t83)-GME.*vE_2.*(t47+t76-t3.*t48.*3.0+t39.*t73.*3.0)-GMS.*vS_2.*(t49+t77-t6.*t50.*3.0+t42.*t75.*3.0)];
mt2 = [GME.*vE_1.*(t55-t79)+GME.*vE_2.*(t56-t80)+GMS.*vS_1.*(t58-t82)+GMS.*vS_2.*(t59-t83)-GME.*vE_3.*(t47+t76-t4.*t48.*3.0+t40.*t73.*3.0)-GMS.*vS_3.*(t49+t77-t7.*t50.*3.0+t43.*t75.*3.0)];
jac_t_dxdt = [mt1;mt2];
if nargout > 1
    t84 = GME.*t78;
    t85 = GME.*t79;
    t86 = GME.*t80;
    t87 = GMS.*t81;
    t88 = GMS.*t82;
    t89 = GMS.*t83;
    t90 = t61+t78;
    t91 = t62+t79;
    t92 = t63+t80;
    t93 = t64+t81;
    t94 = t65+t82;
    t95 = t66+t83;
    t96 = t67+t84+t87;
    t97 = t68+t85+t88;
    t98 = t69+t86+t89;
    jac_x_dxdt = reshape([0.0,0.0,0.0,t60+GMM.*t8.*t52.*3.0-GME.*t73.*(t2.*-2.0+t3+t4+t9+t10+t24+t25+t29+rE_1.*rsc_1.*4.0)-GMS.*t75.*(t5.*-2.0+t6+t7+t9+t10+t27+t28+t29+rS_1.*rsc_1.*4.0),t96,t97,0.0,0.0,0.0,t96,t60+GMM.*t9.*t52.*3.0-GME.*t73.*(t2-t3.*2.0+t4+t8+t10+t23+t25+t30+rE_2.*rsc_2.*4.0)-GMS.*t75.*(t5-t6.*2.0+t7+t8+t10+t26+t28+t30+rS_2.*rsc_2.*4.0),t98,0.0,0.0,0.0,t97,t98,t60+GMM.*t10.*t52.*3.0-GME.*t73.*(t2+t3-t4.*2.0+t8+t9+t23+t24+t31+rE_3.*rsc_3.*4.0)-GMS.*t75.*(t5+t6-t7.*2.0+t8+t9+t26+t27+t31+rS_3.*rsc_3.*4.0),1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0],[6,6]);
end
