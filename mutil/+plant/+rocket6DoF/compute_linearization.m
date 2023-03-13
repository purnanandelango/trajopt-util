function [A,B,S,w] = compute_linearization(z_,u_,s_,c_ax,c_ayz,JBvec,gI,rho,SA,rTB,rcpB,alphmdt,betmdt)
% Linearization of RHS of the nonlinear ODE describing airborne rigid vehicle in lower atmosphere
% Normalize quaternion
z_(8:11) = z_(8:11)./norm(z_(8:11),2);

A = plant.rocket6DoF.compute_A(z_,u_,s_,c_ax,c_ayz,JBvec,rho,SA,rcpB);
B = plant.rocket6DoF.compute_B(z_,u_,s_,JBvec,rTB,alphmdt);
S = plant.rocket6DoF.compute_S(z_,u_,c_ax,c_ayz,JBvec,gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
w = plant.rocket6DoF.dyn_func(z_,u_,s_,c_ax,c_ayz,JBvec,gI,rho,SA,rTB,rcpB,alphmdt,betmdt) - A*z_ - B*u_ - S*s_; 

end