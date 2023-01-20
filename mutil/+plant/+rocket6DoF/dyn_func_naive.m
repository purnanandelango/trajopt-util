function zprime_ = dyn_func_naive(z_,TB,s,CA,JB,gI,rho,SA,rTB,rcpB,alphmdt,betmdt)
% RHS of the nonlinear ODE describing airborne rigid vehicle in lower atmosphere
% Naive function without code optimization
%
% z_ = [m;rI;vI;q;omegaB]
% CA = diag([c_ax,c_ayz,c_ayz])
% JB = diag([J_B1,J_B2,J_B3])

m = z_(1);
% rI = z_(2:4);
vI = z_(5:7);
q = z_(8:11); % scalar-last

% Normalize quaternion
% q = q./norm(q)

omegaB = z_(12:14);

JBinv = diag(1./diag(JB));

% Direction cosine matrix I to B (assumes scalar-first convention)
% CBI = diag([q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2, q(1)^2 + q(3)^2 - q(2)^2 - q(4)^2, ...
%             q(1)^2 + q(4)^2 - q(2)^2 - q(3)^2]);
% CBI(1,2) = 2*( q(2)*q(3) + q(4)*q(1) );
% CBI(1,3) = 2*( q(2)*q(4) - q(3)*q(1) );
% CBI(2,1) = 2*( q(3)*q(2) - q(4)*q(1) );
% CBI(2,3) = 2*( q(3)*q(4) + q(2)*q(1) );
% CBI(3,1) = 2*( q(4)*q(2) + q(3)*q(1) );
% CBI(3,2) = 2*( q(4)*q(3) - q(2)*q(1) );

CBI = qlib.q_dcm(q);

mprime = s*(-alphmdt*sqrt(sum(TB.*TB)) - betmdt);
rIprime = s*vI;

AB = -0.5*rho*SA*sqrt(sum(vI.*vI))*(CA*CBI)*vI;
FI = (CBI')*(TB + AB);
vIprime = s*FI/m + s*gI;

qprime = 0.5*s*qlib.q_mul(q,[omegaB;0]);  
         
omegaBprime = s*JBinv*(cross(rTB,TB) + cross(rcpB,AB) - cross(omegaB,JB*omegaB));

zprime_ = [mprime;rIprime;vIprime;qprime;omegaBprime];

end