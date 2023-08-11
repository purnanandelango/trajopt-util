function [h,dh] = q_aoa_cnstr(vI,qBI,vmax,cosaoamax)
% Air-speed triggered angle-of-attack constraint parameters

    CBI         = qlib.q_dcm(qBI);
    dCBI_dqBI   = qlib.q_dcm_jacobian(qBI);
    nrmvI       = norm(vI);
    e1          = [1;0;0]; 

    g   = vmax^2 - nrmvI^2;                                                                                 % Trigger function
    dg  = [zeros(1,4), -2*vI', zeros(1,7)];                                                                 % Trigger gradient
    c   = cosaoamax*nrmvI + e1'*CBI*vI;                                                                     % Constraint function
    dc  = [zeros(1,4), (cosaoamax*vI'/nrmvI+e1'*CBI), e1'*linalg.sltimes(dCBI_dqBI,vI), zeros(1,3)];     % Constraint gradient
    
    gminus = min(0,g);

    h   = (gminus^2)*c;                                                                                     % STC function
    dh  = 2*gminus*c*dg + (gminus^2)*dc;                                                                    % STC gradient        
end