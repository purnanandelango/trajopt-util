function [h,dh] = q_aoa_cnstr(vI,qBI,vmax,cosaoamax,flag)
% Air-speed triggered angle-of-attack constraint parameters
% Two version of STC can be selected

    CBI         = qlib.q_dcm(qBI);
    dCBI_dqBI   = qlib.q_dcm_jacobian(qBI);
    nrmvI       = norm(vI);
    e1          = [1;0;0]; 

    g   = vmax^2 - nrmvI^2;                                                                                 % Trigger function
    dg  = [zeros(1,4), -2*vI', zeros(1,7)];                                                                 % Trigger gradient
    c   = cosaoamax*nrmvI + e1'*CBI*vI;                                                                     % Constraint function
    dc  = [zeros(1,4), (cosaoamax*vI'/nrmvI+e1'*CBI), e1'*linalg.sltimes(dCBI_dqBI,vI), zeros(1,3)];        % Constraint gradient
    
    gminus = min(0,g);

    if flag == "v1" % New modification
        h   = (gminus^2)*c;                                                                                 % STC function
        dh  = 2*gminus*c*dg + (gminus^2)*dc;                                                                % STC gradient          
    elseif flag == "v2" % Proposed by Smzuk et al.
        h   = -gminus*c;                                                                                    % STC function
        if g <= 0
            dh = -dg*c - g*dc;                                                                              % STC gradient
        else
            dh = zeros(1,14);                                                                               % STC gradient
        end
    else
        error("Incorrect STC version flag.");
    end

end