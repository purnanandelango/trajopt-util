function AI = compute_dragdir(drag,vI,qBI)
% Compute the orientation of aerodynamic force vector in inertial coordinates
    K = size(vI,2);
    rho = drag.rho;
    SA = drag.SA;
    CA = drag.CA;    
    AI(3,K) = 0;
    for k = 1:K
        CBI = qlib.q_dcm(qBI(:,k));
        AI(:,k) = -0.5*rho*SA*sqrt(sum(vI(:,k).*vI(:,k)))*(CBI')*(CA*CBI)*vI(:,k);
    end
end