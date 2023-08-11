function AoA = compute_AoA(vI,qBI)
% Compute angle of attack from inertial velocity and body-axis quaternion
    N = size(vI,2);
    AoA = zeros(1,N);
    for k = 1:N
        CBI = qlib.q_dcm(qBI(:,k));
        vB = CBI*vI(:,k);
        AoA(k) = acosd(-vB(1)/norm(vB));
    end
end