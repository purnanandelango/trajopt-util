function uBI = compute_thrustdir(uB,qBI)
% Compute orientation of thrust vector (uB) in inertial coordinates : uBI
    K = size(qBI,2);
    uBI = zeros(4,K);
    for k = 1:K
       uBI(:,k) = qlib.q_mul(qlib.q_mul(qBI(:,k),[uB(:,k);0]),qlib.q_conj(qBI(:,k)));
       
       % Equivalent to:
       % uBI(1:3,k) = qlib.q_dcm(qBI(:,k))'*uB(:,k);
       
       % Sanity check
       % display(norm(qlib.q_mul(qlib.q_mul(qBI(:,k),[uB(:,k);0]),qlib.q_conj(qBI(:,k))) - [qlib.q_dcm(qBI(:,k))'*uB(:,k);0]))
    end
    uBI(4,:) = []; % Since uBI is a pure quaternion, remove scalar part
end