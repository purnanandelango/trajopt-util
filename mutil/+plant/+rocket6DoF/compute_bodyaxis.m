function xBI = compute_bodyaxis(qBI)
% Compute body axis orientation (xB) in interial coordinates : xBI
% Note xI = [1;0;0]
K = size(qBI,2);
xBI = zeros(4,K);
for k = 1:K
   xBI(:,k) = qlib.q_mul(qlib.q_mul(qBI(:,k),[1;0;0;0]),qlib.q_conj(qBI(:,k)));
   
   % Equivalent to:
   % xBI(1:3,k) = qlib.q_dcm(qBI(:,k))'*[1;0;0];
   
   % Sanity check
   % display(norm(qlib.q_mul(qlib.q_mul(qBI(:,k),[1;0;0;0]),qlib.q_conj(qBI(:,k))) - [qlib.q_dcm(qBI(:,k))'*[1;0;0];0]))
end
xBI(4,:) = []; % Since xBI is a pure quaternion, remove scalar part
end