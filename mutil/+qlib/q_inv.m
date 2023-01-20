function qinv = q_inv(q)
% Compute quaternion inverse
% SCALAR-LAST convention

qinv = q_conj(q)/norm(q);

% for unit quaternions q_conj and q_inv are the same thing
end