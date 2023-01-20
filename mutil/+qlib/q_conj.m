function q_str = q_conj(q)
% Compute conjugate of quaternion
% SCALAR_LAST convention

q0 = q(4);
qv = q(1:3);

q_str = [-qv;q0];
end