function q_cross_mat = q_cross(q)
% Compute quaternion cross-product matrix
% SCALAR-LAST convention

qv = q(1:3);
q0 = q(4);

q_cross_mat = [q0*eye(3) + skew(qv), qv ; zeros(1,4)];

end