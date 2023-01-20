function q_cross_star_mat = q_cross_star(q)
% Compute quaternion cross-product star matrix
% SCALAR-LAST convention

qv = q(1:3);
q0 = q(4);

q_cross_star_mat = [q0*eye(3) - skew(qv), qv ; zeros(1,4)];

end