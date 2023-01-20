function r = q_mul(q,p)
% Compute product of two quaternions
% SCALAR-LAST convention

q0 = q(4);
qv = q(1:3);

p0 = p(4);
pv = p(1:3);

r = [q0*pv + p0*qv + qlib.skew(qv)*pv ; q0*p0 - qv'*pv];

end