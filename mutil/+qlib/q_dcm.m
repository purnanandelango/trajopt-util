function dcm = q_dcm(q)
% Compute the direction cosine matrix (DCM)
% SCALAR-LAST convention

q0 = q(4);
q1 = q(1);
q2 = q(2);
q3 = q(3);

dcm = [q0^2+q1^2-q2^2-q3^2,  2*(q1*q2-q0*q3), 2*(q0*q2+q1*q3);
       2*(q0*q3+q1*q2), q0^2-q1^2+q2^2-q3^2, 2*(q2*q3-q0*q1);
       2*(q1*q3-q0*q2), 2*(q0*q1+q2*q3), q0^2-q1^2-q2^2+q3^2]';
   
   
% verification:
% dcm1 = q_skew_star(a)*q_skew(q_conj(a))
% dcm = dcm1(1:3,1:3)

end