function [q_x_star] = q_skew_star(q)
% Compute skew operation star matrix
% SCALAR-LAST convention
q0 = q(4);
q1 = q(1);
q2 = q(2);
q3 = q(3);

q_x_star = [q0, q3,  -q2, q1;
           -q3,  q0,  q1, q2;
            q2, -q1,  q0, q3;
           -q1, -q2, -q3, q0];
end