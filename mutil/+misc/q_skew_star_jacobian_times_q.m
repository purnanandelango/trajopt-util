function q_out = q_skew_star_jacobian_times_q(q)
% SCALAR-LAST convention is followed
%
% The output of this function is the result of the product:
%
%       jacobian(q_skew_star([omg;0]),omg)*q
%
% where omg is any vector in R^3
%

dOmg_omg_vec = [0	0	0;
                0	0  -1;
                0	1	0;
               -1	0	0;
                0	0	1;
                0	0	0;
               -1	0	0;
                0  -1	0;
                0  -1	0;
                1	0	0;
                0	0	0;
                0	0  -1;
                1	0	0;
                0	1	0;
                0	0	1;
                0	0	0];

q_out = zeros(4,3);
for j = 1:3
    q_out(:,j) = reshape(dOmg_omg_vec(:,j),[4,4])*q;
end

end