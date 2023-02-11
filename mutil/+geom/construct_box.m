function [H,g] = construct_box(z_min,z_max)
% Construct a linear inequality representation of infinity norm ball
% { z | z_min <= z <= z_max }
    n = length(z_min);
    In = eye(n);
    z_min = reshape(z_min,[n,1]);
    z_max = reshape(z_max,[n,1]);
    H = [In;-In];
    g = [z_max;-z_min];
end
