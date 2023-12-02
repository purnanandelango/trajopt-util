function [data,cones,param] = ecos2scs(c,G,h,dims,A,b)
% Construct SCS input using ECOS input
    data = struct;
    cones = struct;
    param = struct;

    data.c = c;
    data.A = [A;G];
    data.b = [b;h];
    
    cones.z = size(A,1);    % Size of zero cone
    cones.l = dims.l;       % Size of non-negative orthant cone
    cones.q = dims.q;       % Sizes of second order cones

    param.eps_abs = 1e-8;
    param.eps_rel = 1e-8;
end