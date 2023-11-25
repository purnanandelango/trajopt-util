function dx = dyn_func(x,u,s,n,c_d,accl)
% RHS of ODE describing 2D or 3D double integrator    
% c_d is the drag coefficient
% accl is a constant external acceleration (e.g. gravity)

    assert(length(x) == 2*n,"State dimension does not match n");
    assert(ismember(n,[2,3]));    
    assert(s > 0,"Dilation factor should be positive.");
    assert(c_d >= 0,"Drag coefficient should be positive.");


    % r = x(1:n);       % Position
    v = x(n+1:2*n);     % Velocity

    dx = s*[v;
            u + accl - c_d*norm(v)*v];

end