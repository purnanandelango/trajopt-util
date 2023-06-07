function [dx,fun1,fun2] = dyn_func(t,x,u,T,n)
% RHS of ODE describing 2D or 3D time-varying double integrator    
% T is the total time of maneuver

    assert(length(x) == 2*n,"State dimension does not match n");
    assert(ismember(n,[2,3]));

    %%% TIME VARIATION
    tau         = pi*t/T;
    fun1        = (1+cos(2*tau))/4;
    fun2        = 1+(sin(4*tau))/6; 
    %%%

    % r = x(1:n);       % Position
    v = x(n+1:2*n);     % Velocity

    dx = [v;
          fun1*v + fun2*u];

end