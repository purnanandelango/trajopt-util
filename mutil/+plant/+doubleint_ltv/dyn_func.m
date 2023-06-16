function dx = dyn_func(t, x, u, ...
                       T, n, accl)
% RHS of ODE describing 2D or 3D time-varying double integrator    
% T is the total time of maneuver
% accl is a constant external acceleration vector

    assert(length(x) == 2*n,"State dimension does not match n.");
    assert(ismember(n,[2,3]),"Only 2D and 3D double integrators are supported.");

    dx = plant.doubleint_ltv.Amat(t,T,n)*x + plant.doubleint_ltv.Bmat(t,T,n)*u + accl;

end