function [A,B,w] = compute_linearization(t,x,u,T,n)
% Linearization of RHS of ODE describing 2D or 3D time-varying double integrator

    f = plant.doubleint_ltv.dyn_func(t,x,u,T,n);
    A = plant.doubleint_ltv.Amat(t,T,n);
    B = plant.doubleint_ltv.Bmat(t,T,n);

    w = f - A*x - B*u;

end