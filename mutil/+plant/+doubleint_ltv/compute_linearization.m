function [A,B,w] = compute_linearization(t,x,u,T,n)
% Linearization of RHS of ODE describing 2D or 3D time-varying double integrator

    [f,fun1,fun2] = plant.doubleint_ltv.dyn_func(t,x,u,T,n);

    A = [zeros(n) eye(n);
         zeros(n) eye(n)*fun1];

    B = [zeros(n);
         eye(n)*fun2];

    w = f - A*x - B*u;

end