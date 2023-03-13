function [A,B,S,w] = compute_linearization(x,u,s,n,coeff_drag,g)
% Linearization of RHS of ODE describing 2D or 3D double integrator

    f = plant.doubleint.dyn_func(x,u,s,n,coeff_drag,g);

    v = reshape(x(n+1:2*n),[n,1]);

    A = [zeros(n) s*eye(n);
         zeros(n) -s*coeff_drag*(norm(v)*eye(n) + v*v'/norm(v))];

    B = [zeros(n);
         s*eye(n)];

    S = [v;
         u + g - coeff_drag*norm(v)*v];

    w = f - A*x - B*u - S*s;

end