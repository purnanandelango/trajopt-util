function [A,B,S,w] = compute_linearization(x,u,s,n,c_d,accl)
% Linearization of RHS of ODE describing 2D or 3D double integrator

    f = plant.doubleint.dyn_func(x,u,s,n,c_d,accl);

    v = reshape(x(n+1:2*n),[n,1]);

    if norm(v) > 1e-7
        term_v = v*v'/norm(v);
    else
        term_v = 0;
    end
    A = [zeros(n) s*eye(n);
         zeros(n) -s*c_d*(norm(v)*eye(n) + term_v)];

    B = [zeros(n);
         s*eye(n)];

    S = [v;
         u + accl - c_d*norm(v)*v];

    w = f - A*x - B*u - S*s;

end