% Compute STM of the LTV approximation of a system about a known trajectory
function STM_t1t2 = stm_ltv(t1,t2,x0,dyn_func,dyn_func_jac,varargin)
    n = length(x0);
    z0 = [x0;reshape(eye(n),[n*n,1])];
    if nargin == 6 % If ODE solver and options are specified
        [~,z] = feval(varargin{1}{1},@(t,z) ode_fun(t,z,dyn_func,dyn_func_jac,n),[t1,t2],z0,varargin{1}{2});
    else
        [~,z] = ode45(@(t,z) ode_fun(t,z,dyn_func,dyn_func_jac,n),[t1,t2],z0,odeset('RelTol',1e-5,'AbsTol',1e-7));
    end
    STM_t1t2 = reshape(z(end,n+1:end)',[n,n]);
end
function dz = ode_fun(t,z,fun,dfun,n)
    A = dfun(t,z(1:n));
    Phi = reshape(z(n+1:end),[n,n]);
    dPhi = A*Phi; % Apply scaling to convert from non-dimensional states to rendezvous scale
    dz = [fun(t,z(1:n));dPhi(:)];
end