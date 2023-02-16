% Compute STM of the LTV approximation of a system about a known trajectory
function STM_t1t2 = stm_ltv(t1,t2,x0,dyn_func,dyn_func_jac)
    n = length(x0);
    z0 = [x0;reshape(eye(n),[n*n,1])];
    [~,z] = ode113(@(t,z) ode_fun(t,z,dyn_func,dyn_func_jac,n),[t1,t2],z0);
    STM_t1t2 = reshape(z(end,n+1:end)',[n,n]);
end
function dz = ode_fun(t,z,fun,dfun,n)
    A = dfun(t,z(1:n));
    Phi = reshape(z(n+1:end),[n,n]);
    dPhi = A*Phi; % Apply scaling to convert from non-dimensional states to rendezvous scale
    dz = [fun(t,z(1:n));dPhi(:)];
end