% Compute STM of the LTV approximation of the about a known trajectory
function STM_t1t2 = stm_ltv_SRP(t1,t2,ppbar,astro)
    z0 = reshape(eye(6),[36,1]);
    [~,z] = ode113(@(t,z) ode_fun(t,z,ppbar,astro),[t1,t2],z0);
    STM_t1t2 = reshape(z(end,:)',[6,6]);
end
function dz = ode_fun(t,z,ppbar,astro)
    xbar = ppval(ppbar,t); % Call piece-wise polynomial interpolant
    [~,A] = plant.sclunar.dyn_func_inert_SRP_jac(t,xbar,astro);    
    Phi = reshape(z,[6,6]);
    dPhi = astro.invSrdv*A*astro.Srdv*Phi; % Apply scaling to convert from non-dimensional states to rendezvous scale
    dz = dPhi(:);
end