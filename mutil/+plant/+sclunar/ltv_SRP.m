% LTV approximation of the system about a known trajectory
% represented by a piece-wise polynomial
function dy = ltv_SRP(t,y,ppbar,astro)
    xbar = ppval(ppbar,t); % Call piece-wise polynomial interpolant
    [~,A] = plant.sclunar.dyn_func_inert_SRP_jac(t,xbar,astro);
    dy = astro.invSrdv*A*astro.Srdv*y; % Apply scaling to convert from non-dimensional states to rendezvous scale
end