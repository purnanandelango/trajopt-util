% LTV approximation of the system about a known trajectory represented by a piece-wise polynomial
function dy = ltv_SRP(t,y,ppbar,astro,scl_flg)
    xbar = ppval(ppbar,t); % Call piece-wise polynomial interpolant
    % Apply scaling to convert from non-dimensional states to rendezvous scale if scl_flg is true
    A = plant.sclunar.dyn_func_inert_SRP_jac(t,xbar,astro,scl_flg);
    dy = A*y; 
end