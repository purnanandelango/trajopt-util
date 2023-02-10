% Integrate equations of motion in Moon-centered inertial frame with initial condition specified in Moon-centered inertial frame
function [x_sol_inert, t] = propagate_dyn_func_inert(x_init,tspan,astro,srp_flg,mJ2_flg)
    
    % Integrate with MATLAB ODE solver
    opts            = odeset('Reltol',1e-8,'AbsTol',1e-8);
    if srp_flg && mJ2_flg
        [t,x_sol_inert] = ode113(@(t,x) plant.sclunar.dyn_func_inert_SRP_J2(t,x,astro),tspan,x_init,opts);        
    elseif srp_flg
        [t,x_sol_inert] = ode113(@(t,x) plant.sclunar.dyn_func_inert_SRP(t,x,astro),tspan,x_init,opts);
    else
        [t,x_sol_inert] = ode113(@(t,x) plant.sclunar.dyn_func_inert(t,x,astro),tspan,x_init,opts);
    end

    x_sol_inert = x_sol_inert';

end
