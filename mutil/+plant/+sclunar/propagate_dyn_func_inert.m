% Integrate equations of motion in Moon-centered inertial frame with initial condition specified in Moon-centered inertial frame
function [x_sol_inert, t] = propagate_dyn_func_inert(x_init,tspan,astro)
    
    % Integrate with MATLAB ODE solver
    opts            = odeset('Reltol',1e-8,'AbsTol',1e-8);
    [t,x_sol_inert] = ode113(@(t,x) plant.sclunar.dyn_func(t,x,astro),tspan,x_init,opts);

end
