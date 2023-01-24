% Integrate equations of motion in Moon-centered inertial frame with initial condition specified in Moon-centered rotating frame
function [x_sol_rot, t] = propagate_dyn_func_rot(x_init,tspan,astro)
    
    % Convert rotational to inertial J2000 equatorial 
    x_init_inert = EM_L2.rot_to_inert(x_init,tspan(1),astro);
    
    % Integrate with MATLAB ODE solver
    opts            = odeset('Reltol',1e-8,'AbsTol',1e-8);
    [t,x_sol_inert] = ode45(@(t,x) EM_L2.dyn_func_inert(t,x,astro),tspan,x_init_inert,opts);

    % Convert inertial to rotational frame 
    x_sol_rot = EM_L2.inert_to_rot(x_sol_inert,t,astro);

end