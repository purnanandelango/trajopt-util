% Integrate equations of motion in Moon-centered inertial frame with initial condition specified in Moon-centered rotating frame
function [x_sol_rot, t] = propagate_dyn_func_rot(x_init,tspan,astro,srp_flg,mJ2_flg)
    
    % Convert rotational to inertial J2000 equatorial 
    x_init_inert = plant.sclunar.rot_to_inert(x_init,tspan(1),astro);
    
    % Integrate with MATLAB ODE solver
    opts            = odeset('Reltol',1e-8,'AbsTol',1e-8);
    if srp_flg && mJ2_flg
        [t,x_sol_inert] = ode113(@(t,x) plant.sclunar.dyn_func_inert_SRP_J2(t,x,astro),tspan,x_init_inert,opts);        
    elseif srp_flg
        [t,x_sol_inert] = ode113(@(t,x) plant.sclunar.dyn_func_inert_SRP(t,x,astro),tspan,x_init_inert,opts);
    else
        [t,x_sol_inert] = ode113(@(t,x) plant.sclunar.dyn_func_inert(t,x,astro),tspan,x_init_inert,opts);
    end    
    % Note that the output of ode113 stacks the states along rows

    % Convert inertial to rotational frame 
    x_sol_rot = plant.sclunar.inert_to_rot(x_sol_inert',t,astro);

end
