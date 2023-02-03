% Transform Moon-centered rotating frame state to Moon-centered inertial frame
function [x_inert,T_rot_to_inert] = rot_to_inert(x_rot,t,astro)
    
    N = length(t);
    
    x_inert = zeros(1,6);
    for j=1:N
        ephm_state = plant.sclunar.ephemeris_state(t(j)*astro.t_star/(24*3600),astro.start_JD,astro.Moon,astro.Earth);
        pos = ephm_state(1,1:3)';
        vel = ephm_state(1,4:6)';

        l_star_inst = norm(pos);
        % t_star_inst = sqrt(l_star_inst^3/(astro.GM_Moon+astro.GM_Earth));
        % v_star_inst = l_star_inst/t_star_inst;

        H = cross(pos,vel);
        mag_H = norm(H);

        x_hat = pos/l_star_inst;
        z_hat = H/mag_H;
        y_hat = cross(z_hat,x_hat);

        Omega_dot = H/l_star_inst^2;

        x_rot_pos_dim   = x_rot(j,1:3)'*astro.l_star; % l_star_inst
        T_rot_to_inert  = [x_hat,y_hat,z_hat];

        x_inert_pos_dim = T_rot_to_inert*x_rot_pos_dim;

        x_inert_vel_dim = T_rot_to_inert*x_rot(j,4:6)'*astro.v_star + cross(Omega_dot,x_inert_pos_dim); % v_star_inst
        x_inert(j,:) = [x_inert_pos_dim/astro.l_star; x_inert_vel_dim/astro.v_star]';                   % l_star, v_star
        
    end

end
