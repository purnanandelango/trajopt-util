% Transform Moon-centered inertial frame state to Moon-centered rotating frame
function [x_rot,T_inert_to_rot] = inert_to_rot(x_inert,t,astro)

    N = length(t);
    
    x_rot = zeros(N,6);
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

        x_inert_pos_dim  = x_inert(j,1:3)'*astro.l_star; % l_star
        T_inert_to_rot   = [x_hat,y_hat,z_hat]';
        x_inert_vel_dim  = x_inert(j,4:6)'*astro.v_star; % (l_star/t_star)
        
        x_rot_vel_dim_i  = x_inert_vel_dim - cross(Omega_dot,x_inert_pos_dim);
        
        x_rot1  = [T_inert_to_rot zeros(3,3); zeros(3,3) T_inert_to_rot]*[x_inert_pos_dim;x_rot_vel_dim_i]; 
        x_rot(j,:) = [x_rot1(1:3,1)'/astro.l_star, x_rot1(4:6,1)'/astro.v_star]; % l_star_inst, v_star_inst
        
    end

end
