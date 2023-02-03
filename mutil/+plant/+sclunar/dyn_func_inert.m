% Compute derivative of state according to equations of motion in N body problem in the J2000 inertial frame (centered at Moon) with SRP and Moon J2 considered
% Output: 
%   dx          = derivative of the state vector
% Input:
%   t           = non dimensional time 
%   x           = non dimensional state vector [x pos, y pos, z pos, x velo, y velo, z velo]
%   astro       = structure of astrodynamics constants
% 
% MIF == Moon-centered inertial frame
function dx = dyn_func_inert(t,x,astro)

    t_ephm = astro.start_JD + t*astro.t_star/(24*3600);
    
    dx = zeros(6,1);

    rsc = x(1:3);               % Normalized spacecraft position in MIF
    vsc = x(4:6);               % Normalized spacecraft position in MIF
    sc_dist = norm(rsc);        % Spacecraft distance to Moon
    
    dx(1:3) = vsc;    

    % Primary body gravitational force: Moon
    dx(4:6) = -astro.ndGM_Moon*rsc/sc_dist^3;
    
    % Perturbing body gravitatinal force: Earth
    ephm_state  = plant.sclunar.ephemeris_state(0,t_ephm,astro.Earth,astro.Moon);
    rE          = ephm_state(1:3)'/astro.l_star;            % Normalized Earth position in MIF
    rEsc        = rE - rsc;                                 % Normalized spacecraft to Earth vector in MIF
    dx(4:6)     = dx(4:6) + astro.ndGM_Earth*(rEsc/norm(rEsc)^3 - rE/norm(rE)^3);
    
    % Perturbing body gravitational force: Sun
    ephm_state  = plant.sclunar.ephemeris_state(0,t_ephm,astro.Sun,astro.Moon);
    rS          = ephm_state(1:3)'/astro.l_star;            % Normalized Sun position in MIF
    rSsc        = rS - rsc;                                 % Normalized spacecraft to Sun vector in MIF
    dx(4:6)     = dx(4:6) + astro.ndGM_Sun*(rSsc/norm(rSsc)^3 - rS/norm(rS)^3);

    % Solar radiation pressure
    a_SRP       = 0.5*astro.Cr_beta*astro.ndGM_Sun*rSsc/(norm(rSsc)^3);
    dx(4:6)     = dx(4:6) + a_SRP;

    % Moon J2
    ephm_state  = plant.sclunar.ephemeris_state(0,t_ephm,astro.Earth,astro.Moon); 
    perp_z_vec  = cross(ephm_state(1:3), ephm_state(4:6));
    perp_vec    = -cross(ephm_state(1:3), perp_z_vec);
    M_sc_vec    = rsc'*astro.l_star;
    M_sc        = norm(M_sc_vec);
    M_sc_perp   = (perp_vec/norm(perp_vec)) * dot(M_sc_vec(1:3),perp_vec)/norm(perp_vec);
    M_sc_plane  = M_sc_vec(1:3) - M_sc_perp;
    Angle1      = acosd(dot(M_sc_plane, ephm_state(1:3))/(norm(ephm_state(1:3))*norm(M_sc_plane))) + 6.68;
    lambda      = Angle1*pi/180;

    dx(4:6)      = dx(4:6) + (astro.ndGM_Moon*rsc/(sc_dist^3))*((3*astro.Moon_J2*(astro.R_Moon/M_sc)^2)*(1.5*(sin(lambda))^2 - 0.5));

end
