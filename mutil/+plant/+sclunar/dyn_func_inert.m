% Compute derivative of state according to equations of motion in N body problem in the J2000 MEME inertial frame (centered at Moon) with SRP and Moon J2 considered
% Output: 
%   dx          = derivative of the state vector
% Input:
%   t           = non dimensional time 
%   x           = non dimensional state vector [x pos, y pos, z pos, x velo, y velo, z velo]
%   astro       = structure of astrodynamics constants
function dx = dyn_func_inert(t,x,astro)

    t_ephm = astro.start_JD + t*astro.t_star/(24*3600);
    
    dx = zeros(6,1);

    rsc = x(1:3);               % Normalized spacecraft position
    vsc = x(4:6);               % Normalized spacecraft position
    sc_dist = norm(rsc);        % Spacecraft distance to Moon
    
    dx(1:3) = vsc;    

    % Primary body gravitational force: Moon
    dx(4:6) = -astro.ndGM_Moon*rsc/sc_dist^3;
    
    % Perturbing body gravitatinal force: Earth
    ephm_state  = plant.sclunar.ephemeris_state(0,t_ephm,astro.Earth,astro.Moon);
    rE          = ephm_state(1:3)'/astro.l_star;            % Normalized Earth position
    vE          = ephm_state(4:6)'/astro.v_star;            % Normalized Earth velocity
    rEsc        = rE - rsc;                                 % Normalized spacecraft to Earth vector
    dx(4:6)     = dx(4:6) + astro.ndGM_Earth*(rEsc/norm(rEsc)^3 - rE/norm(rE)^3);
    
    % Perturbing body gravitational force: Sun
    ephm_state  = plant.sclunar.ephemeris_state(0,t_ephm,astro.Sun,astro.Moon);
    rS          = ephm_state(1:3)'/astro.l_star;            % Normalized Sun position in MIF
    rSsc        = rS - rsc;                                 % Normalized spacecraft to Sun vector in MIF
    dx(4:6)     = dx(4:6) + astro.ndGM_Sun*(rSsc/norm(rSsc)^3 - rS/norm(rS)^3);

    % Solar radiation pressure
    a_SRP       = astro.SRP*rSsc/(norm(rSsc)^3);
    dx(4:6)     = dx(4:6) - a_SRP;

    % Moon J2 zonal harmonic
    vEprp = -cross(rE,cross(rE,vE));
    rscEM = rsc - dot(rsc,vEprp)*vEprp/norm(vEprp)^2;
    lam_sc = acos( dot(rE,rscEM)/(norm(rE)*norm(rscEM)) ) + astro.Moon_eqincl;
    a_J2 = 1.5*astro.ndGM_Moon*astro.Moon_J2*((astro.R_Moon/astro.l_star)^2)*( 3*sin(lam_sc)^2 - 1 )*rsc/sc_dist^5;

    % perp_z_vec  = cross(rE, vE);
    % perp_vec    = -cross(rE, perp_z_vec);
    % M_sc_vec    = rsc;
    % M_sc        = norm(M_sc_vec);
    % M_sc_perp   = (perp_vec/norm(perp_vec)) * dot(M_sc_vec,perp_vec)/norm(perp_vec);
    % M_sc_plane  = M_sc_vec - M_sc_perp;
    % lam1        = acos(dot(M_sc_plane, rE)/(norm(rE)*norm(M_sc_plane))) + astro.Moon_eqincl;
    % a_J22       = (astro.ndGM_Moon*rsc/(sc_dist^3))*((3*astro.Moon_J2*(astro.R_Moon/(astro.l_star*M_sc))^2)*(1.5*(sin(lam1))^2 - 0.5))

    dx(4:6)      = dx(4:6) + a_J2;

end
