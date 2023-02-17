% Assemble a structure of astrodynamics constants relevant in cis-lunar space
function astro = astro_constants()
    astro               = struct;
    
    astro.Sun           = 10;   % NAIF ID
    astro.Earth         = 399;  % NAIF ID
    astro.Moon          = 301;  % NAIF ID
    astro.GM_Sun        = 132712197035.766; % [km^3/s^2]
    astro.GM_Earth      = 398600.432896939; % [km^3/s^2]
    astro.GM_Moon       = 4902.80058214776; % [km^3/s^2]
    astro.muMoon        = astro.GM_Moon/(astro.GM_Earth+astro.GM_Moon);
    
    astro.R_Moon        = 1737.1;       % Moon radius [km]
    astro.Moon_J2       = 2.024e-4;     % Moon zonal J2 value: http://web.gps.caltech.edu/classes/ge131/notes2016/Ch11.pdf
    astro.Moon_eqincl   = 6.68*pi/180;  % Moon equitorial inclination [rad]

    astro.l_star        = 385692.5;         % [km]
    astro.t_star        = 377084.152667039; % [s]
    astro.v_star        = astro.l_star/astro.t_star;    % [km/s]

    ksc         = 1.15;                                             % Coefficient of reflectivity                    
    S0          = 1360;                                             % Solar irradiance at 1 AU [W/m^2 or kg/s^3]
    r0          = 1495978707;                                       % 1 AU [km] 
    c           = 299792.458;                                       % Speed of light [km/s]
    Msc         = 15000;                                            % Mass of spacecraft [kg]
    Asc         = 16*1e-6;                                          % Cross-sectional area of spacecraft [km^2] 
    S0bycMsc    = S0/(Msc*c);                                       % [km^-1 s^-2]
    astro.SRP   = ksc*Asc*S0bycMsc*r0^2;                            % [km^3/s^2]
    astro.SRP   = astro.SRP*(astro.t_star^2)/(astro.l_star^3);      % [ndim]

    astro.start_JD      = 2459957.5; % Jan 13,2023 (Julian date) [day]

    % Nondimensional GM
    astro.ndGM_Sun      = astro.GM_Sun*(astro.t_star^2)/(astro.l_star^3);
    astro.ndGM_Earth    = astro.GM_Earth*(astro.t_star^2)/(astro.l_star^3);
    astro.ndGM_Moon     = astro.GM_Moon*(astro.t_star^2)/(astro.l_star^3);

    % Non-dimensional initial condition of a high-fidelity NRHO baseline in Moon-centered frame
    % Rotating frame
    astro.nrho_init_rot = [1.024839452754877  -0.001469610071726  -0.176331432119678  -0.002193318212420  -0.107740373416905  -0.013148385667988]' ...
                          + [astro.muMoon-1;zeros(5,1)];                                
    % Inertial frame
    plant.sclunar.ephem('load');
    astro.nrho_init_inert = plant.sclunar.rot_to_inert(astro.nrho_init_rot,0,astro);
    plant.sclunar.ephem('unload');

    % Scaling matrices for converting between rendezvous scale and EM CR3BP non-dimensionalization
    astro.Srdv = diag([1e-5*ones(1,3),5e-4*ones(1,3)]); % [rdv] -> [nd] 
    astro.invSrdv = inv(astro.Srdv);                    % [nd] -> [rdv]

    % Scaling matrices for converting between EM CR3BP non-dimensionalization and dimensional quantities
    astro.Snd = diag([astro.l_star*ones(1,3),astro.v_star*1000*ones(1,3)]); % [nd] -> [km,m/s]
    astro.invSnd = inv(astro.Snd);                                          % [km,m/s] -> [nd]

    % Scaling matrices for converting between rendezvous scale and dimensional quantities
    astro.S = astro.Snd*astro.Srdv; % [rdv] -> [km,m/s]
    astro.invS = inv(astro.S);      % [km,m/s] -> [rdv]

    % Unit conversions
    astro.sec2nd = 1/astro.t_star;   % [s] -> [nd]
    astro.nd2sec = 1/astro.sec2nd;   % [nd] -> [s]    

    astro.min2nd = 60/astro.t_star;  % [min] -> [nd]
    astro.nd2min = 1/astro.min2nd;   % [nd] -> [min]        

    astro.hr2nd = 3600/astro.t_star; % [hr] -> [nd]
    astro.nd2hr = 1/astro.hr2nd;     % [nd] -> [hr]

    astro.day2nd = 24*3600/astro.t_star; % [day] -> [nd]
    astro.nd2day = 1/astro.day2nd;       % [nd] -> [day]

end
