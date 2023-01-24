% Assemble a structure of relevant astrodynamics constants
function astro = astro_constants()
    astro               = struct;
    
    astro.Sun           = 10;   % NAIF ID
    astro.Earth         = 399;  % NAIF ID
    astro.Moon          = 301;  % NAIF ID
    astro.GM_Sun        = 132712197035.766; % [km^3/s^2]
    astro.GM_Earth      = 398600.432896939; % [km^3/s^2]
    astro.GM_Moon       = 4902.80058214776; % [km^3/s^2]
    astro.muMoon        = astro.GM_Moon/(astro.GM_Earth+astro.GM_Moon);
    
    Cr          = 1.15;                                             % Coefficient of reflectivity                    
    W_star      = 1368;                                             % Solar irradiance [kg/s^3]
    c           = 299792458;                                        % Speed of light [m/s]
    P_star      = W_star/c;                                         % [N/m^2]
    MAR_star    = 2*P_star*(149597870.700)^2/astro.GM_Sun/1000;     % 2*P_star*L_star^2/GM_Sun [kg/m^2]
    MAR         = 15000/16;                                         % [kg/m^2]
    beta        = MAR_star/MAR;
    astro.Cr_beta       = Cr*beta;

    astro.l_star        = 385692.5;         % [km]
    astro.t_star        = 377084.152667039; % [s]
    astro.v_star        = astro.l_star/astro.t_star;    % [km/s]

    astro.start_JD      = 2459957.5; % Jan 13,2023 (Julian date) [day]

    % Nondimensional GM
    astro.ndGM_Sun      = astro.GM_Sun*(astro.t_star^2)/(astro.l_star^3);
    astro.ndGM_Earth    = astro.GM_Earth*(astro.t_star^2)/(astro.l_star^3);
    astro.ndGM_Moon     = astro.GM_Moon*(astro.t_star^2)/(astro.l_star^3);

    % Initial condition of a high-fidelity NRHO baseline in Moon-centered rotating frame
    astro.nrho_init     = [1.024839452754877  -0.001469610071726  -0.176331432119678  -0.002193318212420  -0.107740373416905  -0.013148385667988] ...
                          + [astro.muMoon-1,zeros(1,5)];

end