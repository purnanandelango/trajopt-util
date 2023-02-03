clearvars
close all
clc

astro = plant.sclunar.astro_constants();

t = sym('t');
x = sym('x_%d',[6,1]);
dx = sym('dx_%d',[6,1]);
x_M2E = sym('xM2E_%d',[6,1]);

t_ephm = astro.start_JD + t*astro.t_star/(24*3600);

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

