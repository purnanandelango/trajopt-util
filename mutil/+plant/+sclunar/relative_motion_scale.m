% Examine the scale of relative motion on NRHO over a span of few days
clc
clearvars
close all

% Load astrodynamics constants
astro = plant.sclunar.astro_constants();

% Load ephemeris data
plant.sclunar.ephem('load');

hour_final = 72; % [hr]
tspan = (0:300:hour_final*3600)/astro.t_star; % Grid with 5 min intervals

% Convert nd initial condition from rotational to inertial frame
x_init_inert = astro.nrho_init_inert;

[xbar,tbar] = plant.sclunar.propagate_dyn_func_inert(x_init_inert,tspan,astro,1,0);

% Constuct piecewise polynomial
ppbar = spline(tbar,xbar);

N = numel(tbar);
Abar = cell(1,N);
inf_norm(N) = 0;
for j = 1:N
    Abar{j} = plant.sclunar.dyn_func_inert_SRP_jac(tbar(j),xbar(:,j),astro,true);
    inf_norm(j) = max(max(abs(Abar{j})));
end

plot(tbar*astro.t_star/3600,inf_norm); xlabel("hr"); title("RDV scale")

rel_dist = 5; % [km]
rel_speed = 5; % [m/s]
unit_vec1 = randn(3,1); unit_vec1 = unit_vec1/norm(unit_vec1);
unit_vec2 = randn(3,1); unit_vec2 = unit_vec2/norm(unit_vec2);
y0 = astro.invS*[rel_dist*unit_vec1;rel_speed*unit_vec2];
[~,y] = ode113(@(t,y) plant.sclunar.ltv_SRP(t,y,ppbar,astro,true),tspan,y0);

% y = (astro.S*(y'))';

figure
subplot(1,2,1)
plot3(y(:,1),y(:,2),y(:,3),'-r');
title('Relative position');
grid on
axis equal

subplot(1,2,2)
plot3(y(:,4),y(:,5),y(:,6),'-b');
title('Relative velocity');
grid on
axis equal


% Clear ephemeris data from memory
plant.sclunar.ephem('unload');
