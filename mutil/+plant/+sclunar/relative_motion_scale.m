clc
clearvars
close all

% Load astrodynamics constants
astro = plant.sclunar.astro_constants();

% Load ephemeris data
cspice_furnsh( { 'naif0011.tls.pc',...
                 'de421.bsp',...
                 'pck00010.tpc' } );

hour_final = 6;
tspan = (0:300:hour_final*3600)/astro.t_star; % Grid with 5 min intervals

% Convert rotational to inertial 
x_init_inert = plant.sclunar.rot_to_inert(astro.nrho_init,tspan(1),astro);

[xbar,tbar] = plant.sclunar.propagate_dyn_func_inert(x_init_inert,tspan,astro,1,0);

% Constuct piecewise polynomial
ppbar = spline(tbar,xbar);

N = numel(tbar);
Abar = cell(1,N);
inf_norm(N) = 0;
for j = 1:N
    [~,Abar{j}] = plant.sclunar.dyn_func_inert_SRP_jac(tbar(j),xbar(:,j),astro);
    astro.invsclx*Abar{j}*astro.sclx
    inf_norm(j) = max(max(abs(astro.invsclx*Abar{j}*astro.sclx)));
end

plot(tbar*astro.t_star/3600,inf_norm); xlabel("hr");

rel_dist = 5; % [km]
rel_speed = 5; % [m/s]
y0 = astro.invsclx*[rel_dist*ones(3,1)/astro.l_star;rel_speed*ones(3,1)/(1000*astro.v_star)];
[~,y] = ode113(@(t,y) ltv_sys(t,y,ppbar,astro),tspan,y0);

y = (diag([astro.l_star*ones(1,3),1000*astro.v_star*ones(1,3)])*astro.sclx*(y'))';

figure
subplot(1,2,1)
plot3(y(:,1),y(:,2),y(:,3),'-r');
title('Relative position');

subplot(1,2,2)
plot3(y(:,4),y(:,5),y(:,6),'-b');
title('Relative velocity');

% Clear ephemeris data from memory
cspice_kclear

function dy = ltv_sys(t,y,ppbar,astro)
    xbar = ppval(ppbar,t);
    [~,A] = plant.sclunar.dyn_func_inert_SRP_jac(t,xbar,astro);
    dy = astro.invsclx*A*astro.sclx*y;
end
