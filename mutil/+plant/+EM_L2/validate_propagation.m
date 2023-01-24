% Validate the propagation functions against the result of existing C++ code
clc
clearvars
close all

addpath('../../../../../OpNav/util');
addpath('../../../../../OpNav/util_mex');

% Load astrodynamics constants
astro = plant.EM_L2.astro_constants();

% Load ephemeris data
cspice_furnsh( { 'naif0011.tls.pc',...
                 'de421.bsp',...
                 'pck00010.tpc' } );

tspan = [0,14*24*3600/astro.t_star];

tic
[x2,t2] = plant.EM_L2.propagate_dyn_func_rot(astro.nrho_init,tspan,astro);
x2 = x2 + repmat([1-astro.muMoon,zeros(1,5)],[length(t2),1]);
fprintf("MATLAB propagator time: %5.3d ms\n",1000*toc);

tic
IC = astro.nrho_init + [1-astro.muMoon,zeros(1,5)];
[x1,t1] = f_Integrate_SRP_J2(IC,tspan,astro.muMoon,astro.start_JD,[astro.Earth,astro.Moon,astro.Sun],...
                             [astro.GM_Earth,astro.GM_Moon,astro.GM_Sun],astro.l_star,astro.t_star,false,astro.Cr_beta,astro.GM_Sun,'mex');
fprintf("MEX propagator time: %5.3d ms\n",1000*toc);

% Clear ephemeris data from memory
cspice_kclear

plot3(x2(:,1),x2(:,2),x2(:,3),'-b');
hold on
plot3(x1(:,1),x1(:,2),x1(:,3),'--r','LineWidth',3);
axis equal
ax = gca;
ax.PlotBoxAspectRatioMode = 'manual';
ax.PlotBoxAspectRatio = [1,1,1];
ax.DataAspectRatioMode = 'manual';
ax.DataAspectRatio = [1,1,1];
ax.View = [42,22];

rmpath('../../../../../OpNav/util');
rmpath('../../../../../OpNav/util_mex');
