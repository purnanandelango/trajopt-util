% Assemble backward reachable set

clc
clearvars
close all

% Load astrodynamics constants
astro = plant.sclunar.astro_constants();

% Load ephemeris data
cspice_furnsh( { 'naif0011.tls.pc',...
                 'de421.bsp',...
                 'pck00010.tpc' } );

% Non-dimensional initial condition in inertial frame
t0 = 0;
x0 = plant.sclunar.rot_to_inert(astro.nrho_init_rot,t0,astro);

delta_t = 5 *astro.min2nd; % Sampling time
 
tend   = 12  *astro.hr2nd; % Final time
ts_end = 6   *astro.hr2nd; % Safety horizon
tp_end = 4   *astro.hr2nd; % Planning horizon

% N = 1 + (tend-t0)/delta_t; % Grid size
% tspan = linspace(t0,tend,N);

Ns = 1 + round((ts_end-0)/delta_t);
Np = 1 + round((tp_end-t0)/delta_t);
ts = linspace(0,ts_end,Ns);
tp = linspace(t0,tp_end,Np);

ymax = astro.invS*[0.1*ones(3,1); 0.1*ones(3,1)]; % [km,km/s] -> [nd]
[H,h] = geom.construct_box(-ymax,ymax);

xbar = plant.sclunar.propagate_dyn_func_inert(x0,tp,astro,1,0);

BRS = cell(Ns,Np);

for j = 1:Np
    BRS{1,j} = H;
    for k = 2:Ns
        STM = disc.stm_ltv(tp(j),tp(j) + ts(k),xbar(:,j),@(t,x) plant.sclunar.dyn_func_inert_SRP(t,x,astro), ...
                                                         @(t,x) plant.sclunar.dyn_func_inert_SRP_jac(t,x,astro));
        BRS{k,j} = H*STM;
    end
end


% Clear ephemeris data from memory
cspice_kclear