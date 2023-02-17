% Assemble backward reachable set

clc
clearvars
close all

% Load astrodynamics constants
astro = plant.sclunar.astro_constants();

% Load ephemeris data
plant.sclunar.ephem('load');

% Non-dimensional initial condition in inertial frame
t0 = 0;
x0 = astro.nrho_init_inert;

delta_t = 15 *astro.min2nd; % Sampling time
 
tend   = 12  *astro.hr2nd; % Final time
ts_end = 6   *astro.hr2nd; % Safety horizon
tp_end = 4   *astro.hr2nd; % Planning horizon

Ns = 1 + round((ts_end-0)/delta_t);
Np = 1 + round((tp_end-t0)/delta_t);
ts = linspace(0,ts_end,Ns);
tp = linspace(t0,tp_end,Np);

ymax = astro.invS*[0.1*ones(3,1); 100*ones(3,1)]; % [km,m/s] -> [rdv]
[H,h] = geom.construct_box(-ymax,ymax);

xbar = plant.sclunar.propagate_dyn_func_inert(x0,tp,astro,1,0);

BRS = cell(Ns,Np);
A = zeros(6,6,Np-1);

for j = 1:Np
    if j<Np
        A(:,:,j) = disc.stm_ltv(tp(j),tp(j+1),xbar(:,j),@(t,x) plant.sclunar.dyn_func_inert_SRP(t,x,astro), ...
                                                        @(t,x) plant.sclunar.dyn_func_inert_SRP_jac(t,x,astro,true));
    end
    BRS{1,j} = H;
    for k = 2:Ns
        STM = disc.stm_ltv(tp(j),tp(j) + ts(k),xbar(:,j),@(t,x) plant.sclunar.dyn_func_inert_SRP(t,x,astro), ...
                                                         @(t,x) plant.sclunar.dyn_func_inert_SRP_jac(t,x,astro,true));
        BRS{k,j} = H*STM;
    end
end

B = [zeros(3,3);eye(3)];

% Generate a relative motion trajectory
y0 = astro.invS*[5*ones(3,1); 5*ones(3,1)]; % [km,m/s] -> [rdv];
y = [y0,zeros(6,Np-1)];
for j = 1:Np-1
    y(:,j+1) = A(:,:,j)*y(:,j);
end

% Compute signed-distance and projection to BRS
sdist = zeros(Ns,Np);
projvec = cell(Ns,Np);
for j = 1:Np
    for k = 1:Ns
        [projvec{k,j},sdist(k,j)] = geom.sign_dist_polyhed(y(:,j),BRS{k,j},h);
    end
end

% Clear ephemeris data from memory
plant.sclunar.ephem('unload');