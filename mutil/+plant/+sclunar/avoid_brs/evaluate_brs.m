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
 
ts_end = 6   *astro.hr2nd; % Safety horizon
tp_end = 4   *astro.hr2nd; % Planning horizon

Ns = 1 + round((ts_end-0)/delta_t);
Np = 1 + round((tp_end-t0)/delta_t);
ts = linspace(0,ts_end,Ns);
tp = linspace(t0,tp_end,Np);

ymax_dim = [0.1*ones(3,1); 100*ones(3,1)]; % [km,m/s] 
[Hdim,hdim] = geom.construct_box(-ymax_dim,ymax_dim);

ymax = astro.invS*ymax_dim; % [km,m/s] -> [rdv]
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

uscl = 10;
umax = 10;
B = [zeros(3,3);eye(3)*uscl];

% Generate a relative motion trajectory
y0 = astro.invS*[1*ones(3,1); 0.1*ones(3,1)]; % [km,m/s] -> [rdv];
y = [y0,zeros(6,Np-1)];
for j = 1:Np-1
    y(:,j+1) = A(:,:,j)*y(:,j);
end

yend = astro.invS*[-0.12;0;0;0;0;0];

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

% [Hl,hl] = geom.project_polyhed2dims(H,h,1:3);
% AvoidSet = Polyhedron('A',Hl,'b',hl);
[Hl,hl] = geom.project_polyhed2dims(Hdim,hdim,1:3);
AvoidSet_dim = Polyhedron('A',Hl,'b',hl);


prb = struct;
prb.BRS = BRS;
prb.x0 = x0;
prb.t0 = t0;
prb.y0 = y0;
prb.yend = yend;
prb.yguess = grid.ends2interp(y0,yend,linspace(0,1,Np),'poly',1);
prb.uguess = zeros(3,Np-1);
prb.uscl = uscl;
prb.umax = umax;
prb.delta_t = delta_t;
prb.ts = ts;
prb.tp = tp;
prb.Ns = Ns;
prb.Np = Np;
prb.umax = umax;
prb.A = A;
prb.B = B;
prb.H = H;
prb.h = h;
% prb.AvoidSet = AvoidSet;
prb.AvoidSet_dim = AvoidSet_dim;

save('prob_data','prb');