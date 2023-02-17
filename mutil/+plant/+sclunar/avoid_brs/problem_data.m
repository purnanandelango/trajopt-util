function prb = problem_data(brs_flg)
% Load astrodynamics constants
astro = plant.sclunar.astro_constants();

% Load ephemeris data
plant.sclunar.ephem('load');

% Non-dimensional initial condition in inertial frame
t0 = 0;
x0 = astro.nrho_init_inert;

delta_t = 15 *astro.min2nd; % Sampling time
 
ts_end = 6   *astro.hr2nd; % Safety horizon
tp_end = 6   *astro.hr2nd; % Planning horizon

Ns = 1 + round((ts_end-0)/delta_t);
Np = 1 + round((tp_end-t0)/delta_t);
ts = linspace(0,ts_end,Ns);
tp = linspace(t0,tp_end,Np);

xbar = plant.sclunar.propagate_dyn_func_inert(x0,[tp,tp(end)+ts(2:end)],astro,1,0);

% Construct avoid set around the target
box_size = 0.1;
ymax_dim = [box_size*ones(3,1); 100*ones(3,1)]; % [km,m/s] 
[Hdim,hdim] = geom.construct_box(-ymax_dim,ymax_dim);

ymax = astro.invS*ymax_dim; % [km,m/s] -> [rdv]
[H,h] = geom.construct_box(-ymax,ymax);

if brs_flg
    % Compute BRS
    [BRS,A] = plant.sclunar.compute_brs_grid(x0,tp,ts,H,astro);
else
    load('brs_data');
end

uscl = 0.1;
umax = 10;
B = [zeros(3,3);eye(3)*uscl];

y0 = astro.invS*[1; 0; 0; 0.1*ones(3,1)]; % [km,m/s] -> [rdv];
yguess = [y0,zeros(6,Np-1)];
for j = 1:Np-1
    yguess(:,j+1) = A(:,:,j)*yguess(:,j);
end
yend = astro.invS*[0.12;0;0;0;0;0];

% Clear ephemeris data from memory
plant.sclunar.ephem('unload');

[Hpos_dim,hpos_dim] = geom.project_polyhed2dims(Hdim,hdim,1:3);
AvoidSet_dim = Polyhedron('A',Hpos_dim,'b',hpos_dim);

prb = struct;
prb.BRS = BRS;
prb.xbar = xbar;
prb.x0 = x0;
prb.t0 = t0;
prb.y0 = y0;
prb.yend = yend;
prb.yguess = yguess;
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
prb.AvoidSet_dim = AvoidSet_dim;

prb.sdist_min = 0.01;
prb.stg_penalty = 0.5;
prb.tr_penalty = 1;
prb.obs_penalty = 20;

prb.cost_term = @(z) z'*z;      % Quadratic
% prb.cost_term = @(z) norm(z);   % Second-order cone

save('problem_data','prb');

end