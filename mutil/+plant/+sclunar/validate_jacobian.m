% Validate jacobian of dyn_func_inert with centeral difference estimate

clc
clearvars
close all

% Load astrodynamics constants
astro = plant.sclunar.astro_constants();

% Load ephemeris data
plant.sclunar.ephem('load');

t0 = 0;
x0 = plant.sclunar.rot_to_inert(astro.nrho_init_rot,t0,astro);

day_final = 3;
tfinal = day_final*3600*24/astro.t_star;
[x,t] = plant.sclunar.propagate_dyn_func_inert(x0,[t0,tfinal],astro,1,0);

fun = @(z) plant.sclunar.dyn_func_inert_SRP(z(1),z(2:7),astro);

N = length(t);
err(N) = 0;
for j = 1:N

    % Analytic Jacobian
    [dft,dfx] = plant.sclunar.dyn_func_inert_SRP_jac(t(j),x(:,j),astro);
    df1 = [dft,dfx];

    % Finite difference Jacobian
    df2 = misc.num_jacobian(fun,[t(j);x(:,j)]);
    
    % Error in estimates
    err(j) = 100*norm(df1-df2)/norm(df1);

    norm(df1-df2);

end

% plot(t*astro.t_star/(24*3600),err,'-b'); xlabel("[day]")
% ylabel("Error \%")


% Compute STM

% via integration of continuous-time linearization
STM_v1 = disc.stm_ltv(t0,tfinal,x0,@(t,x) plant.sclunar.dyn_func_inert_SRP(t,x,astro), ...
                                   @(t,x) plant.sclunar.dyn_func_inert_SRP_jac(t,x,astro));


% via integration of continuous-time nonlinear dynamics with perturbed initial conditions
STM_v2 = eye(6);
pertval = 1e-4;
for j = 1:6
    pert = max(pertval,pertval*abs(x0(j))) * STM_v2(:,j);
    xp = plant.sclunar.propagate_dyn_func_inert(x0+pert,[t0,tfinal],astro,1,0); xp = xp(:,end);
    xm = plant.sclunar.propagate_dyn_func_inert(x0-pert,[t0,tfinal],astro,1,0); xm = xm(:,end);
    STM_v2(:,j) = 0.5*(xp-xm)/pert(j); 
end

norm(STM_v1(:)-STM_v2(:))

% Clear ephemeris data from memory
plant.sclunar.ephem('unload');
