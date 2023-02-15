% Validate jacobian of dyn_func_inert with centeral difference estimate

clc
clearvars
close all

% Load astrodynamics constants
astro = plant.sclunar.astro_constants();

% Load ephemeris data
cspice_furnsh( { 'naif0011.tls.pc',...
                 'de421.bsp',...
                 'pck00010.tpc' } );

t0 = 0;
x0 = plant.sclunar.rot_to_inert(astro.nrho_init_rot,t0,astro);

day_final = 15;
[x,t] = plant.sclunar.propagate_dyn_func_inert(x0,[t0,day_final*3600*24/astro.t_star],astro,1,0);

fun = @(z) plant.sclunar.dyn_func_inert_SRP(z(1),z(2:7),astro);

N = length(t);
err(N) = 0;
for j = 1:N

    [dft,dfx] = plant.sclunar.dyn_func_inert_SRP_jac(t(j),x(:,j),astro);
    df1 = [dft,dfx];
    df2 = misc.num_jacobian(fun,[t(j);x(:,j)]);
    
    err(j) = 100*norm(df1-df2)/norm(df1);

    norm(df1-df2);

end

plot(t*astro.t_star/(24*3600),err,'-b');xlabel("[day]")
ylabel("Error \%")

% Clear ephemeris data from memory
cspice_kclear
