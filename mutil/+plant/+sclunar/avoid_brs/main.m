clc
close all
clearvars

% Load problem data
prb = problem_data(0);

% Load constants
astro = plant.sclunar.astro_constants;

plant.sclunar.ephem('load');

% Initialize reference solution
ybar = prb.yguess;
ubar = prb.uguess;

for itr = 1:15

    % Compute signed-distance and projection to BRS
    sdist = zeros(prb.Ns,prb.Np);
    projvec = cell(prb.Ns,prb.Np);
    for j = 1:prb.Np
        for k = 1:prb.Ns
            [projvec{k,j},sdist(k,j)] = geom.sign_dist_polyhed(ybar(:,j),prb.BRS{k,j},prb.h);
        end
    end
    
    yalmip clear
    
    y = sdpvar(6,prb.Np);
    u = sdpvar(3,prb.Np-1);
    objval = 0.0;
    cnstr = [];

    objval_stg = 0.0;    
    for j = 1:prb.Np-1
        objval_stg = objval_stg + prb.stg_penalty*prb.cost_term(u(:,j));
        
        cnstr = [cnstr;
                 y(:,j+1) == prb.A(:,:,j)*y(:,j) + prb.B*u(:,j)]; 
    end
    cnstr = [cnstr;
             y(:,1) == prb.y0;
             y(1:3,prb.Np) == prb.yend(1:3)];

    objval = objval + objval_stg;

    objval_tr = 0.0;
    for j = 1:prb.Np-1
        objval_tr = objval_tr + prb.tr_penalty*prb.cost_term(y(:,j)-ybar(:,j)) ...
                              + prb.tr_penalty*prb.cost_term(u(:,j)-ubar(:,j));
    end
    objval_tr = objval_tr + prb.tr_penalty*prb.cost_term(y(:,prb.Np)-ybar(:,prb.Np));    
    
    objval = objval + objval_tr;

    slk = sdpvar(prb.Ns,prb.Np);
    dL = sdpvar(prb.Ns,prb.Np);
    objval_obs = 0.0;
    for j = 1:prb.Np
        for k = 1:prb.Ns
            dL(k,j) = sdist(k,j) + dot(ybar(:,j) - projvec{k,j},y(:,j)-ybar(:,j))/sdist(k,j) - prb.sdist_min;
            cnstr = [cnstr;
                     slk(k,j) >= 0;
                     -dL(k,j) <= slk(k,j)];    
            objval_obs = objval_obs + prb.obs_penalty*slk(k,j);
        end
    end
    objval = objval + objval_obs;
    
    optimize(cnstr,objval,sdpsettings('solver','osqp','verbose',0));
    
    ubar = value(u);
    ybar = value(y);

    fprintf("%2d TR  = %.2e\n",itr,value(objval_tr));
    fprintf("   OBS = %.2e\n\n",min(min(sdist)));

end

%% Compute failure trajectories
yfail_dim = zeros(6,prb.Ns,prb.Np);
for j = 1:prb.Np
    xfail = plant.sclunar.propagate_dyn_func_inert(prb.xbar(:,j) + astro.Srdv*ybar(:,j),prb.ts,astro,1,0);
    yfail_dim(:,:,j) = astro.Snd*(xfail - prb.xbar(:,j:j+prb.Ns-1));
end

plant.sclunar.ephem('unload');


%%
figure
% subplot(1,2,1)
ybar_dim = astro.S*ybar;
prb.AvoidSet_dim.plot('alpha',0.5)
hold on
plot3(ybar_dim(1,:),ybar_dim(2,:),ybar_dim(3,:),'o-b');
for j = 1:prb.Np
    plot3(yfail_dim(1,:,j),yfail_dim(2,:,j),yfail_dim(3,:,j),'.-','Color',[0.2,0.2,0.3,0.5]);
end
xlabel("$x$ [km]");
ylabel("$y$ [km]");
zlabel("$z$ [km]");
axis equal
ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
ax.DataAspectRatio = [1,1,1];
view(0,90);

figure
% subplot(1,2,2)
ubar_dim = 100*astro.S(4:6,4:6)*prb.uscl*ubar; % [cm/s]
norm_ubar_dim(prb.Np-1) = 0;
for j = 1:prb.Np-1
    norm_ubar_dim(j) = norm(ubar_dim(:,j));
end
plot(prb.tp(1:prb.Np-1) *astro.nd2hr,norm_ubar_dim,'.-r');
ylabel("[cm/s]");
xlabel("[hr]");
title("Control magnitude");