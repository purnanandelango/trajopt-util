clc
close all
clearvars

% Load problem data
load prob_data.mat

astro = plant.sclunar.astro_constants;

ybar = prb.yguess;
ubar = prb.uguess;

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
obj_val = 0.0;
cnstr = [];
for j = 1:prb.Np-1
    obj_val = obj_val + u(:,j)'*u(:,j);
    % obj_val = obj_val + norm(u(:,j));
    
    cnstr = [cnstr;
             y(:,j+1) == prb.A(:,:,j)*y(:,j) + prb.B*u(:,j)];

    
end
cnstr = [cnstr;
         y(:,1) == prb.y0;
         y(1:3,prb.Np) == prb.yend(1:3)];
optimize(cnstr,obj_val,sdpsettings('solver','osqp'));
ubar = value(u);
ybar = value(y);

figure
subplot(1,2,1)
ybar_dim = astro.S*ybar;
prb.AvoidSet_dim.plot
hold
plot3(ybar_dim(1,:),ybar_dim(2,:),ybar_dim(3,:),'.-b');
xlabel("$x$ [km]");
xlabel("$y$ [km]");
xlabel("$z$ [km]");
axis equal
ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
ax.DataAspectRatio = [1,1,1];
view(-150,26);

subplot(1,2,2)
ubar_dim = 100*astro.S(4:6,4:6)*prb.uscl*ubar; % [cm/s]
norm_ubar_dim(prb.Np-1) = 0;
for j = 1:prb.Np-1
    norm_ubar_dim(j) = norm(ubar_dim(:,j));
end
plot(prb.tp(1:prb.Np-1) *astro.nd2hr,norm_ubar_dim,'.-r');
ylabel("[cm/s]");
xlabel("[hr]");