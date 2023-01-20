clearvars
close all

load recent_solution

figure
% (mplot,nplot)
mplot = 2;
nplot = 4;

% Trajectory, vehicle axis, thrust vector and drag
subplot(mplot,nplot,[1,2,5,6])
plant.rocket6DoF.plot_vehicle_forces(u,rI,vI,qBI,0.4,0.2,struct('scl',0.4,'rho',prb.rho,'SA',prb.SA,'CA',prb.CA),{[2,3,1],{'y','z','x'},'x'});

% Thrust
subplot(mplot,nplot,[3,4])
plt.plot_vec_nrm(tvec,u,[prb.Tmin,prb.Tmax],2,'$\|T_{\mathcal{B}}(t)\|_2$');

% Speed
subplot(mplot,nplot,[7,8])
plt.plot_vec_nrm(tvec,vI,prb.Vmax,2,'$\|v_{\mathcal{I}}(t)\|_2$');
