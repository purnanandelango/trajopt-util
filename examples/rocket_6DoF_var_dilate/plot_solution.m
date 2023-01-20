clearvars
close all

load recent_solution

figure
% (mplot,nplot)
mplot = 3;
nplot = 4;

% Trajectory, vehicle axis, thrust vector and drag
subplot(mplot,nplot,[1,2,5,6])
plant.rocket6DoF.plot_vehicle_forces(u(1:3,:),rI,vI,qBI,0.4,0.2,struct('scl',0.4,'rho',prb.rho,'SA',prb.SA,'CA',prb.CA),{[2,3,1],{'y','z','x'},'x'});
plot3(xbar(3,:),xbar(4,:),xbar(2,:),'--','Color',[0.8,0.8,0.2],'LineWidth',2,'DisplayName','SCP solution');

% Thrust
subplot(mplot,nplot,[3,4])
plt.plot_vec_nrm(tvec,u(1:3,:),[prb.Tmin,prb.Tmax],2,'Single shot: $\|T_{\mathcal{B}}(t)\|_2$');
plt.plot_vec_nrm(tvecbar,ubar(1:3,:),[prb.Tmin,prb.Tmax],2,'SCP: $\|T_{\mathcal{B}}(t)\|_2$','linear',[0.3,0.4,0.5]);

% Speed
subplot(mplot,nplot,[7,8])
plt.plot_vec_nrm(tvec,vI,prb.Vmax,2,'Single shot: $\|v_{\mathcal{I}}(t)\|_2$');
plt.plot_vec_nrm(tvecbar,xbar(5:7,:),prb.Vmax,2,'SCP: $\|v_{\mathcal{I}}(t)\|_2$','linear',[0.3,0.4,0.5]);

% Dilation factors
subplot(mplot,nplot,[9,10])
plot(1:prb.Kfine,u(4,:),'-b','DisplayName','Dilation Factors');
hold on
plot(1:prb.Kfine,tvec,'-r','DisplayName','Time');
legend('Location','best');

%%% Diagnostics %%%
% figure
% subplot(2,2,1)
% plot3(rI(1,:),rI(2,:),rI(3,:),'-r');
% hold on
% plot3(x2(2,:),x2(3,:),x2(4,:),'--b');
% title("Position");
% axis equal
% 
% subplot(2,2,2)
% plot3(vI(1,:),vI(2,:),vI(3,:),'-r');
% hold on
% plot3(x2(5,:),x2(6,:),x2(7,:),'--b');
% title("Velocity");
% axis equal
% 
% subplot(2,2,3)
% plot3(qBI(1,:),qBI(2,:),qBI(3,:),'-r');
% hold on
% plot3(x2(8,:),x2(9,:),x2(10,:),'--b');
% title("Quaternion [1:3]");
% axis equal
%%%%%%