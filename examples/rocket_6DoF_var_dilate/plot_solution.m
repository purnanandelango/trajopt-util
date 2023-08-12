clearvars
close all

load recent_solution

figure
% (mplot,nplot)
mplot = 3;
nplot = 4;

% Trajectory, vehicle axis, thrust vector and drag
subplot(mplot,nplot,[1,2,5,6])
n = 5;
plant.rocket6DoF.plot_vehicle_forces(u(1:3,1:n:end),rI(:,1:n:end),vI(:,1:n:end),qBI(:,1:n:end),0.4,0.2,struct('scl',0.4,'rho',prb.rho,'SA',prb.SA,'CA',prb.CA),{[2,3,1],{'y','z','x'},'x'});
plot3(xbar(3,:),xbar(4,:),xbar(2,:),'o','Color',[0.7,0.1,0.7],'LineWidth',0.5,'DisplayName','SCP solution');
legend('Position',[0.279,0.442,0.236,0.167])
grid on

% Thrust
subplot(mplot,nplot,[3,4])
plt.plot_vec_nrm(tvec,u(1:3,:),[prb.Tmin,prb.Tmax],2,'$\|T_{\mathcal{B}}(t)\|_2$','linear','blue');
nrm_Tbar = misc.compute_vec_norm(ubar(1:3,:));
hold on 
plot(tvecbar,nrm_Tbar,'ob');

nrm_vI = misc.compute_vec_norm(vI);
[~,idx] = find(abs(nrm_vI-prb.Vmax_STC)<1e-1);
t_trig = tvec(min(idx));

% Speed
subplot(mplot,nplot,[7,8])
plt.plot_vec_nrm(tvec,vI,prb.Vmax,2,'$\|v_{\mathcal{I}}(t)\|_2$','linear','blue',"Single shot");
nrm_vIbar = misc.compute_vec_norm(xbar(5:7,:));
hold on 
legend('AutoUpdate','on');
plot(tvecbar,nrm_vIbar,'ob','DisplayName','SCP');
legend('AutoUpdate','off','Position',[0.773,0.507,0.122,0.075]);
plot(t_trig*ones(1,100),linspace(0,prb.Vmax),'-k');
plot(tvec,prb.Vmax_STC*ones(1,prb.Kfine),'-k');

% Dilation factors
subplot(mplot,nplot,[9,10])
plot(tau,u(4,:),'-b','DisplayName','Dilation Factors');
hold on
plot(tau,tvec,'-r','DisplayName','Time');
legend('Location','best','AutoUpdate','off');
plot(taubar,ubar(4,:),'ob');
plot(taubar,tvecbar,'or');
xlim([0,1]);

% Angle of attack
AoA = misc.compute_vec_norm(plant.rocket6DoF.compute_AoA(vI,qBI));
AoAbar = misc.compute_vec_norm(plant.rocket6DoF.compute_AoA(xbar(5:7,:),xbar(8:11,:)));

subplot(mplot,nplot,[11,12])
plot(tvec,AoA,'-b');
hold on
plot(tvecbar,AoAbar,'ob');
title("Angle of attack [${}^\circ$]");

plot(t_trig*ones(1,100),linspace(0.9*min(AoA),1.1*max(AoA)),'-k');
plot(tvec,acosd(prb.cosaoamax)*ones(1,prb.Kfine),'-k');
xlim([0,tvec(end)]);
ylim([0.9*min(AoA),1.1*max(AoA)]);

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