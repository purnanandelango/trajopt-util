% visualization

scl_val = -0.01; % thrust vector scaling

figure
traj = plot3(x(1,:),x(2,:),x(3,:),'o-b','MarkerSize',7);
hold on
Nlim = N-1;
for i = 1:Nlim
    thrust = quiver3(x(1,i),x(2,i),x(3,i),scl_val*T(1,i),scl_val*T(2,i),scl_val*T(3,i),'-r','LineWidth',2);
end
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
title('3D point-mass soft landing')

% plot glideslope cone
r_gs = linspace(0,5*max(abs(x0(1:3))),100);
th_gs = linspace(0,2*pi,100);
[R_gs,T_gs] = meshgrid(r_gs,th_gs);
X = R_gs.*cos(T_gs);
Y = R_gs.*sin(T_gs);
Z = R_gs/tangamgs;
cone = surf(X,Y,Z);
set(cone,'EdgeAlpha',0,'FaceColor',[0.7,0.5,1],'FaceAlpha',0.5)


% % plot gimbal cone
% r_pnt = linspace(0,200,100);
% th_pnt = linspace(0,2*pi,100);
% [R_pnt,T_pnt] = meshgrid(r_pnt,th_pnt);
% for i = 1:3
%     X = R_pnt.*cos(T_pnt) + x(1,i);
%     Y = R_pnt.*sin(T_pnt) + x(2,i);
%     Z = -R_pnt/tan(theta_pnt) + x(3,i);
%     cone2 = surf(X,Y,Z);
%     set(cone2,'EdgeAlpha',0,'FaceColor',[0.5,0.7,1],'FaceAlpha',0.5)
% end

legend([traj,thrust,cone],{'Trajectory','Thrust vector','Glide slope cone'},...
    ...% 'Location',[0.5,0.5,0.03,0.02])
    'Location','northeast')

axis square
ylim([-100,2.5*abs(x0(2))])
zlim([-100,2.5*abs(x0(3))])
ax = gca;
ax.DataAspectRatio = [1,1,1];
ax.PlotBoxAspectRatio = [1,1,1];
view(90,0)

figure
subplot(1,2,1)
plot(tvec,Vmax*ones(1,N),'--k','LineWidth',2);
hold on
plot(tvec,sqrt(x(4,:).^2 + x(5,:).^2 + x(6,:).^2),'o-r')
ylim([0,1.1*Vmax])
xlim([0,tvec(end)])
xlabel('$t$');
title('Velocity magnitude');
legend('$V_{\max}$')

subplot(1,2,2)
hold on
p_Tmax = plot(tvec(1:Nlim),rho1*ones(1,Nlim),'--k','LineWidth',2);
p_Tmin = plot(tvec(1:Nlim),rho2*ones(1,Nlim),'--r','LineWidth',2);
plot(tvec(1:Nlim),lcvx_test,'sm','LineWidth',1);
stairs(tvec(1:Nlim),sqrt(T(1,:).^2 + T(2,:).^2 + T(3,:).^2),'o-b','MarkerSize',10)
xlabel('$t$');
title('Thrust magnitude $\|T\|_2$');
legend([p_Tmax,p_Tmin],{'$T_{\max}$','$T_{\min}$'})