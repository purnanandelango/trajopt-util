function plot_soln(x,u,N,tvec,gamgs,thetmax,Vmax,umax,umin)
% visualization

scl_val = -0.01; % thrust vector scaling

figure
traj = plot3(x(1,:),x(2,:),x(3,:),'o-b','MarkerSize',7);
hold on
for i = 1:N-1
    thrust = quiver3(x(1,i),x(2,i),x(3,i),scl_val*u(1,i),scl_val*u(2,i),scl_val*u(3,i),'-r','LineWidth',2);
end
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
title('3D point-mass soft landing')

% plot glideslope cone
r = linspace(0,30,100);
th = linspace(0,2*pi,100);
[R,T] = meshgrid(r,th);
X = R.*cos(T);
Y = R.*sin(T);
Z = R/tan(gamgs);
cone = surf(X,Y,Z);
set(cone,'EdgeAlpha',0,'FaceColor',[0.7,0.5,1],'FaceAlpha',0.5)


% plot gimbal cone
r = linspace(0,1,100);
th = linspace(0,2*pi,100);
[R,T] = meshgrid(r,th);
for i = 1:3
    X = R.*cos(T) + x(1,i);
    Y = R.*sin(T) + x(2,i);
    Z = -R/tan(thetmax) + x(3,i);
    cone2 = surf(X,Y,Z);
    set(cone2,'EdgeAlpha',0,'FaceColor',[0.5,0.7,1],'FaceAlpha',0.5)
end

legend([traj,thrust,cone,cone2],{'Trajectory','Thrust vector','Glide slope cone','Gimbal cone'},'Location',[0.3,0.3,0.03,0.02])

axis square
ax = gca;
ax.DataAspectRatio = [1,1,1];
ax.PlotBoxAspectRatio = [1,1,1];
view(-71,14)

figure
subplot(1,2,1)
plot(tvec,Vmax*ones(1,N),'--k','LineWidth',2);
hold on
plot(tvec,sqrt(x(4,:).^2 + x(5,:).^2 + x(6,:).^2),'o-r')
ylim([0,1.1*Vmax])
xlabel('$t$');
title('Velocity magnitude');
legend('$V_{\max}$')

subplot(1,2,2)
hold on
p_umax = plot(tvec,umax*ones(1,N),'--k','LineWidth',2);
p_umin = plot(tvec,umin*ones(1,N),'--r','LineWidth',2);
stairs(tvec(1:size(u,2)),sqrt(u(1,:).^2 + u(2,:).^2 + u(3,:).^2),'o-b','MarkerSize',10)
xlabel('$t$');
title('Thrust magnitude $\|u\|_2$');
legend([p_umax,p_umin],{'$u_{\max}$','$u_{\min}$'})

end
