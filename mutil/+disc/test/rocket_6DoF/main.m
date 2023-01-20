clearvars
close all
clc

foh_type = "v3"; % v1, v2, v3


% Plot figures
flg0 = true;

% Test rk4 march?
flg = true;

% State vector : x = [m;rI;vI;q;omegaB]
% Control input : u = TB

% System parameters
g0 = 9.81;
c_ax = 0.5; c_ayz = 1; CA = diag([c_ax,c_ayz,c_ayz]);
JB = 0.168*diag([2e-2,1,1]);
gI = [-g0;0;0];
rho = 1.225;
SA = 0.5;
rTB = -0.25*[1;0;0];
rcpB = 0.05*[1;0;0];
alphmdt = 1/(30*g0);
betmdt = 0;

% Convenient functions for accessing RHS of nonlinear and linearized ODE
sixDoF_rocket_dyn = @(t,x,u,s) plant.rocket6DoF.dyn_func(x,u,s,c_ax,c_ayz,diag(JB),gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
sixDoF_rocket_dyn_linearize = @(tbar,xbar,ubar,sbar) plant.rocket6DoF.compute_linearization(xbar,ubar,sbar,c_ax,c_ayz,diag(JB),gI,rho,SA,rTB,rcpB,alphmdt,betmdt);

nx = 14;
nu = 3;
tspan = [0,1]; % Dilated time interval

% Baseline solution
x0 = [2;[0;0;0];[100;10;0];[1;0;0;0];[0;0;0]];  % Initial condition
ubarfunc = @(t) [1./(t+0.1);exp(-t.^2);t.*0];   % Control input profile
sbar = 3;

% Compute baseline solution of nonlinear ODE
h = 0.01*diff(tspan);
tbar_ = tspan(1):h:tspan(2);
ubar_ = ubarfunc(tbar_);
[~,xbar_] = disc.rk4_march(sixDoF_rocket_dyn,tspan,x0,h,ubarfunc,sbar);
if flg
    [~,xbar2_] = ode45(@(t,x) sixDoF_rocket_dyn(t,x,ubarfunc(t),sbar),tbar_,x0);
end

h = 0.01*diff(tbar_(1:2)); % Step size for integration that computes FOH matrices

fprintf('Computation of FOH matrices:\n')
tic
[Ak,Bmk,Bpk,Sk,wk] = feval("disc.compute_foh_"+foh_type,tbar_,xbar_,ubar_,sbar,h,sixDoF_rocket_dyn,sixDoF_rocket_dyn_linearize);
toc

N = size(xbar_,2);
z_ = zeros(nx,N); z_(:,1) = x0;
v_ = ubarfunc(tbar_);
for k = 1:N-1
    z_(:,k+1) = Ak(:,:,k)*z_(:,k) + Bmk(:,:,k)*v_(:,k) + Bpk(:,:,k)*v_(:,k+1) + Sk(:,k)*sbar + wk(:,k);
end

errVal = abs(z_-xbar_);
fprintf("\nError in FOH and xbar:\n%.2e\n",max(max(errVal)));

%%

if flg0
    if flg
    % Validation of results for ode45 and rk4_march
        figure
        plot3(xbar2_(:,3),xbar2_(:,4),xbar2_(:,2),'--k')
        hold on
        plot3(xbar_(3,:),xbar_(4,:),xbar_(2,:),'-.r')
        axis equal
        title('Comparison of \texttt{ode45} and \texttt{rk4\_march}');
    end

    % Comparison of xbar and foh solution
    figure
    subplot(2,3,[1,4])
    plot3(xbar_(3,:),xbar_(4,:),xbar_(2,:),'-k')
    hold on
    plot3(z_(3,:),z_(4,:),z_(2,:),'--r')
    axis equal
    xlabel('$y$');
    ylabel('$z$');
    zlabel('$x$')
    legend('$\bar{x}$','FOH');
    title({'Comparison of FOH and $\bar{x}$','Position'});
    
    subplot(2,3,2)
    plot3(xbar_(10,:),xbar_(11,:),xbar_(9,:),'--k')
    hold on
    plot3(z_(10,:),z_(11,:),z_(9,:),'-.r')
    xlabel('$q_1$');
    ylabel('$q_2$');
    zlabel('$q_3$');
    title('Quaternion')
    
    subplot(2,3,5)
    plot(tbar_,xbar_(8,:),'--k');
    hold on
    plot(tbar_,z_(8,:),'-.r');
    title('$q_0$');
    xlabel('$t$');
    
    subplot(2,3,3)
    plot3(xbar_(6,:),xbar_(7,:),xbar_(5,:),'--k')
    hold on
    plot3(z_(6,:),z_(7,:),z_(5,:),'-.r')
    xlabel('$v_y$');
    ylabel('$v_z$');
    zlabel('$v_x$');
    title('$v_{\mathcal{I}}$')     
    
    subplot(2,3,6)
    plot3(xbar_(12,:),xbar_(13,:),xbar_(14,:),'--k')
    hold on
    plot3(z_(12,:),z_(13,:),z_(14,:),'-.r')
    xlabel('$\omega_y$');
    ylabel('$\omega_z$');
    zlabel('$\omega_x$');
    title('$\omega_{\mathcal{B}}$') 
    
end