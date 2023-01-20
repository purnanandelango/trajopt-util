clearvars
close all
clc

rng default

% Plot figure?
flg0 = true;

% Test rk4_march?
flg = true;

% Type of discretization
% discflg = "FOH";
discflg = "FOH w/o param.";
% discflg = "ZOH";
% discflg = "ZOH w/o param.";

% Instance of compute_foh
foh_type = "v3"; % v1, v2, v3

nx = 2;
nu = 2;
np = 2;
tspan = [0,1]; % Dilated time interval

% Baseline solution
x0 = [1;0]; % Initial condition
ubarfunc = @(t) [1./(t+0.1);exp(-t.^2)]; % Control input profile
pbar = [4;1];

% Compute baseline solution of nonlinear ODE
h = 0.01*diff(tspan);
tbar_ = tspan(1):h:tspan(2);
ubar_ = ubarfunc(tbar_);
[~,xbar_] = disc.rk4_march(@vdp_mod,tspan,x0,h,ubarfunc,pbar);
if flg
    [~,xbar2_] = ode45(@(t,x) vdp_mod(t,x,ubarfunc(t),pbar),tbar_,x0,odeset('AbsTol',1e-8,'RelTol',1e-8));
end

h = 0.01*diff(tbar_(1:2)); % step size for integration that computes FOH matrices

N = size(xbar_,2);
randn_pertx = 0.01*randn(nx,1);
randn_pertu = 0.01*randn(nu,1);
randn_pertp = 0.01*randn(np,1);
z_ = zeros(nx,N); z_(:,1) = x0+randn_pertx;
v_ = ubarfunc(tbar_) + randn_pertu;
p_ = pbar+randn_pertp;

fprintf("Computation of "+discflg+" matrices and linear propagation:\n")
if discflg == "FOH"
    tic
    [Ak,Bmk,Bpk,Sk,wk] = feval("disc.compute_foh_"+foh_type,tbar_,xbar_,ubar_,pbar,h,@vdp_mod,@vdp_mod_linearize);
    for k = 1:N-1
        z_(:,k+1) = Ak(:,:,k)*z_(:,k) + Bmk(:,:,k)*v_(:,k) + Bpk(:,:,k)*v_(:,k+1) + Sk(:,:,k)*p_ + wk(:,k);
    end
    toc
elseif discflg == "FOH w/o param."
    vdp_mod_noparam = @(t,x,u) vdp_mod(t,x,u,pbar);
    tic
    [Ak,Bmk,Bpk,wk] = feval("disc.compute_foh_noparam_"+foh_type,tbar_,xbar_,ubar_,h,vdp_mod_noparam,@(t,x,u) vdp_mod_linearize_noparam(t,x,u,pbar));
    for k = 1:N-1
        z_(:,k+1) = Ak(:,:,k)*z_(:,k) + Bmk(:,:,k)*v_(:,k) + Bpk(:,:,k)*v_(:,k+1) + wk(:,k);
    end
    toc    
elseif discflg == "ZOH"
    tic
    [Ak,Bk,Sk,wk] = disc.compute_zoh(tbar_,xbar_,ubar_,pbar,h,@vdp_mod,@vdp_mod_linearize);
    for k = 1:N-1
        z_(:,k+1) = Ak(:,:,k)*z_(:,k) + Bk(:,:,k)*v_(:,k) + Sk(:,:,k)*p_ + wk(:,k);
    end
    toc
elseif discflg == "ZOH w/o param."
    vdp_mod_noparam = @(t,x,u) vdp_mod(t,x,u,pbar);
    tic
    [Ak,Bk,wk] = disc.compute_zoh_noparam(tbar_,xbar_,ubar_,h,vdp_mod_noparam,@(t,x,u) vdp_mod_linearize_noparam(t,x,u,pbar));
    for k = 1:N-1
        z_(:,k+1) = Ak(:,:,k)*z_(:,k) + Bk(:,:,k)*v_(:,k) + wk(:,k);
    end
    toc        
end

fprintf("Error between FOH linear propagation and nonlinear propagation:\n%.3e\n",norm(xbar2_'-z_,inf));

%%

if flg0
    if flg
    % Validation of results for ode45 and rk4_march
        figure
        plot(xbar2_(:,1),xbar2_(:,2),'-k')
        hold on
        plot(xbar_(1,:),xbar_(2,:),'--r')
        axis equal
        title('Comparison of \texttt{ode45} and \texttt{rk4\_march}');
    end

    % Comparison of xbar and foh solution
    figure
    plot(xbar_(1,:),xbar_(2,:),'-k')
    hold on
    plot(z_(1,:),z_(2,:),'--r')
    axis equal
    legend({'$\bar{x}$',discflg});
    title(horzcat('Comparison of ',discflg,' and $\bar{x}$'));    
end