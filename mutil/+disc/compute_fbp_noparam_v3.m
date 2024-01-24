function [Ak,Bk,wk,defect_traj,xbarprop] = compute_fbp_noparam_v3(tbar,xbar,ubar,t_burn,func,func_linz,ode_solver)
% Compute finite-burn pulse discretization of a nonlinear system for N-1 time intervals with intial conditions xbar(:,1:N-1) and control inputs ubar
%   No linearization with respect to system parameters
%   Note that this approach completely avoids matrix inversions by reformulating the initial value problem for computing the discrete-time B matrices
%
%   tbar          : 1 x N
%   xbar          : nx x N
%   ubar          : nu x N 
%   t_burn        : control input pulse duration
%   func          : rhs of system ODE 
%                   func(x,u)
%   func_linz     : linearization of rhs of system ODE 
%                   [A,B,w] = func_linz(x,u)
%   ode_solver    : specify in-built MATLAB ode solver and its options
%
%   Ak            : nx x nx x N-1 
%   Bk            : nx x nu x N-1 
%   wk            : nx x N-1
%   defect_traj   : 1  x N-1
%   xbarprop      : nx x N

    [nx,N] = size(xbar);
    nu = size(ubar,1);
    
    Ak  = zeros(nx,nx,N-1);
    Bk  = zeros(nx,nu,N-1);
    wk  = zeros(nx,N-1);
    % vk  = zeros(nx,N-1);
    defect_traj(N-1) = 0.0; % records the defect between the propagated and reference traj. 
    xbarprop = zeros(nx,N); xbarprop(:,1) = xbar(:,1);
    
    nx2 = nx^2;
    nxnu = nx*nu;
    Inx = eye(nx);
    ABw_init = [Inx(:);zeros(nxnu,1);zeros(nx,1)];

    for k = 1:N-1

        zk = [xbar(:,k);ABw_init];
        tspan = tbar(k) + [0, t_burn];
        
        [~,z_] = ode45(@(t,z) fbp_ode(t,z,ubar(:,k),1,func,func_linz,nx,nu,nx2,nxnu),tspan,zk,odeset('AbsTol',1e-5,'RelTol',1e-4));

        tspan = [tbar(k) + t_burn, tbar(k+1)]; 
        
        [~,z_] = feval(ode_solver{1},@(t,z) fbp_ode(t,z,zeros(nu,1),0,func,func_linz,nx,nu,nx2,nxnu),tspan,z_(end,:)',ode_solver{2});

        zkp1 = z_(end,:)';
        
        defect_traj(k)  = norm(zkp1(1:nx) - xbar(:,k+1));
        xbarprop(:,k+1) = zkp1(1:nx);
        
        Akmat       = reshape(zkp1(nx+1         :nx+nx2),        [nx,nx]);
        Bk(:,:,k)   = reshape(zkp1(nx+nx2+1     :nx+nx2+nxnu),   [nx,nu]);
        wk(:,k)     =         zkp1(nx+nx2+nxnu+1:nx+nx2+nxnu+nx);
        % vk(:,k)     = zkp1(1:nx) - Akmat*xbar(:,k) - Bmk(:,:,k)*ubar(:,k) - Bpk(:,:,k)*ubar(:,k+1);
        Ak(:,:,k)   = Akmat;

        % norm(wk(:,k)-vk(:,k))
    end

end

function f = fbp_ode(t,z,u,burn_active,func,func_linz,nx,nu,nx2,nxnu)
% z = [x;PhiA(:);PhiB(:);Phiw];
% burn_active is either 0 or 1

    x     =         z(1               : nx);    
    PhiA  = reshape(z(nx+1            : nx+nx2),[nx,nx]);
    PhiB  = reshape(z(nx+nx2+1        : nx+nx2+nxnu),[nx,nu]);
    Phiw  =         z(nx+nx2+nxnu+1   : nx+nx2+nxnu+nx);    

    [A,B,w] = func_linz(t,x,u);

    f1 = A*PhiA;
    f2 = A*PhiB + B*burn_active;
    f3 = A*Phiw + w;    

    f = [func(t,x,u);f1(:);f2(:);f3];

end