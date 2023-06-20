function [Ak,Bk,wk,defect_traj] = compute_zoh_noparam(tbar,xbar,ubar,h,func,func_linz,varargin)
% Compute ZOH discretization of a nonlinear system for N-1 time intervals with intial conditions xbar(:,1:N-1) and control inputs ubar
%   No linearization with respect to system parameters
%
%   tbar          : 1 x N
%   xbar          : nx x N
%   ubar          : nu x N-1 
%   h             : step-size for integration that computes FOH matrices of sub-intervals
%   func          : rhs of system ODE 
%                   func(x,u)
%   func_linz     : linearization of rhs of system ODE 
%                   [A,B,w] = func_linz(x,u)
%   varargin{1}   : structure containing system parameters (optional)
%
%   Ak            : nx x nx x N-1 
%   Bk            : nx x nu x N-1 
%   wk            : nx x N-1


    [nx,N] = size(xbar);
    nu = size(ubar,1);
    
    Ak  = zeros(nx,nx,N-1);
    Bk  = zeros(nx,nu,N-1);
    wk  = zeros(nx,N-1);
    vk  = zeros(nx,N-1);
    defect_traj(N-1) = 0.0; % records the defect between the propagated and reference traj. 
    
    nx2 = nx^2;
    nxnu = nx*nu;
    Inx = eye(nx);
    ABw_init = [Inx(:);zeros(nxnu,1);zeros(nx,1)];

    for k = 1:N-1
        zk = [xbar(:,k);ABw_init];
        tspan = [tbar(k);tbar(k+1)];
        ufunc = @(t) ubar(:,k);

        % Ensure that the integration step is not too small
        h_step = max((1/40)*diff(tspan),h);        
        
        if nargin == 7
           [~,z_] = disc.rk4_march(@(t,z,u) zoh_ode(t,z,u,func,func_linz,nx,nx2,varargin{1}),tspan,zk,h_step,ufunc);
        elseif nargin == 6
            [~,z_] = disc.rk4_march(@(t,z,u) zoh_ode(t,z,u,func,func_linz,nx,nx2),tspan,zk,h_step,ufunc);
        else
            error("Incorrect no. of arguments for compute_foh.")
        end
        zkp1 = z_(:,end);
        
        defect_traj(k) = norm(zkp1(1:nx) - xbar(:,k+1));
        
        Akmat = reshape(zkp1(nx+1:nx+nx2),[nx,nx]);
        Bk(:,:,k) = Akmat*reshape(zkp1(nx+nx2+1:nx+nx2+nxnu),[nx,nu]);
        wk(:,k) = Akmat*zkp1(nx+nx2+nxnu+1:2*nx+nx2+nxnu,1);
        vk(:,k) = zkp1(1:nx) - Akmat*xbar(:,k) - Bk(:,:,k)*ubar(:,k);
        Ak(:,:,k) = Akmat;

        % norm(wk(:,k)-vk(:,k))
    end

end

function f = zoh_ode(t,z,u,func,func_linz,nx,nx2,varargin)
% z = [x;PhiA(:);Btil(:);Stil(:);wtil];
% varargin{1} : structure containing system parameters (optional)

    flg = false;
    if length(varargin)==1
        sys_struct = varargin{1};
        flg = true;
    end
    
    x = z(1:nx);
    PhiA = reshape(z(nx+1:nx+nx2),[nx,nx]);
    PhiAinv = PhiA\eye(nx);
    
    if flg
        [A,B,w] = func_linz(t,x,u,sys_struct);
    else
        [A,B,w] = func_linz(t,x,u);
    end
    
    f1 = A*PhiA;
    f2 = PhiAinv*B;
    f3 = PhiAinv*w;
    
    if flg
        f = [func(t,x,u,sys_struct);f1(:);f2(:);f3];
    else
        f = [func(t,x,u);f1(:);f2(:);f3];    
    end

end