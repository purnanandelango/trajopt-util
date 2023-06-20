function [Ak,Bmk,Bpk,wk,defect_traj] = compute_foh_noparam_v3(tbar,xbar,ubar,h,func,func_linz,varargin)
% Compute FOH discretization of a nonlinear system for N-1 time intervals with intial conditions xbar(:,1:N-1) and control inputs ubar
%   No linearization with respect to system parameters
%   Note that this approach completely avoids matrix inversions by reformulating the initial value problem for computing the discrete-time B matrices
%
%   tbar          : 1 x N
%   xbar          : nx x N
%   ubar          : nu x N 
%   h             : step-size for integration that computes FOH matrices of sub-intervals
%   func          : rhs of system ODE 
%                   func(x,u)
%   func_linz     : linearization of rhs of system ODE 
%                   [A,B,w] = func_linz(x,u)
%   varargin{1}   : structure containing system parameters (optional)
%
%   Ak            : nx x nx x N-1 
%   Bmk           : nx x nu x N-1 
%   Bpk           : nx x nu x N-1 
%   wk            : nx x N-1


    [nx,N] = size(xbar);
    nu = size(ubar,1);
    
    Ak  = zeros(nx,nx,N-1);
    Bmk = zeros(nx,nu,N-1);
    Bpk = zeros(nx,nu,N-1);
    wk  = zeros(nx,N-1);
    vk  = zeros(nx,N-1);
    defect_traj(N-1) = 0.0; % records the defect between the propagated and reference traj. 
    
    nx2 = nx^2;
    nxnu = nx*nu;
    Inx = eye(nx);
    ABmBpw_init = [Inx(:);zeros(nxnu,1);zeros(nxnu,1);zeros(nx,1)];

    for k = 1:N-1
        zk = [xbar(:,k);ABmBpw_init];
        tspan = [tbar(k);tbar(k+1)];
        ufunc = @(t) ( ubar(:,k)*(tbar(k+1)-t) + ubar(:,k+1)*(t-tbar(k)) )/( tbar(k+1) - tbar(k) );

        % Ensure that the integration step is not too small
        % h_step = max((1/40)*diff(tspan),h);
        h_step = h;
        
        if nargin == 7
           [~,z_] = disc.rk4_march(@(t,z,u,p) foh_ode(t,z,u,p,func,func_linz,nx,nu,nx2,nxnu,varargin{1}),tspan,zk,h_step,ufunc,tspan);
        elseif nargin == 6
            [~,z_] = disc.rk4_march(@(t,z,u,p) foh_ode(t,z,u,p,func,func_linz,nx,nu,nx2,nxnu),tspan,zk,h_step,ufunc,tspan);
        else
            error("Incorrect no. of arguments for compute_foh.")
        end
        zkp1 = z_(:,end);
        
        defect_traj(k) = norm(zkp1(1:nx) - xbar(:,k+1));
        
        Akmat       = reshape(zkp1(nx+1:nx+nx2),[nx,nx]);
        Bmk(:,:,k)  = reshape(zkp1(nx+nx2+1:nx+nx2+nxnu),[nx,nu]);
        Bpk(:,:,k)  = reshape(zkp1(nx+nx2+nxnu+1:nx+nx2+2*nxnu),[nx,nu]);
        wk(:,k)     = zkp1(nx+nx2+2*nxnu+1:2*nx+nx2+2*nxnu,1);
        vk(:,k)     = zkp1(1:nx) - Akmat*xbar(:,k) - Bmk(:,:,k)*ubar(:,k) - Bpk(:,:,k)*ubar(:,k+1);
        Ak(:,:,k)   = Akmat;

        % norm(wk(:,k)-vk(:,k))
    end

end

function f = foh_ode(t,z,u,p,func,func_linz,nx,nu,nx2,nxnu,varargin)
% z = [x;PhiA(:);PhiBm(:);PhiBp(:);Phiw];
% p = [tk;tkp1]
% varargin{1} : structure containing system parameters (optional)

    flg = false;
    if length(varargin)==1
        sys_struct = varargin{1};
        flg = true;
    end
    
    x     =         z(1               : nx);    
    PhiA  = reshape(z(nx+1            : nx+nx2),[nx,nx]);
    PhiBm = reshape(z(nx+nx2+1        : nx+nx2+nxnu),[nx,nu]);
    PhiBp = reshape(z(nx+nx2+nxnu+1   : nx+nx2+2*nxnu),[nx,nu]);
    Phiw  =         z(nx+nx2+2*nxnu+1 : nx+nx2+2*nxnu+nx);    
    
    if flg
        [A,B,w] = func_linz(t,x,u,sys_struct);
    else
        [A,B,w] = func_linz(t,x,u);
    end
    
    tkp1 = p(end);
    tk = p(end-1);
    delta_t = tkp1-tk;

    f1 = A*PhiA;
    f2 = A*PhiBm + B*(tkp1-t)/delta_t;
    f3 = A*PhiBp + B*(t-tk)/delta_t;
    f4 = A*Phiw  + w;    
    
    if flg
        f = [func(t,x,u,sys_struct);f1(:);f2(:);f3(:);f4];
    else
        f = [func(t,x,u);f1(:);f2(:);f3(:);f4];    
    end

end