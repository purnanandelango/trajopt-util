function [Ak,Bk,Sk,wk,defect_traj] = compute_zoh(tbar,xbar,ubar,pbar,h,func,func_linz,varargin)
% Compute ZOH discretization of a nonlinear system for N-1 time intervals with intial conditions xbar(:,1:N-1) and control inputs ubar
%   Linearization wrt system parameters pbar (including the time dilation factor if applicable) is considered
%   This function doesn't discriminate between the time dilation factor and other system parameters
%
%   tbar          : 1 x N
%   xbar          : nx x N
%   ubar          : nu x N-1 
%   pbar          : np x N
%   h             : step-size for integration that computes FOH matrices of sub-intervals
%   func          : rhs of system ODE 
%                   func(x,u,p)
%   func_linz     : linearization of rhs of system ODE 
%                   [A,B,S,w] = func_linz(x,u,p)
%   varargin{1}   : structure containing system parameters (optional)
%
%   Ak            : nx x nx x N-1 
%   Bk            : nx x nu x N-1 
%   Sk            : nx x np x N-1 
%   wk            : nx x N-1


    [nx,N] = size(xbar);
    nu = size(ubar,1);
    np = length(pbar);
    
    Ak  = zeros(nx,nx,N-1);
    Bk  = zeros(nx,nu,N-1);
    Sk  = zeros(nx,np,N-1);
    wk  = zeros(nx,N-1);
    vk  = zeros(nx,N-1);
    defect_traj(N-1) = 0.0; % records the defect between the propagated and reference traj. 
    
    nx2 = nx^2;
    nxnu = nx*nu;
    nxnp = nx*np;
    Inx = eye(nx);
    ABSw_init = [Inx(:);zeros(nxnu,1);zeros(nxnp,1);zeros(nx,1)];

    for k = 1:N-1
        zk = [xbar(:,k);ABSw_init];
        tspan = [tbar(k);tbar(k+1)];
        ufunc = @(t) ubar(:,k);

        % Ensure that the integration step is not too small
        h_step = max((1/40)*diff(tspan),h);        
        
        if nargin == 8
           [~,z_] = disc.rk4_march(@(t,z,u,p) zoh_ode(t,z,u,p,func,func_linz,nx,nx2,varargin{1}),tspan,zk,h_step,ufunc,pbar);
        elseif nargin == 7
            [~,z_] = disc.rk4_march(@(t,z,u,p) zoh_ode(t,z,u,p,func,func_linz,nx,nx2),tspan,zk,h_step,ufunc,pbar);
        else
            error("Incorrect no. of arguments for compute_foh.")
        end
        zkp1 = z_(:,end);
        
        defect_traj(k) = norm(zkp1(1:nx) - xbar(:,k+1));
        
        Akmat = reshape(zkp1(nx+1:nx+nx2),[nx,nx]);
        Bk(:,:,k) = Akmat*reshape(zkp1(nx+nx2+1:nx+nx2+nxnu),[nx,nu]);
        Sk(:,:,k) = Akmat*reshape(zkp1(nx+nx2+nxnu+1:nx+nx2+nxnu+nxnp,1),[nx,np]); 
        wk(:,k) = Akmat*zkp1(nx+nx2+nxnu+nxnp+1:2*nx+nx2+nxnu+nxnp,1);
        vk(:,k) = zkp1(1:nx) - Akmat*xbar(:,k) - Bk(:,:,k)*ubar(:,k) - Sk(:,:,k)*pbar;
        Ak(:,:,k) = Akmat;

        % norm(wk(:,k)-vk(:,k))
    end

end

function f = zoh_ode(t,z,u,p,func,func_linz,nx,nx2,varargin)
% z = [x;PhiA(:);Btil(:);Stil(:);wtil];
% p = param
% varargin{1} : structure containing system parameters (optional)

    flg = false;
    if length(varargin)==1
        sys_struct = varargin{1};
        flg = true;
    end
    
    x = z(1:nx);
    param = p;
    PhiA = reshape(z(nx+1:nx+nx2),[nx,nx]);
    PhiAinv = PhiA\eye(nx);
    
    if flg
        [A,B,S,w] = func_linz(t,x,u,param,sys_struct);
    else
        [A,B,S,w] = func_linz(t,x,u,param);
    end
    
    f1 = A*PhiA;
    f2 = PhiAinv*B;
    f3 = PhiAinv*S;
    f4 = PhiAinv*w;
    
    if flg
        f = [func(t,x,u,param,sys_struct);f1(:);f2(:);f3(:);f4];
    else
        f = [func(t,x,u,param);f1(:);f2(:);f3(:);f4];    
    end

end