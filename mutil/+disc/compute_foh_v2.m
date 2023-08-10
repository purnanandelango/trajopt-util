function [Ak,Bmk,Bpk,Sk,wk,defect_traj] = compute_foh_v2(tbar,xbar,ubar,pbar,h,func,func_linz,varargin)
% Compute FOH discretization of a nonlinear system for N-1 time intervals with intial conditions xbar(:,1:N-1) and control inputs ubar
%   Linearization wrt system parameters pbar (including the time dilation factor if applicable) is considered
%   This function doesn't discriminate between the time dilation factor and other system parameters
%   Note that the inverse of the state transition matrix is directly integrated to avoid inversion operation inside foh_ode
%
%   tbar          : 1 x N
%   xbar          : nx x N
%   ubar          : nu x N 
%   pbar          : np x N
%   h             : step-size for integration that computes FOH matrices of sub-intervals
%   func          : rhs of system ODE 
%                   func(x,u,p)
%   func_linz     : linearization of rhs of system ODE 
%                   [A,B,S,w] = func_linz(x,u,p)
%   varargin{1}   : specify in-built MATLAB ode solver (optional)
%
%   Ak            : nx x nx x N-1 
%   Bmk           : nx x nu x N-1 
%   Bpk           : nx x nu x N-1 
%   Sk            : nx x np x N-1 
%   wk            : nx x N-1


    [nx,N] = size(xbar);
    nu = size(ubar,1);
    np = length(pbar);

    if length(h) == 1 % Same integration step for each subinterval
        h = h*ones(1,N-1);
    end    
    
    Ak  = zeros(nx,nx,N-1);
    Bmk = zeros(nx,nu,N-1);
    Bpk = zeros(nx,nu,N-1);
    Sk  = zeros(nx,np,N-1);
    wk  = zeros(nx,N-1);
    vk  = zeros(nx,N-1);
    defect_traj(N-1) = 0.0; % records the defect between the propagated and reference traj. 
    
    nx2 = nx^2;
    nxnu = nx*nu;
    nxnp = nx*np;
    Inx = eye(nx);
    ABmBpSw_init = [Inx(:);zeros(nxnu,1);zeros(nxnu,1);zeros(nxnp,1);zeros(nx,1)];

    for k = 1:N-1
        zk = [xbar(:,k);ABmBpSw_init];
        tspan = [tbar(k);tbar(k+1)];
        ufunc = @(t) ( ubar(:,k)*(tbar(k+1)-t) + ubar(:,k+1)*(t-tbar(k)) )/( tbar(k+1) - tbar(k) );

        % Ensure that the integration step is not too small
        % h_step = max((1/40)*diff(tspan),h(k));    
        h_step = h(k);
        
        if nargin == 8
           [~,z_tmp] = feval(varargin{1},@(t,z) foh_ode(t,z,ufunc(t),[pbar;tspan],func,func_linz,nx,nx2),tspan,zk);
           z_ = z_tmp'; 
        elseif nargin == 7
            [~,z_] = disc.rk4_march(@(t,z,u,p) foh_ode(t,z,u,p,func,func_linz,nx,nx2),tspan,zk,h_step,ufunc,[pbar;tspan]);
        else
            error("Incorrect no. of arguments for compute_foh_v2.");
        end
        zkp1 = z_(:,end);
        
        defect_traj(k) = norm(zkp1(1:nx) - xbar(:,k+1));
        
        Akmat = reshape(zkp1(nx+1:nx+nx2),[nx,nx])\Inx;
        Bmk(:,:,k) = Akmat*reshape(zkp1(nx+nx2+1:nx+nx2+nxnu),[nx,nu]);
        Bpk(:,:,k) = Akmat*reshape(zkp1(nx+nx2+nxnu+1:nx+nx2+2*nxnu),[nx,nu]);
        Sk(:,:,k) = Akmat*reshape(zkp1(nx+nx2+2*nxnu+1:nx+nx2+2*nxnu+nxnp,1),[nx,np]); 
        wk(:,k) = Akmat*zkp1(nx+nx2+2*nxnu+nxnp+1:2*nx+nx2+2*nxnu+nxnp,1);
        vk(:,k) = zkp1(1:nx) - Akmat*xbar(:,k) - Bmk(:,:,k)*ubar(:,k) - Bpk(:,:,k)*ubar(:,k+1) - Sk(:,:,k)*pbar;
        Ak(:,:,k) = Akmat;

        % norm(wk(:,k)-vk(:,k))
    end

end

function f = foh_ode(t,z,u,p,func,func_linz,nx,nx2)
% z = [x;invPhiA(:);Btilm(:);Btilp(:);Stil(:);wtil];
% p = [param;tk;tkp1]
    
    x = z(1:nx);
    param = p(1:end-2);
    invPhiA = reshape(z(nx+1:nx+nx2),[nx,nx]);
    
    [A,B,S,w] = func_linz(t,x,u,param);
    
    tkp1 = p(end);
    tk = p(end-1);
    delta_t = tkp1-tk;
    
    f1 = -invPhiA*A;
    f2 = invPhiA*B*(tkp1-t)/delta_t;
    f3 = invPhiA*B*(t-tk)/delta_t;
    f4 = invPhiA*S;
    f5 = invPhiA*w;
    
    f = [func(t,x,u,param);f1(:);f2(:);f3(:);f4(:);f5];    

end