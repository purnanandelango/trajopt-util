function [Ak,Bk,wk,defect_traj,xbarprop] = compute_zoh_noparam_v3_parallel(tbar,xbar,ubar,h,func,func_linz,varargin)
% Compute ZOH discretization of a nonlinear system for N-1 time intervals with intial conditions xbar(:,1:N-1) and control inputs ubar
% Computation across the N-1 intervals is parallelized with parfor
%   No linearization with respect to system parameters; no system parameters are allowed
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
%   varargin{1}   : specify in-built MATLAB ode solver and its options (optional)
%
%   Ak            : nx x nx x N-1 
%   Bk            : nx x nu x N-1 
%   wk            : nx x N-1
%   defect_traj   : 1  x N-1
%   xbarprop      : nx x N

    [nx,N] = size(xbar);
    nu = size(ubar,1);
    
    Ak  = zeros(nx,nx,N-1);
    Bk = zeros(nx,nu,N-1);
    wk  = zeros(nx,N-1);
    % vk  = zeros(nx,N-1);
    defect_traj(N-1) = 0.0; % records the defect between the propagated and reference traj. 
    xbarprop = zeros(nx,N-1);
    
    nx2 = nx^2;
    nxnu = nx*nu;
    Inx = eye(nx);
    ABw_init = [Inx(:);zeros(nxnu,1);zeros(nx,1)];

    if length(h) == 1 % Same integration step for each subinterval
        h = h*ones(1,N-1);
    end
    
    % SETUP
    tspan = cell(1,N-1);
    ufunc = cell(1,N-1);
    for k = 1:N-1
        tspan{k} = [tbar(k);tbar(k+1)];
        ufunc{k} = @(t) ubar(:,k);
    end

    if nargin == 7
        solver_name = varargin{1}{1};            
        solver_options = varargin{1}{2};

        parfor k = 1:N-1
            zk = [xbar(:,k);ABw_init];

            [~,z_tmp] = feval(solver_name,@(t,z) zoh_ode(t,z,ufunc{k}(t),func,func_linz,nx,nu,nx2,nxnu),tspan{k},zk,solver_options);
            z_ = z_tmp'; 
            zkp1 = z_(:,end);               
           
            xbarprop(:,k) = zkp1(1:nx);
            
            Akmat       = reshape(zkp1(nx+1:nx+nx2),[nx,nx]);
            Bk(:,:,k)   = reshape(zkp1(nx+nx2+1:nx+nx2+nxnu),[nx,nu]);
            wk(:,k)     = zkp1(nx+nx2+nxnu+1:2*nx+nx2+nxnu,1);
            % vk(:,k)     = zkp1(1:nx) - Akmat*xbar(:,k) - Bk(:,:,k)*ubar(:,k);
            Ak(:,:,k)   = Akmat;
    
            % norm(wk(:,k)-vk(:,k))
        end            

    elseif nargin == 6

        parfor k = 1:N-1
            zk = [xbar(:,k);ABw_init];

            [~,z_] = disc.rk4_march(@(t,z,u) zoh_ode(t,z,u,func,func_linz,nx,nu,nx2,nxnu),tspan{k},zk,h(k),ufunc{k});
            zkp1 = z_(:,end);
            
            xbarprop(:,k) = zkp1(1:nx);
            
            Akmat       = reshape(zkp1(nx+1:nx+nx2),[nx,nx]);
            Bk(:,:,k)   = reshape(zkp1(nx+nx2+1:nx+nx2+nxnu),[nx,nu]);
            wk(:,k)     = zkp1(nx+nx2+nxnu+1:2*nx+nx2+nxnu,1);
            % vk(:,k)     = zkp1(1:nx) - Akmat*xbar(:,k) - Bk(:,:,k)*ubar(:,k);
            Ak(:,:,k)   = Akmat;
    
            % norm(wk(:,k)-vk(:,k))
        end            

    else 
        error("Incorrect number of arguments passed.")
    end

    xbarprop = [xbar(:,1),xbarprop];
    for k=1:N-1
        defect_traj(k) = norm(xbarprop(:,k+1)-xbar(:,k+1));
    end

end

function f = zoh_ode(t,z,u,func,func_linz,nx,nu,nx2,nxnu)
% z = [x;PhiA(:);PhiB(:);Phiw];
    
    x     =         z(1               : nx);    
    PhiA  = reshape(z(nx+1            : nx+nx2),[nx,nx]);
    PhiB  = reshape(z(nx+nx2+1        : nx+nx2+nxnu),[nx,nu]);
    Phiw  =         z(nx+nx2+nxnu+1   : nx+nx2+nxnu+nx);    
    
    [A,B,w] = func_linz(t,x,u);

    f1 = A*PhiA;
    f2 = A*PhiB + B;
    f3 = A*Phiw + w;    
    
    f = [func(t,x,u);f1(:);f2(:);f3];    

end
