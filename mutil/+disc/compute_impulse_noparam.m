function [Ak,wk,defect_traj,xbarprop] = compute_impulse_noparam(tbar,xbar,ubar,Eu2x,h,func,func_linz,varargin)
% Compute impulse input discretization of a nonlinear system for N-1 time intervals with intial conditions xbar(:,1:N-1) and control inputs ubar
%   No linearization with respect to system parameters
%
%   tbar          : 1 x N
%   xbar          : nx x N
%   ubar          : nu x N 
%   Eu2x          : input-to-state matrix 
%   h             : step-size for integration that computes FOH matrices of sub-intervals
%   func          : rhs of system ODE 
%                   func(x,u)
%   func_linz     : linearization of rhs of system ODE 
%                   [A,B,w] = func_linz(x,u)
%   varargin{1}   : specify in-built MATLAB ode solver and its options (optional)
%
%   Ak            : nx x nx x N-1 
%   wk            : nx x N-1


    [nx,N] = size(xbar);
    nu = size(ubar,1);

    if length(h) == 1 % Same integration step for each subinterval
        h = h*ones(1,N-1);
    end    
    
    Ak  = zeros(nx,nx,N-1);
    wk  = zeros(nx,N-1);
    vk  = zeros(nx,N-1);
    defect_traj(N-1) = 0.0; % records the defect between the propagated and reference traj. 
    xbarprop = zeros(nx,N); xbarprop(:,1) = xbar(:,1);
    
    zeros_nu = zeros(nu,1);
    nx2      = nx^2;
    Inx      = eye(nx);
    Aw_init  = [Inx(:);zeros(nx,1)];

    for k = 1:N-1
        zk = [xbar(:,k) + Eu2x*ubar(:,k);
              Aw_init];
        tspan = [tbar(k);tbar(k+1)];

        % Ensure that the integration step is not too small
        % h_step = max((1/40)*diff(tspan),h(k)); 
        h_step = h(k);
        
        if nargin == 8
            [~,z_tmp] = feval(varargin{1}{1},@(t,z) impulse_ode(t,z,zeros_nu,func,func_linz,nx,nx2),tspan,zk,varargin{1}{2});
            z_ = z_tmp';
        elseif nargin == 7
            [~,z_] = disc.rk4_march(@(t,z,u) impulse_ode(t,z,zeros_nu,func,func_linz,nx,nx2),tspan,zk,h_step,@(t) zeros_nu);
        else
            error("Incorrect no. of arguments for compute_impulse_noparam.");
        end
        zkp1 = z_(:,end);
        
        defect_traj(k)  = norm(zkp1(1:nx) - xbar(:,k+1));
        xbarprop(:,k+1) = zkp1(1:nx);
        
        Akmat       = reshape(zkp1(nx+1:nx+nx2),[nx,nx]);
        wk(:,k)     = zkp1(nx+nx2+1:2*nx+nx2,1);
        vk(:,k)     = zkp1(1:nx) - Akmat*(xbar(:,k) + Eu2x*ubar(:,k));
        Ak(:,:,k)   = Akmat;

        % norm(wk(:,k)-vk(:,k))
    end

end

function f = impulse_ode(t,z,u,func,func_linz,nx,nx2)
% z = [x;PhiA(:);Phiw];
% u = zeros(nu,1)
    
    x     =         z(1               : nx);    
    PhiA  = reshape(z(nx+1            : nx+nx2),[nx,nx]);
    Phiw  =         z(nx+nx2+1        : nx+nx2+nx);    
    
    [A,~,w] = func_linz(t,x,u);

    f1 = A*PhiA;
    f2 = A*Phiw  + w;    
    
    f = [func(t,x,u);f1(:);f2];

end