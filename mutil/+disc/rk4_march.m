function [t_,x_] = rk4_march(f,tspan,x0,h,u,varargin)
% Integrate nonlinear system using fixed-step 4th order Runge-Kutta method
%   [t_,x_] = rk4_march(f,tspan,x0,h,u,varargin)
%   f       : vector-valued function, RHS of ODE: dx(t)/dt = f(t,x(t),u(t),p) (f has same dimension as x0)
%   tspan   : [t_initial, t_final]
%   x0      : initial condition x(t_initial) (n x 1 vector)
%   h       : RK4 step-size
%   u       : vector-valued control input function accepted by f
%   p       : vector-valued parameter accepted by f (varargin{1})
%   x_      : n x N matrix of integrated trajectory (N is determined by h)
%   t_      : 1 x N vector of discrete time grid over which x_ is defined

% assert(mod(diff(tspan),h) == 0,"Time interval length should be an integral multiple of step length h.");

N = round(diff(tspan)/h)+1;

% Initialization
t_ = linspace(tspan(1),tspan(2),N);
x_(size(x0,1),N) = 0; x_(:,1) = x0;

% RK4 march
if nargin == 6
    p = varargin{1};

    for k = 1:N-1
        k1 = h*f(t_(k),x_(:,k),u(t_(k)),p);
        k2 = h*f(t_(k)+0.5*h,x_(:,k)+0.5*k1,u(t_(k)+0.5*h),p);
        k3 = h*f(t_(k)+0.5*h,x_(:,k)+0.5*k2,u(t_(k)+0.5*h),p);
        k4 = h*f(t_(k)+h,x_(:,k)+k3,u(t_(k)+h),p);
        x_(:,k+1) = x_(:,k) + k1/6 + k2/3 + k3/3 + k4/6;
    end

else % RHS of ODE does not require parameters

    for k = 1:N-1
        k1 = h*f(t_(k),x_(:,k),u(t_(k)));
        k2 = h*f(t_(k)+0.5*h,x_(:,k)+0.5*k1,u(t_(k)+0.5*h));
        k3 = h*f(t_(k)+0.5*h,x_(:,k)+0.5*k2,u(t_(k)+0.5*h));
        k4 = h*f(t_(k)+h,x_(:,k)+k3,u(t_(k)+h));
        x_(:,k+1) = x_(:,k) + k1/6 + k2/3 + k3/3 + k4/6;
    end    

end

end