function x_ = rk4_nocontrol(f,t_,x0)
% Integrate a dynamical system (with no control inputs) using fixed-step 4th order Runge-Kutta method
%   x_      = rk4_march(f,t_,x0)
%   f       : vector-valued function, RHS of ODE: dx(t)/dt = f(t,x(t)) (f has same dimension as x0)
%   t_      : grid of t with potentially nonuniform spacing 
%   x0      : initial condition x(t_(1)) (n x 1 vector)

% Initialization
N = length(t_);
x_(size(x0,1),N) = 0; x_(:,1) = x0;

% RK4 march
for k = 1:N-1
    h  = t_(k+1)-t_(k);
    k1 = h*f(t_(k),x_(:,k));
    k2 = h*f(t_(k)+0.5*h,x_(:,k)+0.5*k1);
    k3 = h*f(t_(k)+0.5*h,x_(:,k)+0.5*k2);
    k4 = h*f(t_(k)+h,x_(:,k)+k3);
    x_(:,k+1) = x_(:,k) + k1/6 + k2/3 + k3/3 + k4/6;
end    

end