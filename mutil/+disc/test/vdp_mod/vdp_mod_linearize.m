function [A,B,S,w] = vdp_mod_linearize(t,x,u,p)
% Linearize about (x,u,p)

A = [ 0,                                  p(1)*u(1); 
     -p(1)*p(2)-2*p(1)*u(2)*x(1)*x(2),    p(1)*u(2)*(1-x(1)^2)];

B = [p(1)*x(2),     0; 
      0,            p(1)*x(2)-p(1)*x(2)*x(1)^2];

S = [x(2)*u(1),     0;
     u(2)*x(2)-u(2)*x(2)*x(1)^2-p(2)*x(1),  -p(1)*x(1)];

w = vdp_mod(t,x,u,p) - A*x - B*u - S*p;

end