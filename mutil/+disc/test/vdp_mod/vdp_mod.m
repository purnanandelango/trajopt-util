function f = vdp_mod(t,x,u,p)
    f = [p(1)*x(2)*u(1);
         p(1)*u(2)*x(2)*(1-x(1)^2) - p(1)*p(2)*x(1)];
end