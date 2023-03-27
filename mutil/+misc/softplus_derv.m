function y = softplus_derv(x,gam)
% Derivative of single-variable softplus function with sharpness parameter gamma
    y = exp(gam*x)/(exp(gam*x)+1);
end