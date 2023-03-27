function y = softplus(x,gam)
% Single-variable softplus function with sharpness parameter gamma
    y = log(exp(gam*x)+1)/gam;
end