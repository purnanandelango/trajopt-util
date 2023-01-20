function x = ends2interp(x1,x2,tau,flg,m)
% Construct an N-point interpolation between x1 and x2
% tau must be a grid in [0,1] with tau(1) == 0 and tau(N) == 1
    n = length(x1);
    N = length(tau);
    x = zeros(n,N);
    switch flg
        case 'poly'
            for k = 1:N
                x(:,k) = x1*(1-tau(k)^m) + x2*(tau(k)^m);
            end
        case 'exp'
            exp1 = exp(tau(1)*m);
            expN = exp(tau(N)*m); 
            dexp = expN - exp1;
            for k = 1:N
                x(:,k) = x1*(expN - exp(tau(k)*m))/dexp + x2*(exp(tau(k)*m) - exp1)/dexp;
            end
        case 'log'
            pert = 1e-5;
            tau = tau + pert;
            log1 = log(tau(1)*m);
            logN = log(tau(N)*m);
            dlog = logN - log1;
            for k = 1:N
                x(:,k) = x1*(logN - log(tau(k)*m))/dlog + x2*(log(tau(k)*m) - log1)/dlog;
            end
        otherwise
            error("Incorrect flag for interpolation.");
    end
end