function t = time_grid(disc,tau,s)
% Compute time grid, given the dilation factor and choice of discretization 
    N = length(tau);
    t = zeros(1,N);
    switch disc
        case "ZOH"
            for k = 1:N-1
                t(k+1) = t(k) + (tau(k+1)-tau(k))*s(k);
            end
        case "FOH"
            for k = 1:N-1
                t(k+1) = t(k) + 0.5*(s(k+1)+s(k))*(tau(k+1)-tau(k));
            end
    end
end