function t = time_grid_via_b(tau,b)
% Compute time grid, given the path velocity b and choice of discretization 
    t = zeros(size(tau));
    for k = 1:length(tau)-1
        t(k+1) = t(k) + (1/b(k))*diff(tau(k:k+1));
    end
end