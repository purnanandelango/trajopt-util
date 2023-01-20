function t = time_of_maneuver_via_b(tau,b)
% Compute time of maneuver, given the path velocity b and choice of discretization 
    t = 0;
    for k = 1:length(tau)-1
        t = t + (1/b(k))*diff(tau(k:k+1));
    end
end