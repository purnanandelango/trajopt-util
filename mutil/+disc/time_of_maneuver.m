function ToF = time_of_maneuver(disc,tau,s)
% Compute time of maneuver, given the dilation factor and choice of discretization 
    K = length(tau);
    ToF = 0;
    switch disc
        case "ZOH"
            for k = 1:K-1
                ToF = ToF + (tau(k+1)-tau(k))*s(k);
            end
        case "FOH"
            for k = 1:K-1
                ToF = ToF + 0.5*(s(k+1)+s(k))*(tau(k+1)-tau(k));
            end
    end
end