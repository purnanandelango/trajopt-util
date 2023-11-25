function cnstr = time_cnstr(s,dtau,bounds,disc)
% Create constraints on time interval length and time of maneuver
% s is a YALMIP sdpvar for unscaled dilation factor

    dtmin = bounds{1};
    dtmax = bounds{2};

    Km1 = length(dtau); % K-1
    cnstr = [];

    ToF = 0;
    switch disc
        case "ZOH"
            for k = 1:Km1
                ToF = ToF + dtau(k)*s(k);
                cnstr = [cnstr; dtmin <= dtau(k)*s(k) <= dtmax]; 
            end
        case "FOH"
            for k = 1:Km1
                ToF = ToF + 0.5*dtau(k)*(s(k+1)+s(k));
                cnstr = [cnstr; dtmin <= 0.5*dtau(k)*(s(k+1)+s(k)) <= dtmax];
            end
    end    

    % Time of maneuver upper bound
    if length(bounds) == 3
        ToFmax = bounds{3};
        cnstr = [cnstr;ToF <= ToFmax];        
    end
    
end