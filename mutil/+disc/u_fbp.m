function u = u_fbp(tt,uvec,tvec,t_burn)
% Finite-burn pulse interpolation of a signal defined at discrete nodes
%   u = u_fbp(tt,uvec,tvec,t_burn)
%   u is FBP interpolation at query node(s) tt of signal uvec defined at discrete nodes tvec

    % N = length(tvec);
    M = length(tt);
    u = zeros(size(uvec,1),M);

    for k = 1:M 
        t = tt(k);

        [~,idx] = mink(abs(t-tvec),2);
        if tvec(idx(1)) < tvec(idx(2))
            tj = tvec(idx(1));
            % tjp1 = tvec(idx(2));
            j = idx(1);
        else
            tj = tvec(idx(2));
            % tjp1 = tvec(idx(1));
            j = idx(2);
        end

        if t <= tj+t_burn
            u(:,k) = uvec(:,j);
        end

    end    
end