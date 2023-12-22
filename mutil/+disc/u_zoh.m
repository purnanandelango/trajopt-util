function u = u_zoh(tt,uvec,tvec)
% Zero-order hold interpolation of a signal defined at discrete nodes
%   u = u_zoh(tt,uvec,tvec)
%   u is the ZOH interpolation at query node(s) tt of signal uvec defined at discrete nodes tvec

    N = length(tvec);
    M = length(tt);
    u = zeros(size(uvec,1),M);

    for k = 1:M 
        t = tt(k);

        j=1;
        while tvec(j)<=t && j<N
            j=j+1;
        end
        j=j-1;
        
        u(:,k) = uvec(:,j); 
    end
    
end