function u = u_foh(tt,uvec,tvec)
% First-order hold interpolation of a signal defined at discrete nodes
%   u = u_foh(tt,uvec,tvec)
%   u is the FOH interpolation at query node(s) tt of signal uvec defined at discrete nodes 

    N = length(tvec);
    M = length(tt);
    nu = size(uvec,1);
    u = zeros(nu,M);

    for k = 1:M 
        t = tt(k);

        j=1;
        while tvec(j)<=t && j<N
            j=j+1;
        end
        j=j-1;
        
        if abs(tvec(j)-t)>1e-8
            t1 = tvec(j);
            t2 = tvec(j+1);
            u1 = uvec(:,j);
            u2 = uvec(:,j+1);
            u(:,k) = ( (t-t1)*u2 + (t2-t)*u1 )/(t2-t1); % Linear interpolation
        else
            u(:,k) = uvec(:,j);
        end
    
    end
    
end