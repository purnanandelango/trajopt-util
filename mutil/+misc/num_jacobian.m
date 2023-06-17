function dfun = num_jacobian(fun,z)
% Compute central difference Jacobian of a function
    dfun = [];
    n = length(z);
    In = eye(n);
    epsabs = 1e-8;
    epsrel = 1e-6; 

    pert = max(epsabs,epsrel*z);    
    for j = 1:n        
        fp = fun(z+pert(j)*In(:,j));
        fm = fun(z-pert(j)*In(:,j));
        dfpm = fp-fm;
        if norm(dfpm) < epsabs 
            dfun = [dfun, zeros(size(dfpm))];
        else
            dfun = [dfun,0.5*(fp-fm)/pert(j)];
        end
    end
end
