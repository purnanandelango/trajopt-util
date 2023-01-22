function [y,dist_val] = project_halfspace(x,a,b)
% Compute projection of x onto the halfspace { z | a^T z <= b }
    
    constr_viol = a'*x - b;

    if constr_viol <= 0
        y = x;
    else
        y = x - constr_viol*a/norm(a)^2;
    end

    dist_val = norm(x-y);
    
end
