function [y,dist_val] = sign_dist_ellip_exterior(x,A,varargin)
% Compute projection of an exterior point x onto ellipsoid { z | norm(A\z) <= 1 } 
% Solve a convex program via YALMIP for computing the projection
% This function will return an error when x is inside the ellipsoid
    
% The ellipsoid can also be defined by { z | z'*inv(P)*z <= 1 }, where
%         P = inv( inv(A)'*inv(A) ) = A*A'

    yalmip clear
    y = sdpvar(length(x),1);        % Projection point
    slk = sdpvar(1,1);              % Slack variable

    assert(norm(A\x) > 1,'Query point x cannot be inside the elliposid.')

    % Should diagnostic information be printed?
    if nargin == 2
        verbose = false;
    else
        verbose = varargin{1};
    end

    constr = [norm(A\y) <= 1;       % Ellipsoid constraintv (SOC)
              norm(x-y) <= slk];    % Penalize distance x (SOC)
    obj_fun = slk;                  % Objective function   

    yalmip_out = optimize(constr,obj_fun,sdpsettings('solver','ECOS','verbose',verbose));
   
    if yalmip_out.problem ~= 0
        error('Failed to solve convex program for projecting onto elliopsoid.');
    end

    y = value(y);
    dist_val = norm(x-y);
    
end
