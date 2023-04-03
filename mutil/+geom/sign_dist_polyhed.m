function [y,dist_val] = sign_dist_polyhed(x,H,h,varargin)
% Compute projection of x on polyhedron defined by { z | H*z - h <= 0 }
% Solve a QP for computing the projection

    verbose = false;
    if nargin > 3
        verbose = varargin{1};
    end

    if max(H*x-h) > 0 % x is outside the polyhedron

        % Solve QP to compute projection

        % Use OSQP
        n = size(x,1);
        P = 2*eye(n);
        q = -2*x;
        prob = osqp;
        prob.setup(P,q,H,[],h,'alpha',0.1,'verbose',verbose);
        res = prob.solve();
        y = res.x;
        dist_val = norm(y-x);

        % Use YALMIP-ECOS
        % yalmip clear
        % 
        % y = sdpvar(length(x),1);        % Projection point
        % slk = sdpvar(1,1);              % Slack variable
        % 
        % constr = [H*y       <= g;       % Linear inequality constraint
        %           norm(x-y) <= slk];    % Second order cone constraint
        % obj_fun = slk;                  % Objective function
        % 
        % yalmip_out = optimize(constr,obj_fun,sdpsettings('solver','ECOS','verbose',verbose));
        % 
        % if yalmip_out.problem ~= 0
        %     error('Failed to solve convex program for projecting onto polyhedron.');
        % end
        % 
        % y = value(y);
        % dist_val = norm(x-y);

    else % x is inside the polyhedron

        [m,n] = size(H);
        dist_val_halfspace(m,1) = 0;
        y_halfspace(n,m) = 0;

        % Compute projection onto halfspaces which define the complement of the polyhedron
        for j = 1:m
           [y_halfspace(:,j),dist_val_halfspace(j)] = geom.project_halfspace(x,-H(j,:)',-h(j)); 
        end

        % Pick the minimum distance projection point
        [dist_val,idx] = min(dist_val_halfspace);
        y = y_halfspace(:,idx);
        dist_val = -dist_val;

    end

end
