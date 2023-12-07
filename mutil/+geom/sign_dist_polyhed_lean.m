function [y,dist_val,dist_jac] = sign_dist_polyhed_lean(x,H,h)
% Compute projection of x on polyhedron defined by { z | H*z - h <= 0 }
% Solve a QP for computing the projection

    [m,n] = size(H);

    if max(H*x-h) > 0 % x is outside the polyhedron

        % Solve QP to compute projection

        % Use quadprog
        P = 2*eye(n);
        q = -2*x;
        opts = optimoptions('quadprog','Algorithm','active-set',...
               'ConstraintTolerance',1e-7,'OptimalityTolerance',1e-7,'Display','off');
        y = quadprog(P,q,H,h,[],[],[],[],ones(n,1),opts);
        dist_val = norm(y-x);

    else % x is inside the polyhedron

        dist_val_halfspace = zeros(m,1);
        y_halfspace = zeros(n,m);

        % Compute projection onto halfspaces which define the complement of the polyhedron
        for j = 1:m
           [y_halfspace(:,j),dist_val_halfspace(j)] = geom.project_halfspace(x,-H(j,:)',-h(j)); 
        end

        % Pick the minimum distance projection point
        [dist_val,idx] = min(dist_val_halfspace);
        y = y_halfspace(:,idx);
        dist_val = -dist_val;

    end

    % Computation of signed distance Jacobian wrt x
    % Jacobian always has unit norm
    if abs(dist_val) >= 1e-5
        dist_jac = (x-y)/dist_val;
    else % Find a hyperplane of the polyhedron which contains x
        err_val = zeros(1,m);
        for j = 1:m
            err_val(j) = abs(H(j,:)*x-h(j));
        end
        [~,idx] = min(err_val);
        dist_jac = H(idx(1),:)'/norm(H(idx(1),:));
    end

end