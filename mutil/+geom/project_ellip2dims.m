function [A_proj] = project_ellip2dims(A,dims)
% Project an ellipsoid onto specified dimensions
% Ellipsoid boundary representation: y = Au where u is a unit vector
% Representation of ellipsoid used in this code { z | z'*Q*z <= 1 }
% First step is to obtain Q from A
% Final setp is to obtain A_proj from Q_proj
% Inspired from https://math.stackexchange.com/a/2438682

    % Sort the dimensions in ascending order
    dims = sort(dims);

    n = size(A,1);
    ndims = numel(dims);

    % Permute columns of A so that the desired dimensions are to the left
    % In = eye(n);
    A = A([dims,setdiff(1:n,dims)],[dims,setdiff(1:n,dims)]);

    if ndims == n
        A_proj = A;
        return
    elseif ndims > n
        error('Requested dimension is higher than the ambient dimension of the ellipsoid.')
    end        

    Q = inv(A)' * inv(A);

    J = Q(1:ndims,1:ndims);
    L = Q(ndims+1:n,1:ndims);
    K = Q(ndims+1:n,ndims+1:n);

    Q_proj = J - L'*(K\L);

    A_proj = eye(ndims)/sqrtm(Q_proj);

end
