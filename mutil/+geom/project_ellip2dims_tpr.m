function [A_proj] = project_ellip2dims_tpr(A,dims)
% Project an ellipsoid onto specified dimensions (Taylor P. Reynolds method)
% Ellipsoid boundary representation: y = Au where u is a unit vector
% Representation of ellipsoid used in this code { z | z'*(Q\z) <= 1 }
% First step is to obtain Q from A
% Final setp is to obtain A_proj from Q_proj
% Adapted from: https://github.com/tpreynolds/RSS_2020/blob/master/classes/%40cfga/cfga.m

    % Sort the dimensions
    dims = sort(dims);

    n = size(A,1);
    ndims = numel(dims);    

    % Check if the number of requested dimensions is the same as the original dimension of the ellipsoid
    % This function *cannot* handle the case where dims is the permutation of all the original dimensions
    if ndims == n
        A_proj = A;
        return
    end

    Q = inv( inv(A)' * inv(A) );

    if (ndims>n)
        error('Requestion dimension is lower than the ambient dimension of the ellipsoid.')
        % error('Requested dimensions not compatible with ellipsoid.')
    end
    % Build projection matrix
    In = eye(n);
    Im = eye(ndims);
    T  = zeros(n,ndims);
    for k = 1:ndims
        T(:,k) = In(:,dims(k));
    end
    P = In/Q;
    % The LDL decomposition here is more robust than the cholesky
    % since Yalmip can produce some solutions that have a very
    % small negative eigenvalues. Mathematically this procedure is
    % identical, but ldl does not throw an error if P is not pos
    % def. Moreover, the factor L is not necessarily lower
    % triangular since ldl does some pivoting that a cholesky
    % factorization does not need to do. However, the results are
    % the same.
    % L = chol(P,'lower');
    [L_,D_,P_] = ldl(P);
    L = P_ * L_ * sqrtm(D_);
    A = T' * (In/L)';
    [U,S,~] = svd(A,'econ');
    iS = Im/S;
    P_proj  = U * iS * iS * U';
    Q_proj  = Im/P_proj;
    % Qh_proj = sqrtm(Q_proj);

    A_proj = inv(sqrtm(inv(Q_proj)));

end
