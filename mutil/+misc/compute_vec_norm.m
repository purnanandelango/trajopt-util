function nrm_v = compute_vec_norm(v)
% Compute norm of a grid of vectors arranged along column-wise
    N = size(v,2);
    nrm_v = zeros(1,N);
    for k = 1:N
        nrm_v(k) = norm(v(:,k));
    end
end