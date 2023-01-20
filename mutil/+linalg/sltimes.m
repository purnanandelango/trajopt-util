function C = sltimes(A,B)
% Slice-wise product of 3-dimensional matrix A with 2-dimensional matrix B
% Right multiply all slices of A (described by third dimension) with B
[n1,n2,n3] = size(A);
[m1,m2] = size(B);

    assert(min(size(A))>1,"A must be a (nontrivial) three dimensional matrix (tensor).");
    assert(n2==m1,"Check dimensions. It is not possible to left-multiply B with slices of A.");

    if m2 == 1
        C = zeros(n1,n3);
        for j = 1:n3
            C(:,j) = A(:,:,j)*B;
        end
    else
        C = zeros(n1,m2,n3);
        for j = 1:n3
            C(:,:,j) = A(:,:,j)*B;
        end
    end

end