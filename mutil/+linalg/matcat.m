function C = matcat(A,B,dim)
% Concatenate two 2D matrices A and B along first, second or third dimension
[na,ma] = size(A);
[nb,mb] = size(B);
    switch dim
        case 1
            assert(na==nb,"Horizontal concatenation is not possible. Number of rows of A and B should be same.");
            C = [A B];
        case 2
            assert(ma==mb,"Vertical concatenation is not possible. Number of columns of A and B should be same.");
            C = [A; B];
        case 3
            assert(na==nb && ma==mb,"Both dimensions of A and B must match for 3D concatenation.");
            % Columns of A and B are arrange along the third dimension
            C = zeros(na,2,ma);
            for k = 1:ma
                C(:,:,k) = [A(:,k) B(:,k)];
            end
        otherwise
            error("Only upto third dimension concatenation is possible.");
    end
end