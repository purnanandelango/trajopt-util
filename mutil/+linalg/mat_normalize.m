function A = mat_normalize(A,flag,varargin)
% Normalize rows or columns of a matrix with a desired norm

    if nargin == 3
        norm_type = varargin{1};
    else
        norm_type = 'inf';
    end

    switch flag

        case 'row'

            n = size(A,1);
            for j = 1:n
                row_norm = norm(A(j,:),norm_type);
                A(j,:) = A(j,:)/row_norm;
            end            

        case 'column'

            m = size(A,2);
            for j = 1:m
                column_norm = norm(A(:,j),norm_type);
                A(:,j) = A(:,j)/column_norm;
            end                        

    end

end