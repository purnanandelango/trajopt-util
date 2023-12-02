function A = mat_normalize(A,flag)

    switch flag

        case 'row'

            n = size(A,1);
            for j = 1:n
                row_norm = norm(A(j,:),'inf');
                A(j,:) = A(j,:)/row_norm;
            end            

        case 'column'

            m = size(A,2);
            for j = 1:m
                column_norm = norm(A(:,j),'inf');
                A(:,j) = A(:,j)/column_norm;
            end                        

    end

end