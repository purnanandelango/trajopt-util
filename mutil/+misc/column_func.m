function b = column_func(a,func)
% Operate on each column of matrix a matrix
    N = size(a,2);
    b(1,N) = 0;
    for j = 1:N
        b(j) = func(a(:,j));
    end
end