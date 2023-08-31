function b = row_func(a,func)
% Operate on each row of matrix a matrix
    N = size(a,1);
    b(N,1) = 0;
    for j = 1:N
        b(j) = func(a(j,1));
    end
end