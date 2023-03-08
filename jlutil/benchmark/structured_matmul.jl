using LinearAlgebra
using BenchmarkTools

n1 = 4
n = 20

Q11 = randn(n1,n1)
q22 = randn(n-n1)
Qsprs = [Q11 zeros(n1,n-n1);
     zeros(n-n1,n1) Diagonal(q22)]
Qfull = Array(Qsprs)
A = randn(n,n)
A11 = A[1:n1,1:n1]
A12 = A[1:n1,n1+1:n]
A21 = A[n1+1:n,1:n1]
A22 = A[n1+1:n,n1+1:n]
C = zeros(n,n)

function f1(X,Y)
    return X*Y
end

function f2!(Z,X11,x22,Y11,Y12,Y21,Y22,n,n1)
    Z[1:n1,1:n1] .= X11*Y11
    Z[1:n1,n1+1:n] .= X11*Y12
    Z[n1+1:n,1:n1] .= x22 .* Y21
    Z[n1+1:n,n1+1:n] .= x22 .* Y22
end

B = f1(Qfull,A);
f2!(C,Q11,q22,A11,A12,A21,A22,n,n1);

@btime f1($Qfull,$A);

@btime f1($Qsprs,$A);

@btime f2!($C,$Q11,$q22,$A11,$A12,$A21,$A22,$n,$n1);
