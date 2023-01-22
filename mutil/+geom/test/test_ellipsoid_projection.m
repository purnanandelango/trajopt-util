clearvars
clc
close all

n = 6;
A = randn(n);
P = A*A';

% dims = randi([1,6],[1,3]);
dims = [3,5,6]
% dims = 4:6;
% dims = 1:3;

A_proj1 = project_ellipsoid2dims(A,dims)
A_proj2 = project_ellipsoid2dims_tpr(A,dims)

norm(A_proj2-A_proj1)
