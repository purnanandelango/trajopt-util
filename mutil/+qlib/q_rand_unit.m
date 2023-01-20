function a = q_rand_unit()
% Compute a random unit quaternion
% SCALAR-LAST

a = randn(4,1);
a = a./norm(a,2);


end