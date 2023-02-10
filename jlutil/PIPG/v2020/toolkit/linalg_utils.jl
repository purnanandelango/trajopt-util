# Collection of useful linear algebra operations
module la_utils
using LinearAlgebra

# Compute the skew matrix associated with cross product of 3D vector
skew_mat(a) = [0.0  -a[3]  a[2]
			   a[3]  0.0  -a[1]
			  -a[2]  a[1]  0.0] 

# Compute the projection matrix mapping to the kernel of a matrix A ∈ ℜ^(mxn)
# proj_ker(A,n) = Array(Diagonal(ones(n)) .- transpose(A)*inv(A*transpose(A))*A)
proj_ker(A,n) = Array(Diagonal(ones(n)) .- pinv(A)*A)

end