# 03/05/21
# Oscillating masses

const path_to_tools = "../../toolkit/"
macro load_file(str_val,flg) return :( $flg ? include(string(path_to_tools,$str_val)) : include($str_val) ) end

@load_file "proj_funcs.jl" true			  # Library of projection functions	

# Module defining the parameters for example 12 problem
module eX 
using LinearAlgebra, StaticArrays         # Required packages
using JuMP
using ..proj 							  # Library of projection functions

### User-defined algorithm settings

	const kmax_pipg = 5000			   	  # Max PIPG    iterations
	const kmax_powiter = 50				  # Max power   iterations

	const ϵ_abs = 1e-10					  # Small       absolute tolerance
	const ϵ_powiter = 1e-2				  # Power iter. exit tolerance							(relative)

	if isdefined(Main,:load_termination_criteria)
		const ϵ_pd_JuMP = copy(Main.ϵ_pd_JuMP) 
		const ϵ_gap_JuMP = copy(Main.ϵ_gap_JuMP) 
		const ϵ_primal = copy(Main.ϵ_primal)
		const ϵ_dual = copy(Main.ϵ_dual)
		# const ϵ_admm = copy(Main.ϵ_admm)
	else
		const ϵ_pd_JuMP = 1e-9				  # JuMP solver primal and dual feasibility tolerance	
		const ϵ_gap_JuMP = 1e-9 			  # JuMP solver primal and dual objective gap tolerance
		const ϵ_primal = 1e-4				  # PIPG primal exit tolerance							(relative)
		const ϵ_dual = 1e-4		     		  # PIPG dual   exit tolerance							(relative)
		# const ϵ_admm = 1e-6 				  # ADMM    	exit tolerance
	end
	
	const ρ_admm = 2.0					  # ADMM 		step-size

	const solver_JuMP = :mosek 			  # Choice of solver in JuMP (gurobi,mosek,ecos)
	const err_type = :linf  			  # Type of norm for computing error (l2,linf)

### User-defined problem parameters

	const nx = 20						  # No. of states (nx/2 masses)
	const nu = 10 						  # No. of inputs (nu = nx/2; one for each mass)

	const N = 20   					      # No. of discretization points

	const Δ = 0.1						  # Sampling time (s)
	const tvec = [(j-1)*Δ for j in 1:N]   # Time vector (s)

	# System parameters
	const pmax = 2.0					  # Position	
	const vmax = 0.8					  # Velocity		
	const umax = 2.0					  # Forcing	

	# Reference trajectory and boundary conditions
	const y_unscl = [SVector{nx}(cat(ones(Int64(nx/2)),zeros(Int64(nx/2)),dims=1)) for _ in 1:N]
	const x0_unscl = SVector{nx}(1.5 .* y_unscl[1])
	const xf_unscl = SVector{nx}(y_unscl[N])

	# State segments and scaling
		# Segment 1 (position)
		const ix1 = 1:Int64(nx/2)
		const nx1 = length(ix1)
		const scl_x1 = 1.0

		# Segment 2 (velocity)
		const ix2 = Int64((nx/2)+1):nx
		const nx2 = length(ix2)
		const scl_x2 = 1.0

	# Control segments and scaling
		# No segments
		const scl_u = 1.0

	# Scaling matrices
		const scl_x_mat = Array(Diagonal(cat(scl_x1 .* ones(nx1), scl_x2 .* ones(nx2), dims=1)))
		const scl_x_imat = Array(Diagonal(cat((1/scl_x1) .* ones(nx1), (1/scl_x2) .* ones(nx2), dims=1)))
		const scl_u_mat = Array(Diagonal(scl_u .* ones(nu)))
		const scl_u_imat = Array(Diagonal( (1/scl_u) .* ones(nu) ))

	# ZOH discretization: x(t+1) = Ad*x(t) + Bd*u(t) + gd
		# Unscaled
		const Ad_unscl = SMatrix{nx,nx}([zeros(nx1,nx1) Diagonal(ones(nx1));
										 Tridiagonal(ones(nx1-1),-2 .* ones(nx1),ones(nx1-1)) zeros(nx2,nx2)])

		const Bd_unscl = SMatrix{nx,nu}([zeros(nx1,nx1)
										 Diagonal(ones(nx1,nx1))])

		const gd_unscl = SVector{nx}(zeros(nx))

		# Scaled
		const Ad = SMatrix{nx,nx}(scl_x_imat*Ad_unscl*scl_x_mat)
		const Bd = SMatrix{nx,nu}(scl_x_imat*Bd_unscl*scl_u_mat)
		const gd = SVector{nx}(scl_x_imat*gd_unscl)

		const mAdT = SMatrix{nx,nx}(-1 .* Array(transpose(Ad)))
		const mBdT = SMatrix{nu,nx}(-1 .* Array(transpose(Bd)))

		const x0 = SVector{nx}(scl_x_imat*x0_unscl)
		const y = [SVector{nx}(scl_x_imat*y_unscl[t]) for t = 1:N] 
		const uref = nothing

	# LQR weighting matrices
	const Q = SMatrix{nx,nx}(scl_x_mat*Array(Diagonal( cat((1/scl_x1^2) .* ones(nx1),(1/scl_x2^2) .* ones(nx2),dims=1) ))*scl_x_mat)
	const R = SMatrix{nu,nu}(scl_u_mat*Array(Diagonal( (1/scl_u^2) .* ones(nu) ))*scl_u_mat)
	const Qf = SMatrix{nx,nx}(scl_x_mat*Array(Diagonal( cat((1/scl_x1^2) .* ones(nx1),(1/scl_x2^2) .* ones(nx2),dims=1) ))*scl_x_mat)

	# Step-size parameters (min and max eigenvalues of LQR weight matrices)
	const μ = min(min([R[i,i] for i in 1:nu]...),min([Q[i,i] for i in 1:nx]...),min([Qf[i,i] for i in 1:nx]...))
	const λ = max(max([R[i,i] for i in 1:nu]...),max([Q[i,i] for i in 1:nx]...),max([Qf[i,i] for i in 1:nx]...))

	@assert μ>0 && (λ/μ)<10 "Condition number of objective Hessian is too larger: $(λ/μ). Modify Q, R and Qf"

	# Constraints

		# Bound on position (box)
		const l_pmax = -pmax/scl_x1
		const u_pmax = pmax/scl_x1

		# Bound on velocity (box)
		const l_vmax = -vmax/scl_x2
		const u_vmax = vmax/scl_x2

		# Bound on control input (box)
		const l_umax = -umax/scl_u
		const u_umax = umax/scl_u

### Specifications for projection on constraint sets
	# Projection counter
	const proj_count_limit = SVector{2}([2,1])			# Asymptotic projection estimate for x and u	
	const proj_counter = MVector{2}(randn(2))			# [proj_count_x, proj_count_u]

	const κ1_x1 = MVector{nx1}(randn(nx1))
	const κ2_x1 = MVector{nx1}(randn(nx1))
	const κ1_x2 = MVector{nx2}(randn(nx2))
	const κ2_x2 = MVector{nx2}(randn(nx2))
	const x_y_admm = [nothing for _ in 1:N-1]
	const x_u_admm = [nothing for _ in 1:N-1]
	const u_y_admm = [nothing for _ in 1:N-1]
	const u_u_admm = [nothing for _ in 1:N-1]

	function project_x!(x,z,t)

		# Position
		BLAS.blascopy!(nx1,z,1,κ1_x1,1)
		proj.box!(κ2_x1,κ1_x1,l_pmax,u_pmax,nx1)
		for j in ix1
			x[j] = κ2_x1[j]
		end

		# Velocity
		for i in ix2
			κ1_x2[i-nx1] = z[i]
		end
		proj.box!(κ2_x2,κ1_x2,l_vmax,u_vmax,nx2)
		for i in ix2
			x[i] = κ2_x2[i-nx1]
		end
	end

	const project_u!(u,z,t) = proj.box!(u,z,l_umax,u_umax,nu)

	# Diagnostics projection functions
	function project_x_diagnostic!(x,z,t)
		project_x!(x,z,t)
		proj_counter[1] = 2.0
	end
	function project_u_diagnostic!(u,z,t)
		project_u!(u,z,t)
		proj_counter[2] = 1.0
	end

	# Specify constraints in a JuMP problem
	function set_constr_JuMP!(model,x,u)			

		# Bounds on position
		@constraint(model,[t=1:N-1],l_pmax .* ones(nx1) .≤ x[ix1,t] .≤ u_pmax .* ones(nx1))

		# Bounds on velocity
		@constraint(model,[t=1:N-1],l_vmax .* ones(nx2) .≤ x[ix2,t] .≤ u_vmax .* ones(nx2))

		# Bounds on control
		@constraint(model,[t=1:N-1],l_umax .* ones(nu) .≤ u[:,t] .≤ u_umax .* ones(nu))

	end
end	