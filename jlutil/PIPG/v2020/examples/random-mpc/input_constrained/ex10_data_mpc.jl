# 02/12/21
# Random MPC

const path_to_tools = "../pipg_toolkit/"
macro load_file(str_val,flg) return :( $flg ? include(string(path_to_tools,$str_val)) : include($str_val) ) end

@load_file "proj_funcs.jl" true			  # library of projection functions	

# module defining the parameters for example 10 problem
module eX 
using LinearAlgebra, StaticArrays         # required packages
using JuMP
using ..proj 							  # library of projection functions

### User-defined algorithm settings

	const kmax_pipg = 5000			   	  # max PIPG    iterations
	const kmax_powiter = 50				  # max power   iterations

	const ϵ_abs = 1e-10					  # small       absolute tolerance
	const ϵ_powiter = 1e-2				  # power iter. exit tolerance							(relative)

	const ϵ_pd_JuMP = 1e-3				  # JuMP solver primal and dual feasibility tolerance	
	const ϵ_gap_JuMP = 1e-3 			  # JuMP solver primal and dual objective gap tolerance

	const ϵ_primal = 1e-3				  # PIPG primal exit tolerance							(relative)
	const ϵ_dual = 1e-3		     		  # PIPG dual   exit tolerance							(relative)
	# const ϵ_admm = 1e-5 				  # admm    	exit tolerance
	
	const ρ_admm = 2.0					  # admm 		step-size

	const solver_JuMP = :gurobi 		  # choice of solver in JuMP (gurobi,mosek,ecos)
	const err_type = :linf  			  # type of norm for computing error (l2,linf)

### User-defined problem parameters

	const nx = 10						  # no. of states
	const nu = 20						  # no. of inputs

	const N = 10   					      # no. of discretization points (prediction horizon)
	const NN = 30						  # no. of dicretization points for main reference

	const ϵ_mpc = 0.001					  # model uncertainty
	const mpc_constr_update = false		  # flag for updating constraints at each MPC solve step

	const Δ = 0.1						  # sampling time (s)
	const tvec = [(j-1)*Δ for j in 1:N]   # time vector (s)

	# system parameters
	const umax = 10.0

	# reference trajectory and boundary conditions
	const yy_unscl = [SVector{nx}(randn(nx)) for _ in 1:NN]
	const x0_unscl_init = SVector{nx}(yy_unscl[1])
	const xf_unscl = SVector{nx}(yy_unscl[N])

	# state segments and scaling
		# no segments
		scl_x = 1.0

	# control segments and scaling
		# no segments
		scl_u = 1.0

	# scaling matrices
		const scl_x_mat = Array(Diagonal(scl_x .* ones(nx)))
		const scl_x_imat = Array(Diagonal( (1/scl_x) .* ones(nx) ))
		const scl_u_mat = Array(Diagonal(scl_u .* ones(nu)))
		const scl_u_imat = Array(Diagonal( (1/scl_u) .* ones(nu) ))

	# ZOH discretization: x(t+1) = Ad*x(t) + Bd*u(t) + gd
		# unscaled
		Ad_tmp = randn(nx,nx); svd_Ad_tmp = svd(Ad_tmp) 
		const Ad_unscl = SMatrix{nx,nx}(Ad_tmp ./ svd_Ad_tmp.S[1])		# marginally stable s

		const Bd_unscl = SMatrix{nx,nu}(randn(nx,nu))

		const gd_unscl = SVector{nx}(zeros(nx))

		# scaled
		const Ad = SMatrix{nx,nx}(scl_x_imat*Ad_unscl*scl_x_mat)
		const Bd = SMatrix{nx,nu}(scl_x_imat*Bd_unscl*scl_u_mat)
		const gd = SVector{nx}(scl_x_imat*gd_unscl)

		const mAdT = SMatrix{nx,nx}(-1 .* Array(transpose(Ad)))
		const mBdT = SMatrix{nu,nx}(-1 .* Array(transpose(Bd)))

		const x0 = MVector{nx}(scl_x_imat*x0_unscl_init)			# scaled initial state updated after plant propagation at each instant
		# const xf = SVector{nx}(scl_x_imat*xf_unscl)
		const y = [MVector{nx}(scl_x_imat*yy_unscl[t]) for t = 1:N] # scaled reference for prediction horizon
		# mpc_utils has functions for updating x0 and y
		uref = nothing

	# LQR weighting matrices
	const Q = SMatrix{nx,nx}(scl_x_mat*Array(Diagonal( (1/scl_x^2) .* ones(nx) ))*scl_x_mat)
	const R = SMatrix{nu,nu}(scl_u_mat*Array(Diagonal( (1/scl_u^2) .* ones(nu) ))*scl_u_mat)
	const Qf = SMatrix{nx,nx}(scl_x_mat*Array(Diagonal( (1/scl_x^2) .* ones(nx) ))*scl_x_mat)

	# step-size parameters (min and max eigenvalues of LQR weight matrices)
	const μ = min(min([R[i,i] for i in 1:nu]...),min([Q[i,i] for i in 1:nx]...),min([Qf[i,i] for i in 1:nx]...))
	const λ = max(max([R[i,i] for i in 1:nu]...),max([Q[i,i] for i in 1:nx]...),max([Qf[i,i] for i in 1:nx]...))

	@assert μ>0 && (λ/μ)<10 "Condition number of objective Hessian is too larger: $(λ/μ). Modify Q, R and Qf"

	# constraints

		# bound on control input (box)
		const l_umax = -umax/scl_u
		const u_umax = umax/scl_u

### Specifications for projection on constraint sets
	# projection counter
	const proj_count_limit = SVector{2}([0,1])			# asymptotic projection estimate for x and u	
	const proj_counter = MVector{2}(randn(2))			# [proj_count_x, proj_count_u]

	const x_y_admm = [nothing for _ in 1:N-1]
	const x_u_admm = [nothing for _ in 1:N-1]
	const u_y_admm = [nothing for _ in 1:N-1]
	const u_u_admm = [nothing for _ in 1:N-1]

	function project_x!(x,z,t)
		# no projection operation
		BLAS.blascopy!(nx,z,1,x,1)
	end

	const project_u!(u,z,t) = proj.box!(u,z,l_umax,u_umax,nu)

	# diagnostics projection functions
	function project_x_diagnostic!(x,z,t)
		project_x!(x,z,t)
		proj_counter[1] = 0.0
	end
	function project_u_diagnostic!(u,z,t)
		project_u!(u,z,t)
		proj_counter[2] = 1.0
	end

	# specify constraints in a JuMP problem
	function set_constr_JuMP!(model,x,u)			

		# bounds on control
		@constraint(model,[t=1:N-1],l_umax .* ones(nu) .≤ u[:,t] .≤ u_umax .* ones(nu))

	end
end	