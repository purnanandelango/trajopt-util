# 01/29/21
# Grasp force optimization for a 3D block with a two-finger manipulator

include("../../../pipg_toolkit/proj_funcs.jl") 	  # library of projection functions
include("../../../pipg_toolkit/linalg_utils.jl") 	  # library of linear algebra utilities 	

# module defining the parameters for example 9a problem
module eX 
using LinearAlgebra, StaticArrays         # required packages
using JuMP
using ..proj 							  # library of projection functions
using ..la_utils   					      # library of linear algebra utilities 

### User-defined algorithm settings

	const kmax_pipg = 1000			   	  # max PIPG    iterations
	const kmax_powiter = 50				  # max power   iterations

	const ϵ_abs = 1e-10					  # small       absolute tolerance
	const ϵ_powiter = 1e-2				  # power iter. exit tolerance							(relative)
	const ϵ_primal = 1e-6				  # PIPG primal exit tolerance							(relative)
	const ϵ_dual = 1e-6		     		  # PIPG dual   exit tolerance							(relative)
	const ϵ_pd_JuMP = 1e-9				  # JuMP solver primal and dual feasibility tolerance	
	const ϵ_gap_JuMP = 1e-9 			  # JuMP solver primal and dual objective gap tolerance

	const ϵ_admm = 1e-6 				  # admm    	exit tolerance
	const ρ_admm = 2.0					  # admm 		step-size

	const solver_JuMP = :gurobi 			  # choice of solver in JuMP (gurobi,mosek,ecos)
	const err_type = :linf  			  # type of norm for computing error (l2,linf)

### User-defined problem parameters

	const nx = 6						  # no. of states
	const nu = 6 						  # no. of inputs
	const N = 3   					      # no. of discretization points
	const Δ = 0.9						  # sampling time (s)
	const tvec = [(j-1)*Δ for j in 1:N]   # time vector (s)

	# system parameters
	const a_blk = 0.1					  							# half the width of the block (m)
	const r_blk = 2.0					  							# radius of arc traversed by the block in the reference solution (m)
	const m = 2.0 						  							# mass of block (fixed) (kg)
	const ge = 9.806 					  							# magnitude of acceleration due to gravity at Earth sea-level (m s^2)

	const s1 = SVector{3}([-a_blk,0.0,0.0]) 						# position vector first contact point
	const s2 = SVector{3}([a_blk,0.0,0.0]) 							#                 second 
	
	const μ1 = tand(55)					      						# coefficient of friction associated with first contact point
	const μ2 = tand(55)					      						# 										  second
	
	const F1max = 100.0												# upper bound on magnitude of first contact force
	const F2max = 100.0												#							  second	
	
	const vmax = 100.0 												# upper bound on magnitude of velocity 

	# matrix in linear equality constraints due to rotational equilibrium
	const Γrot = [la_utils.skew_mat(s1) la_utils.skew_mat(s2)]  
	# projection matrix for ker Γrot
	const proj_ker_Γrot = la_utils.proj_ker(Γrot,nu)

	# reference trajectory and boundary conditions
	# τfun(t) = r_blk*(2*t/(N-1)-1)									# convenient function for reference trajector definition
	# const y_unscl = [SVector{nx}([τfun(j-1),2*τfun(j-1),sqrt(r_blk*r_blk-τfun(j-1)^2),0,0,0]) for j in 1:N]
	const y_unscl = [SVector{nx}(zeros(nx)) for _ in 1:N]
	const x0_unscl = SVector{nx}(y_unscl[1])
	const xf_unscl = SVector{nx}(y_unscl[N])
	
	# const x0_unscl = SVector{nx}(zeros(nx))
	# const xf_unscl = SVector{nx}([10.0,0.0,0.0,0.0,0.0,0.0])
	# const y_unscl = [SVector{nx}( x0_unscl .* (1-t/(N-1)) .+ xf_unscl .* (t/(N-1)) ) for t in 0:(N-1)]

	# state segments and scaling
		# segment 1 (position)
		const ix1 = 1:3						  			# segment indices
		const nx1 = length(ix1)							# length of segment
		const scl_x1 = 1.0							 	# scaling
		# const scl_x1 = 0.5*norm(x0_unscl[1:3])

		# segemnt 2 (velocity)
		const ix2 = 4:6 				  				# segment indices
		const nx2 = length(ix2)			  				# length of segment	
		const scl_x2 = 1.0 							 	# scaling
		# const scl_x2 = 0.5*norm(x0_unscl[4:6])

	# control segments and scaling
		# segment 1 (finger 1)
		const iu1 = 1:3									# segment indices
		const nu1 = length(iu1)							# length of segment
		const scl_u1 = 1.0 								# scaling

		# segment 2 (finger 2)
		const iu2 = 4:6									# segment indices
		const nu2 = length(iu2)							# length of segment
		const scl_u2 = 1.0 								# scaling

	# scaling matrices
		const scl_x_mat = Array(Diagonal(cat( scl_x1 .* ones(nx1) , scl_x2 .* ones(nx2) , dims=1 )))
		const scl_x_imat = Array(Diagonal(cat( (1/scl_x1) .* ones(nx1) , (1/scl_x2) .* ones(nx2) , dims=1 )))
		const scl_u_mat = Array(Diagonal(cat( scl_u1 .* ones(nu1) , scl_u2 .* ones(nu2) , dims=1 )))
		const scl_u_imat = Array(Diagonal(cat( (1/scl_u1) .* ones(nu1) , (1/scl_u2) .* ones(nu2) , dims=1 )))

	# ZOH discretization: x(t+1) = Ad*x(t) + Bd*u(t) + gd
		# unscaled
		const Ad_unscl = SMatrix{nx,nx}([1.0 0.0 0.0 Δ   0.0 0.0
								   		 0.0 1.0 0.0 0.0 Δ   0.0
								   		 0.0 0.0 1.0 0.0 0.0 Δ  
								   		 0.0 0.0 0.0 1.0 0.0 0.0
								   		 0.0 0.0 0.0 0.0 1.0 0.0
								  		 0.0 0.0 0.0 0.0 0.0 1.0])
		const Bd_unscl = SMatrix{nx,nu}([0.5*(Δ^2)/m  0.0          0.0			0.5*(Δ^2)/m  0.0          0.0		
								  		 0.0          0.5*(Δ^2)/m  0.0			0.0          0.5*(Δ^2)/m  0.0		
								  		 0.0          0.0          0.5*(Δ^2)/m  0.0          0.0          0.5*(Δ^2)/m  
								  		 Δ/m          0.0          0.0			Δ/m          0.0          0.0 		
								  		 0.0          Δ/m          0.0			0.0          Δ/m          0.0		
								  		 0.0          0.0          Δ/m 			0.0          0.0          Δ/m])
		const gd_unscl = SVector{nx}([0.0, 0.0, -0.5*(Δ^2)*ge, 0.0, 0.0, -Δ*ge])
		# const gd_unscl = SVector{nx}(zeros(nx))

		# scaled
		const Ad = SMatrix{nx,nx}(scl_x_imat*Ad_unscl*scl_x_mat)
		const Bd = SMatrix{nx,nu}(scl_x_imat*Bd_unscl*scl_u_mat)
		const gd = SVector{nx}(scl_x_imat*gd_unscl)

		const mAdT = SMatrix{nx,nx}(-1 .* Array(transpose(Ad)))
		const mBdT = SMatrix{nu,nx}(-1 .* Array(transpose(Bd)))

		const x0 = SVector{nx}(scl_x_imat*x0_unscl)
		const xf = SVector{nx}(scl_x_imat*xf_unscl)
		const y = [SVector{nx}(scl_x_imat*y_unscl[t]) for t = 1:N] 

	# LQR weighting matrices
	const Q = SMatrix{nx,nx}(scl_x_mat*Array(Diagonal( cat( 1.0 .* ones(nx1) , 1.0 .* ones(nx2) , dims=1 ) ))*scl_x_mat)
	const R = SMatrix{nu,nu}(scl_u_mat*Array(Diagonal( cat( 1.0 .* ones(nu1) , 1.0 .* ones(nu2) , dims=1 ) ))*scl_u_mat)
	const Qf = SMatrix{nx,nx}(scl_x_mat*Array(Diagonal( cat( 1.0 .* ones(nx1) , 1.0 .* ones(nx2) , dims=1 ) ))*scl_x_mat)

	# step-size parameters (min and max eigenvalues of LQR weight matrices)
	const μ = min(min([R[i,i] for i in 1:nu]...),min([Q[i,i] for i in 1:nx]...),min([Qf[i,i] for i in 1:nx]...))
	const λ = max(max([R[i,i] for i in 1:nu]...),max([Q[i,i] for i in 1:nx]...),max([Qf[i,i] for i in 1:nx]...))

	@assert μ>0 && (λ/μ)<10 "Condition number of objective Hessian is too larger: $(λ/μ). Modify Q, R and Qf"

	# constraints

		# friction cone on three segmens of u (soc)
		const α_μ1 = copy(μ1)
		const α_μ2 = copy(μ2)
		# upper bound on magnitude of contact force due to the three fingers (l2-ball)
		const r_F1 = F1max/scl_u1
		const r_F2 = F2max/scl_u2

		# rotational equilibrium (hyperplane)
		# Γrot and the projection matrix mapping to its kernel already defined

		# upper bound on velocity magnitude (l2-ball)
		const r_vmax = vmax/scl_x2

### Specifications for projection on constraint sets
	# projection counter
	const proj_count_limit = SVector{2}([1,2])			# asymptotic projection estimate for x and u	
	const proj_counter = MVector{2}(randn(2))			# [proj_count_x, proj_count_u]

	# inferred admm step-size
	const ν_admm = 1/(ρ_admm+1)

	# temp storage variables
	const γ1 = randn()
	const γ2 = randn()

	const κ1_x2 = MVector{nx2}(randn(nx2))
	const κ2_x2 = MVector{nx2}(randn(nx2))
	const κ1_u1 = MVector{nu1}(randn(nu1))
	const κ2_u1 = MVector{nu1}(randn(nu1))
	const κ1_u2 = MVector{nu2}(randn(nu2))
	const κ2_u2 = MVector{nu2}(randn(nu2))
	const κ1_u = MVector{nu}(randn(nu))
	const κ2_u = MVector{nu}(randn(nu))					# first element stores no. of projections performed in the last call to admm!
	const y_u = [MVector{nu}(randn(nu)) for _ in 1:N-1] # container for admm warm-start
	const u_u = [MVector{nu}(randn(nu)) for _ in 1:N-1] # container for admm warm-start

	function project_x!(x,z,t)
		# no projection operation on segment 1
		BLAS.blascopy!(nx1,z,1,x,1)

		# projection on segment 2
		for j in ix2
			κ2_x2[j-nx1] = z[j]
		end
		proj.ball!(κ1_x2,κ2_x2,r_vmax,nx2)
		for j in ix2
			x[j] = κ1_x2[j-nx1]
		end
	end

	function proj_friction_soc_ball!(u,w)
		κ2_u1[1] = w[2]
		κ2_u1[2] = w[3]
		κ2_u1[3] = w[1]
		proj.soc_ball!(κ1_u1,κ2_u1,α_μ1,r_F1,nu1,γ1)
		u[1] = κ1_u1[3]
		u[2] = κ1_u1[1]
		u[3] = κ1_u1[2]

		κ2_u2[1] = w[5]
		κ2_u2[2] = w[6]
		κ2_u2[3] = -w[4]
		proj.soc_ball!(κ1_u2,κ2_u2,α_μ2,r_F2,nu2,γ1)
		u[4] = -κ1_u2[3]
		u[5] = κ1_u2[1]
		u[6] = κ1_u2[2]
	end

	const fu1(u,z) = proj_friction_soc_ball!(u,z)
	const fu2(u,z) = mul!(u,proj_ker_Γrot,z)
	const project_u!(u,z,t) = proj.admm!(u,y_u[t],u_u[t],z,fu1,fu2,ρ_admm,ν_admm,ϵ_admm,nu,γ1,γ2,κ1_u,κ2_u)

	# diagnostics projection functions
	function project_x_diagnostic!(x,z,t)
		project_x!(x,z,t)
		proj_counter[1] = 1
	end
	function project_u_diagnostic!(u,z,t)
		project_u!(u,z,t)
		proj_counter[2] = κ2_u[1]
	end

	# specify constraints in a JuMP problem
	function set_constr_JuMP!(model,x,u)			
		# friction cone 1
		@constraint(model,[t=1:N-1],[α_μ1*u[1,t],u[2:3,t]...] in SecondOrderCone())

		# friction cone 2
		@constraint(model,[t=1:N-1],[-α_μ2*u[4,t],u[5:6,t]...] in SecondOrderCone())

		# contact force 1 upper bound
		@constraint(model,[t=1:N-1],[r_F1,u[1:3,t]...] in SecondOrderCone())

		# contact force 2 upper bound
		@constraint(model,[t=1:N-1],[r_F2,u[4:6,t]...] in SecondOrderCone())

		# rotational equilibrium
		@constraint(model,[t=1:N-1],Γrot*u[:,t] .== 0)

		# velocity magnitude upper bound
		@constraint(model,[t=1:N-1],[r_vmax,x[4:6,t+1]...] in SecondOrderCone())
	end
end	