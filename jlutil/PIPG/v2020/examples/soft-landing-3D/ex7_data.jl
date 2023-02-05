# Rocket soft landing
# Thrust lower bound is approximated using halfspace
# Mass dynamics is excluded
# Convex thrust pointing constraint is included (i.e. θ ∈ (0,π/2))
# Convec glide-slope constraints is included (i.e. γ ∈ (0,π/2))

const path_to_tools = "../pipg_toolkit/"
macro load_file(str_val,flg) return :( $flg ? include(string(path_to_tools,$str_val)) : include($str_val) ) end

@load_file "proj_funcs.jl" true			  # projection functions library	

# module defining the parameters for example 7 problem
module eX 
using LinearAlgebra, StaticArrays         # required packages
using JuMP
using ..proj 							  # library of projection functions

### User-defined algorithm settings

	const kmax_pipg = 5000			   	  # max PIPG    iterations
	const kmax_powiter = 50				  # max power   iterations

	const ϵ_abs = 1e-10					  # small       absolute tolerance
	const ϵ_powiter = 1e-2				  # power iter. exit tolerance							(relative)

	if isdefined(Main,:load_termination_criteria)
		const ϵ_pd_JuMP = copy(Main.ϵ_pd_JuMP) 
		const ϵ_gap_JuMP = copy(Main.ϵ_gap_JuMP) 
		const ϵ_primal = copy(Main.ϵ_primal)
		const ϵ_dual = copy(Main.ϵ_dual)
		const ϵ_admm = copy(Main.ϵ_admm)
	else
		const ϵ_pd_JuMP = 1e-9				  # JuMP solver primal and dual feasibility tolerance	
		const ϵ_gap_JuMP = 1e-9 			  # JuMP solver primal and dual objective gap tolerance
		const ϵ_primal = 1e-4				  # PIPG primal exit tolerance							(relative)
		const ϵ_dual = 1e-4		     		  # PIPG dual   exit tolerance							(relative)
		const ϵ_admm = 1e-6 				  # admm    	exit tolerance
	end
	
	const ρ_admm = 2.2					  # admm 		step-size

	const solver_JuMP = :mosek 			  # choice of solver in JuMP (gurobi,mosek,ecos)
	const err_type = :linf  			  # type of norm for computing error (l2,linf)


### User-defined problem parameters

	const nx = 6						  # no. of states
	const nu = 3 						  # no. of inputs
	const N = 20   					      # no. of discretization points
	const T = 10.0  			 		  # duration of flight (s)

	const Δ = T/(N-1)					  # sampling time (s)
	const tvec = [(j-1)*Δ for j in 1:N]   # time vector (s)

	# system parameters
	const m = 10.0 						  # mass of lander (fixed) (kg)
	const ge = 9.806 					  # magnitude of acceleration due to gravity at Earth sea-level (m s^2)

	# boundary conditions
	const x0_unscl =SVector{nx}([4.4, 13.5, 23.0, -3.0, 2.0, -5.0])
	const xf_unscl = SVector{nx}(zeros(6))

	# reference trajectory
	const y_unscl = [SVector{nx}((1-t/(N-1)) .* x0_unscl .+ (t/(N-1)) .* xf_unscl) for t=0:(N-1)]


	# state segments and scaling
		# segment 1 (position)
		const ix1 = 1:3						  			# segment indices
		const nx1 = length(ix1)							# length of segment
		const scl_x1 = 0.5*norm(x0_unscl[1:3]) 			# scaling
	
		# segemnt 2 (velocity)
		const ix2 = 4:6 				  				# segment indices
		const nx2 = length(ix2)			  				# length of segment	
		const scl_x2 = 0.5*norm(x0_unscl[4:6]) 			# scaling

	# control segments and scaling
		# no segments
		const scl_u = 0.8*m*ge 							# scaling

	# scaling matrices
		const scl_x_mat = Array(Diagonal(cat( scl_x1 .* ones(nx1) , scl_x2 .* ones(nx2) , dims=1 )))
		const scl_x_imat = Array(Diagonal(cat( (1/scl_x1) .* ones(nx1) , (1/scl_x2) .* ones(nx2) , dims=1 )))
		const scl_u_mat = Array(Diagonal( scl_u .* ones(nu) ))
		const scl_u_imat = Array(Diagonal( (1/scl_u) .* ones(nu) ))

	# ZOH discretization: x(t+1) = Ad*x(t) + Bd*u(t) + gd
		# unscaled
		const Ad_unscl = SMatrix{nx,nx}([1.0 0.0 0.0 Δ   0.0 0.0
								   		 0.0 1.0 0.0 0.0 Δ   0.0
								   		 0.0 0.0 1.0 0.0 0.0 Δ  
								   		 0.0 0.0 0.0 1.0 0.0 0.0
								   		 0.0 0.0 0.0 0.0 1.0 0.0
								  		 0.0 0.0 0.0 0.0 0.0 1.0])
		const Bd_unscl = SMatrix{nx,nu}([0.5*(Δ^2)/m  0.0          0.0
								  		 0.0          0.5*(Δ^2)/m  0.0
								  		 0.0          0.0          0.5*(Δ^2)/m
								  		 Δ/m          0.0          0.0
								  		 0.0          Δ/m          0.0
								  		 0.0          0.0          Δ/m])

		const gd_unscl = SVector{nx}([0.0, 0.0, -0.5*(Δ^2)*ge, 0.0, 0.0, -Δ*ge])

		# scaled
		const Ad = SMatrix{nx,nx}(scl_x_imat*Ad_unscl*scl_x_mat)
		const Bd = SMatrix{nx,nu}(scl_x_imat*Bd_unscl*scl_u_mat)
		const gd = SVector{nx}(scl_x_imat*gd_unscl)

		const mAdT = SMatrix{nx,nx}(-1 .* Array(transpose(Ad)))
		const mBdT = SMatrix{nu,nx}(-1 .* Array(transpose(Bd)))

		const x0 = SVector{nx}(scl_x_imat*x0_unscl)
		const xf = SVector{nx}(scl_x_imat*xf_unscl)
		const y = [SVector{nx}(scl_x_imat*y_unscl[t]) for t = 1:N] 
		const uref = nothing

	# LQR weighting matrices
	const Q = SMatrix{nx,nx}(scl_x_mat*Array(Diagonal( cat( 0.1 .* ones(nx1) , 1.0 .* ones(nx2) , dims=1 ) ))*scl_x_mat)
	const R = SMatrix{nu,nu}(scl_u_mat*Array(Diagonal( 0.005 .* ones(nu) ))*scl_u_mat)
	const Qf = SMatrix{nx,nx}(scl_x_mat*Array(Diagonal( cat( 0.5 .* ones(nx1) , 5.0 .* ones(nx2) , dims=1 ) ))*scl_x_mat)

	# step-size parameters (min and max eigenvalues of LQR weight matrices)
	const μ = min(min([R[i,i] for i in 1:nu]...),min([Q[i,i] for i in 1:nx]...),min([Qf[i,i] for i in 1:nx]...))
	const λ = max(max([R[i,i] for i in 1:nu]...),max([Q[i,i] for i in 1:nx]...),max([Qf[i,i] for i in 1:nx]...))

	@assert μ>0 && (λ/μ)<10 "Condition number of objective Hessian is too larger: $(λ/μ). Modify Q, R and Qf"

	# constraints

		# glide-slope cone on position (assuming xf[1:3] = 0)  
		const γ_gs = 50.0*π/180							# maximum glide-slope angle
		
			# soc type 1 parameters
			const α_gs = tan(γ_gs)

		# upper bound on velocity magnitude
		const vmax = 9.0/scl_x2							# upper bound on speed

			# l2-ball parameters
			const r_vmax = copy(vmax)

		# thrust poining cone 
		const θ_tl = 9.0*π/180							# maximum pointing angle

			# soc type 1 parameters
			const α_tl = tan(θ_tl)

		# thrust upper bound
		const umax = 1.1*m*ge/scl_u 					# thrust magnitude upper bound

			# l2-ball parameters
			const r_umax = copy(umax)

		# thrust lower bound
		const umin = 0.6*m*ge/scl_u 					# thrust third component lower bound

			# half space parameters
			const u0_umin = SVector{nu}([0.0,0.0,-1.0])
			const ζ0_umin = copy(-umin)

### Specifications for projection on constraint sets
	# projection counter
	const proj_count_limit = SVector{2}([2,3])			# asymptotic projection estimate for x and u	
	const proj_counter = MVector{2}(randn(2))			# [proj_count_x, proj_count_u]

	# inferred admm step-size
	const ν_admm = 1/(ρ_admm+1)

	# temp storage variables
	const γ1 = randn()
	const γ2 = randn()
	
	const κ1_x1 = MVector{nx1}(randn(nx1))
	const κ2_x1 = MVector{nx1}(randn(nx1))
	const κ1_x2 = MVector{nx2}(randn(nx2))
	const κ2_x2 = MVector{nx2}(randn(nx2))
	const κ1_u = MVector{nu}(randn(nu))
	const κ2_u = MVector{nu}(randn(nu))			 			 # first element stores no. of projections performed in the last call to admm!
	const u_y_admm = [MVector{nu}(randn(nu)) for _ in 1:N-1] # container for admm warm-start
	const u_u_admm = [MVector{nu}(randn(nu)) for _ in 1:N-1] # container for admm warm-start
	const x_y_admm = [nothing for _ in 1:N-1]
	const x_u_admm = [nothing for _ in 1:N-1]

	function project_x!(x,z,t)
		# projection on segment 1
		BLAS.blascopy!(nx1,z,1,κ2_x1,1)
		proj.soc!(κ1_x1,κ2_x1,α_gs,nx1,γ1)
		BLAS.blascopy!(nx1,κ1_x1,1,x,1)

		# projection on segment 2
		for j in ix2
			κ2_x2[j-nx1] = z[j]
		end
		proj.ball!(κ1_x2,κ2_x2,r_vmax,nx2)
		for j in ix2
			x[j] = κ1_x2[j-nx1]
		end
	end

	const fu1(u,z) = proj.soc_ball!(u,z,α_tl,r_umax,nu,γ1)
	const fu2(u,z) = proj.hlfspace!(u,z,u0_umin,ζ0_umin,nu,γ1)
	const project_u!(u,z,t) = proj.admm!(u,u_y_admm[t],u_u_admm[t],z,fu1,fu2,ρ_admm,ν_admm,ϵ_admm,nu,γ1,γ2,κ1_u,κ2_u)

	# diagnostics projection functions
	function project_x_diagnostic!(x,z,t)
		project_x!(x,z,t)
		proj_counter[1] = 2
	end

	function project_u_diagnostic!(u,z,t)
		project_u!(u,z,t)
		proj_counter[2] = κ2_u[1]+1
	end

	# specify constraints in a JuMP problem
	function set_constr_JuMP!(model,x,u)			
		# thrust upper bound
		@constraint(model,[t=1:N-1],[r_umax,u[:,t]...] in SecondOrderCone()) 

		# thrust lower bound
		@constraint(model,[t=1:N-1],u[3,t] ≥ umin)

		# thrust pointing constraint
		@constraint(model,[t=1:N-1],[α_tl*u[3,t],u[1:2,t]...] in SecondOrderCone())

		# glide-slope constraint 
		@constraint(model,[t=1:N-1],[α_gs*x[3,t+1],x[1:2,t+1]...] in SecondOrderCone())

		# velocity magnitude upper bound
		@constraint(model,[t=1:N-1],[r_vmax,x[4:6,t+1]...] in SecondOrderCone())
	end
end	