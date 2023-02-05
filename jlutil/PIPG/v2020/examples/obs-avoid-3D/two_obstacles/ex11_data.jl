# 03/01/21
# quadrotor (3D double integrator) path planning with 2 obstacles
# keep-out-zone constraint convexified using rotating halfspaces
# thrust pointing constraints

const path_to_tools = "../../pipg_toolkit/"
macro load_file(str_val,flg) return :( $flg ? include(string(path_to_tools,$str_val)) : include($str_val) ) end

@load_file "proj_funcs.jl" true			  # projection functions library

# module defining the parameters for example 11 problem
module eX 
using LinearAlgebra, StaticArrays         # required packages
using JuMP
using ..proj 							  # library of projection functions

### User-defined algorithm settings

	const kmax_pipg = 10000			   	  # max PIPG    iterations
	const kmax_powiter = 50				  # max power   iterations

	const ϵ_abs = 1e-10					  # small       absolute tolerance
	const ϵ_powiter = 1e-2				  # power iter. exit tolerance							(relative)
	const ϵ_pd_JuMP = 1e-9				  # JuMP solver primal and dual feasibility tolerance	
	const ϵ_gap_JuMP = 1e-9 			  # JuMP solver primal and dual objective gap tolerance

	const ϵ_primal = 1e-4				  # PIPG primal exit tolerance							(relative)
	const ϵ_dual = 1e-4		     		  # PIPG dual   exit tolerance							(relative)
	const ϵ_admm = 1e-6 				  # admm    	exit tolerance
	
	const ρ_admm = 2.0					  # admm 		step-size

	const solver_JuMP = :mosek 			  # choice of solver in JuMP (gurobi,mosek,ecos)
	const err_type = :linf  			  # type of norm for computing error (l2,linf)

### User-defined problem parameters

	const nx = 6						  # no. of states
	const nu = 3 						  # no. of inputs
	const N = 25   					      # no. of discretization points

	const Δ = 0.5						  # sampling time (s)
	const tvec = [(j-1)*Δ for j in 1:N]   # time vector (s)

	# system parameters
	const m = 3.0 						      # mass of quadrotor (fixed) (kg)
	const ge = 9.806 					      # magnitude of acceleration due to gravity at Earth sea-level (m s^-2)	

	const vmax = 7.0						  # velocity magnitude upper bound (m s^-1)
	const umax = 30.0 					      # acceleration magnitude upper bound (m s^-2)
	const θ_tl = 35*π/180 					  # maximum thrust pointing angle (rad)

	const p1 = SVector{3}([-10.0,-10.0,0.0])  # center of obstacle #1 
	const p2 = SVector{3}([8.0,-8.0,0.0])	  # center of obstacle #2	
	const p3 = SVector{3}([10.0,10.0,0.0])	  # center of obstacle #3

	const r1 = 7.6						  	  # radius of obstacle #1
	const r2 = 9.0						      # radius of obstacle #2
	const r3 = 7.6                            # radius of obstacle #3 

	# boundary conditions
	const x0_unscl = SVector{nx}(cat(p1 .- [0.0,1.2*r1,0.0],[0.0,0.0,0.0],dims=1))
	const xf_unscl = SVector{nx}(cat(p3 .+ [1.2*r3,0.0,0.0],[0.0,0.0,0.0],dims=1))

	# reference trajectory
	const y_unscl = [SVector{nx}((1-t/(N-1)) .* x0_unscl .+ (t/(N-1)) .* xf_unscl) for t=0:(N-1)]


	# state segments and scaling
		# segment 1 (position)
		const ix1 = 1:3						  			# segment indices
		const nx1 = length(ix1)							# length of segment
		# const scl_x1 = 1.0					 		# scaling
		const scl_x1 = max(r1,r2,r3)

		# segemnt 2 (velocity)
		const ix2 = 4:6 				  				# segment indices
		const nx2 = length(ix2)			  				# length of segment	
		# const scl_x2 = 1.0					 		# scaling
		const scl_x2 = copy(vmax)

	# control segments and scaling
		# no segments
		# const scl_u = 1.0	 							# scaling
		const scl_u = copy(umax)

	# scaling matrices (inferred from scl_x1, scl_x2, ... , scl_u)
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
	const Q = SMatrix{nx,nx}(scl_x_mat*Array(Diagonal( cat( (1.0/scl_x1^2) .* ones(nx1) , (1.0/scl_x2^2) .* ones(nx2) , dims=1 ) ))*scl_x_mat)
	const R = SMatrix{nu,nu}(scl_u_mat*Array(Diagonal( (1.0/scl_u^2) .* ones(nu) ))*scl_u_mat)
	const Qf = SMatrix{nx,nx}(scl_x_mat*Array(Diagonal( cat( (10.0/scl_x1^2) .* ones(nx1) , (10.0/scl_x2^2) .* ones(nx2) , dims=1 ) ))*scl_x_mat) 

	# step-size parameters (min and max eigenvalues of LQR weight matrices)
	const μ = min(min([R[i,i] for i in 1:nu]...),min([Q[i,i] for i in 1:nx]...),min([Qf[i,i] for i in 1:nx]...))
	const λ = max(max([R[i,i] for i in 1:nu]...),max([Q[i,i] for i in 1:nx]...),max([Qf[i,i] for i in 1:nx]...))

	@assert μ>0 && (λ/μ)<=50 "Condition number of objective Hessian is too larger: $(λ/μ). Modify Q, R and Qf"

	# constraints

		# obstacle avoidance

			# linearization
				# initial and final point for obstacle linearization reference
				const ȳ0_unscl = SVector{nx}(cat(p1 .- [0.0,r1,0.0],[0.0,0.0,0.0],dims=1))
				const ȳf_unscl = SVector{nx}(cat(p3 .+ [r3,0.0,0.0],[0.0,0.0,0.0],dims=1))
				# obstructed straight-line reference for obstacle linearization
				const ȳ_obs_unscl = [SVector{nx}( (1-t/(N-1)) .* ȳ0_unscl .+ (t/(N-1)) .* ȳf_unscl ) for t=0:(N-1)]

				# move state to nearest obstacle boundary if penetrated
				function path_unobstructed(z)
					# globals used p1, p2, p3
					yy = zeros(nx)
					yy .= z

					tol_fac = 1.0

					# nrm_y_p1 = norm(yy[1:3]-p1)
					# if nrm_y_p1 ≤ tol_fac*r1
					# 	yy[1:3] .= (yy[1:3] .- p1) .* (tol_fac*r1/nrm_y_p1) .+ p1 	
					# end

					nrm_y_p2 = norm(yy[1:3] .- p2)
					if nrm_y_p2 ≤ tol_fac*r2
						yy[1:3] .= (yy[1:3] .- p2) .* (tol_fac*r2/nrm_y_p2) .+ p2 	
					end

					nrm_y_p3 = norm(yy[1:3]-p3)
					if nrm_y_p3 ≤ tol_fac*r3
						yy[1:3] .= (yy[1:3] .- p3) .* (tol_fac*r3/nrm_y_p3) .+ p3 	
					end
					return yy
				end

				# (scaled) unobstructed straight-line reference for obstacle linearization
				const ȳ_unobs = [SVector{nx}(scl_x_imat*(path_unobstructed(ȳ_obs_unscl[t]))) for t=1:N]

				# (scaled) terms in the obstacle avoidance constraint linearization
				#  - ϕᵢ⋅x + θᵢ ≈ norm(x[1:3]-pᵢ) - rᵢ ≥ 0
				# const θ1 = SVector{N-1}([norm(ȳ_unobs[t+1][1:3] .- (p1 ./ scl_x1)) - r1/scl_x1 + (-1 * norm(ȳ_unobs[t+1][1:3])^2 + (p1 ./ scl_x1)⋅ȳ_unobs[t+1][1:3]) / norm(ȳ_unobs[t+1][1:3] .- (p1 ./ scl_x1)) for t=1:(N-1)])
				const θ2 = SVector{N-1}([norm(ȳ_unobs[t+1][1:3] .- (p2 ./ scl_x1)) - r2/scl_x1 + (-1 * norm(ȳ_unobs[t+1][1:3])^2 + (p2 ./ scl_x1)⋅ȳ_unobs[t+1][1:3]) / norm(ȳ_unobs[t+1][1:3] .- (p2 ./ scl_x1)) for t=1:(N-1)])
				const θ3 = SVector{N-1}([norm(ȳ_unobs[t+1][1:3] .- (p3 ./ scl_x1)) - r3/scl_x1 + (-1 * norm(ȳ_unobs[t+1][1:3])^2 + (p3 ./ scl_x1)⋅ȳ_unobs[t+1][1:3]) / norm(ȳ_unobs[t+1][1:3] .- (p3 ./ scl_x1)) for t=1:(N-1)])			

				# const ϕ1 = [SVector{nx1}(-1 .* (ȳ_unobs[t+1][1:3] .- (p1 ./ scl_x1)) ./ norm(ȳ_unobs[t+1][1:3] .- (p1 ./ scl_x1)) ) for t=1:N-1]
				const ϕ2 = [SVector{nx1}(-1 .* (ȳ_unobs[t+1][1:3] .- (p2 ./ scl_x1)) ./ norm(ȳ_unobs[t+1][1:3] .- (p2 ./ scl_x1)) ) for t=1:N-1]
				const ϕ3 = [SVector{nx1}(-1 .* (ȳ_unobs[t+1][1:3] .- (p3 ./ scl_x1)) ./ norm(ȳ_unobs[t+1][1:3] .- (p3 ./ scl_x1)) ) for t=1:N-1]

			# half space parameters
			# const u1_obs = deepcopy(ϕ1)
			# const ζ1_obs = deepcopy(θ1)

			const u2_obs = deepcopy(ϕ2)
			const ζ2_obs = deepcopy(θ2)

			const u3_obs = deepcopy(ϕ3)
			const ζ3_obs = deepcopy(θ3)
					
		# l2-ball parameters (velocity magnitude upper bound)
		const r_vmax = vmax/scl_x2

		# l2-ball parameters (thrust magnitude upper bound)
		const r_umax = umax/scl_u

		# soc parameters (thrust poining cone) 
		const α_tl = tan(θ_tl)


### Specifications for projection on constraint sets
	# projection counter
	const proj_count_limit = SVector{2}([2,1])			# asymptotic projection estimate for x and u
	const proj_counter = MVector{2}(randn(2))			# [proj_count_x, proj_count_u]

	# inferred admm step-size
	const ν_admm = 1/(ρ_admm+1)

	# temp storage variables
	const γ1 = randn()
	const γ2 = randn()
	const γ3 = randn()

	const κ1_x1 = MVector{nx1}(randn(nx1))
	const κ2_x1 = MVector{nx1}(randn(nx1))
	const κ3_x1 = MVector{nx1}(randn(nx1))
	const κ4_x1 = MVector{nx1}(randn(nx1))

	const κ1_x2 = MVector{nx2}(randn(nx2))
	const κ2_x2 = MVector{nx2}(randn(nx2))
	
	const κ1_u = MVector{nu}(randn(nu))
	const κ2_u = MVector{nu}(randn(nu))
	
	const u_y_admm = [nothing for _ in 1:N]
	const u_u_admm = [nothing for _ in 1:N]
	const x_y_admm = [nothing for _ in 1:N]
	const x_u_admm = [nothing for _ in 1:N]

	function project_x!(x,z,t)
		# projection on segment 1
		BLAS.blascopy!(nx1,z,1,κ2_x1,1)
		proj.hlfspace2!(κ1_x1,κ2_x1,u2_obs[t],u3_obs[t],ζ2_obs[t],ζ3_obs[t],nx1,γ1,γ2,γ3)
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

	function project_u!(u,z,t)
		proj.soc_ball!(u,z,α_tl,r_umax,nu,γ1)
	end

	# diagnostics projection function
	function project_x_diagnostic!(x,z,t)
		project_x!(x,z,t)
		proj_counter[1] = 2
	end

	function project_u_diagnostic!(u,z,t)
		project_u!(u,z,t)
		proj_counter[2] = 1
	end 

	# specify constraints in a JuMP problem
	function set_constr_JuMP!(model,x,u)			
		# acceleration upper bound
		@constraint(model,[t=1:N-1],[r_umax,u[:,t]...] in SecondOrderCone()) 

		# obstacle avoidance constraints constraint 
		# @constraint(model,[t=1:N-1],u1_obs[t]⋅x[1:3,t+1]-ζ1_obs[t] ≤ 0)
		@constraint(model,[t=1:N-1],u2_obs[t]⋅x[1:3,t+1]-ζ2_obs[t] ≤ 0)
		@constraint(model,[t=1:N-1],u3_obs[t]⋅x[1:3,t+1]-ζ3_obs[t] ≤ 0)

		# thrust pointing constraint
		@constraint(model,[t=1:N-1],[α_tl*u[3,t],u[1:2,t]...] in SecondOrderCone())

		# velocity magnitude upper bound
		@constraint(model,[t=1:N-1],[r_vmax,x[4:6,t+1]...] in SecondOrderCone())
	end
end	