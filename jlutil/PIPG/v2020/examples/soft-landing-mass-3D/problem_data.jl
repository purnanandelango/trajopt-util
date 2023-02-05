# Rocket soft landing 
# Thrust lower bound is handled via LCvx
# Mass dynamics is included
# Convex glide-slope constraints is included (i.e. γ ∈ (0,π/2))
# Thrust pointing constraint is excluded

const path_to_tools = "../../toolkit/"
macro load_file(str_val,flg) return :( $flg ? include(string(path_to_tools,$str_val)) : include($str_val) ) end

@load_file "proj_funcs.jl" true			  # Projection functions library	

# Module defining the problem parameters
module eX 
using LinearAlgebra, StaticArrays         # Required packages
using JuMP
using ..proj 							  # Library of projection functions

### User-defined algorithm settings

	const kmax_pipg = 20000			   	  # Max PIPG    iterations
	const kmax_powiter = 50				  # Max power   iterations

	const ϵ_abs = 1e-10					  # Small       absolute tolerance
	const ϵ_powiter = 1e-2				  # Power iter. exit tolerance							(relative)

	if isdefined(Main,:load_termination_criteria)
		const ϵ_pd_JuMP = copy(Main.ϵ_pd_JuMP) 
		const ϵ_gap_JuMP = copy(Main.ϵ_gap_JuMP) 
		const ϵ_primal = copy(Main.ϵ_primal)
		const ϵ_dual = copy(Main.ϵ_dual)
		const ϵ_admm = copy(Main.ϵ_admm)
	else
		const ϵ_pd_JuMP = 1e-7				  # JuMP solver primal and dual feasibility tolerance	
		const ϵ_gap_JuMP = 1e-7 			  # JuMP solver primal and dual objective gap tolerance
		const ϵ_primal = 9e-4				  # PIPG primal exit tolerance							(relative)
		const ϵ_dual = 9e-4		     		  # PIPG dual   exit tolerance							(relative)
		const ϵ_admm = 1e-4 				  # ADMM    	exit tolerance
	end
	
	const ρ_admm = 2.2					  # ADMM 		step-size

	const solver_JuMP = :gurobi 		  # Choice of solver in JuMP (gurobi,mosek,ecos)
	const err_type = :linf  			  # Type of norm for computing error (l2,linf)


### User-defined problem parameters

	const nx = 7						  # No. of states
	const nu = 4 						  # No. of inputs
	const N = 25   					      # No. of discretization points
	const T = 80.0  			 		  # Duration of flight (s)

	const Δ = T/(N-1)					  # Sampling time (s)
	const tvec = [(j-1)*Δ for j in 1:N]   # Time vector (s)

	# System parameters
	const mwet = 2000.0
	const mdry = 150.0
	const Tmax = 24000.0 
	const ρ2 = 0.8*Tmax
	const ρ1 = 0.2*Tmax
	const vmax = 150.0
	const tanγ_gs = tan( 85 * π / 180 )
	const α_m = 0.0005
	const g_mars = 3.71 				  # Magnitude of acceleration due to gravity at Mars "sea"-level (m s^2)

	const z1 = [log(mwet-α_m*ρ1*tvec[t]) for t=1:N]
	const z0 = [log(mwet-α_m*ρ2*tvec[t]) for t=1:N]	
	const μ1 = [ρ1*exp(-z0[t]) for t=1:N]
	const μ2 = [ρ2*exp(-z0[t]) for t=1:N]	

	# Boundary conditions
	const x0_unscl =SVector{nx}([0.0, 2000.0, 1500.0, -0.0, 100.0, -75.0, log(mwet)])
	const xf_unscl = SVector{nx}([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, log(mdry)])

	# Reference trajectory
	const y_unscl = [SVector{nx}((1-t/(N-1)) .* x0_unscl .+ (t/(N-1)) .* xf_unscl) for t=0:(N-1)]


	# State segments and scaling
		# Segment 1 (position)
		const ix1 = 1:3						  			# Segment indices
		const nx1 = length(ix1)							# Length of segment
		const scl_x1 = 500.0
		# const scl_x1 = 0.5*norm(x0_unscl[1:3]) 		# Scaling
	
		# Segment 2 (velocity)
		const ix2 = 4:6 				  				# Segment indices
		const nx2 = length(ix2)			  				# Length of segment	
		# const scl_x2 = 0.5*norm(x0_unscl[4:6]) 		# Scaling
		const scl_x2 = 75.0

		# Segment 3 (log mass)
		const ix3 = 7									# Segment indices
		const nx3 = 1									# Length of segment
		const scl_x3 = 1.0								# Scaling

	# Control segments and scaling
		# Segments 1 (scaled thrust)
		const iu1 = 1:3									# Segment indices
		const nu1 = length(iu1)							# Length of segment			
		const scl_u1 = 1.0	 							# Scaling

		# Segments 2 (scaled thrust)
		const iu2 = 4									# Segment indices
		const nu2 = length(iu2)							# Length of segment			
		const scl_u2 = 10.0	 							# Scaling

	# Scaling matrices
		const scl_x_mat = Array(Diagonal(cat( scl_x1 .* ones(nx1) , scl_x2 .* ones(nx2) , scl_x3 .* ones(nx3) , dims=1 )))
		const scl_x_imat = Array(Diagonal(cat( (1/scl_x1) .* ones(nx1) , (1/scl_x2) .* ones(nx2) , (1/scl_x3) .* ones(nx3) , dims=1 )))
		const scl_u_mat = Array(Diagonal( cat( scl_u1 .* ones(nu1) , scl_u2 .* ones(nu2) , dims=1 )))
		const scl_u_imat = Array(Diagonal( cat( (1/scl_u1) .* ones(nu1) , (1/scl_u2) .* ones(nu2) , dims=1 )))

	# ZOH discretization: x(t+1) = Ad*x(t) + Bd*u(t) + gd
		# Unscaled
		const Ad_unscl = SMatrix{nx,nx}([1.0 0.0 0.0 Δ   0.0 0.0 0.0
								   		 0.0 1.0 0.0 0.0 Δ   0.0 0.0
								   		 0.0 0.0 1.0 0.0 0.0 Δ   0.0
								   		 0.0 0.0 0.0 1.0 0.0 0.0 0.0
								   		 0.0 0.0 0.0 0.0 1.0 0.0 0.0
								  		 0.0 0.0 0.0 0.0 0.0 1.0 0.0
								  		 0.0 0.0 0.0 0.0 0.0 0.0 1.0])
		const Bd_unscl = SMatrix{nx,nu}([0.5*(Δ^2)	  0.0          0.0		 	0.0
								  		 0.0          0.5*(Δ^2)	   0.0			0.0
								  		 0.0          0.0          0.5*(Δ^2)	0.0	
								  		 Δ            0.0          0.0			0.0
								  		 0.0          Δ            0.0			0.0
								  		 0.0          0.0          Δ            0.0
								  		 0.0 		  0.0 		   0.0		   -α_m*Δ])

		const gd_unscl = SVector{nx}([0.0, 0.0, -0.5*(Δ^2)*g_mars, 0.0, 0.0, -Δ*g_mars, 0.0])


		# Scaled
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
	const Q = SMatrix{nx,nx}(scl_x_mat*Array(Diagonal( cat( (10/scl_x1^2) .* ones(nx1) , (10/scl_x2^2) .* ones(nx2) , (1/scl_x3^2) .* ones(nx3) , dims=1 ) ))*scl_x_mat)
	const R = SMatrix{nu,nu}(scl_u_mat*Array(Diagonal( cat( (1/scl_u1^2) .* ones(nu1) , (1/scl_u2^2) .* ones(nu2) , dims=1 ) ))*scl_u_mat)
	const Qf = SMatrix{nx,nx}(scl_x_mat*Array(Diagonal( cat( (20/scl_x1^2) .* ones(nx1) , (20/scl_x2^2) .* ones(nx2) , (1/scl_x3^2) .* ones(nx3) , dims=1 ) ))*scl_x_mat)


	# Step-size parameters (min and max eigenvalues of LQR weight matrices)
	const μ = min(min([R[i,i] for i in 1:nu]...),min([Q[i,i] for i in 1:nx]...),min([Qf[i,i] for i in 1:nx]...))
	const λ = max(max([R[i,i] for i in 1:nu]...),max([Q[i,i] for i in 1:nx]...),max([Qf[i,i] for i in 1:nx]...))

	@assert μ>0 && (λ/μ)≤10000 "Condition number of objective Hessian is too larger: $(λ/μ). Modify Q, R and Qf"

	@show λ/μ

	# Constraints

		# Glide-slope cone on position (assuming xf[1:3] = 0)  
		# SOC type 1 parameters
		const α_gs = tanγ_gs

		const α_Tσ = scl_u2/scl_u1

		# Upper bound on velocity magnitude
		# l2-ball parameters
		const r_vmax = vmax/scl_x2

		# Half space parameters
		const u1_zσ = [SVector{2}([-μ1[t]*scl_x3,-1.0*scl_u2] ./ norm([-μ1[t]*scl_x3,-1.0*scl_u2])) for t in 1:N]
		const ζ1_zσ = [-μ1[t]*(1+z0[t]) / norm([-μ1[t]*scl_x3,-1.0*scl_u2]) for t in 1:N]

		const u2_zσ = [SVector{2}([μ2[t]*scl_x3,1.0*scl_u2] ./ norm([μ2[t]*scl_x3,1.0*scl_u2])) for t in 1:N]
		const ζ2_zσ = [μ2[t]*(1+z0[t]) / norm([μ2[t]*scl_x3,1.0*scl_u2]) for t in 1:N]

		# Box parameters
		const l_z = [z0[t]/scl_x3 for t=1:N]
		const u_z = [z1[t]/scl_x3 for t=1:N]

### Specifications for projection on constraint sets
	# Projection counter
	const proj_count_limit = SVector{2}([2,2])			# Asymptotic projection estimate for x and u	
	const proj_counter = MVector{2}(randn(2))			# [proj_count_x, proj_count_u]

	# Inferred admm step-size
	const ν_admm = 1/(ρ_admm+1)
	const nxnu = nx+nu

	# Temp storage variables
	const γ1 = randn()
	const γ2 = randn()
	const γ3 = randn()
	
	const κ1_zσ = MVector{2}(randn(2))
	const κ2_zσ = MVector{2}(randn(2))
	const κ1_x1 = MVector{nx1}(randn(nx1))
	const κ2_x1 = MVector{nx1}(randn(nx1))
	const κ1_x2 = MVector{nx2}(randn(nx2))
	const κ2_x2 = MVector{nx2}(randn(nx2))
	const κ1_u = MVector{nu}(randn(nu))
	const κ2_u = MVector{nu}(randn(nu))	
	const κ1_xu = MVector{nxnu}(randn(nxnu))
	const κ2_xu = MVector{nxnu}(randn(nxnu))	
	const κ3_xu = MVector{nxnu}(randn(nxnu))
	const κ4_xu = MVector{nxnu}(randn(nxnu))			# First element stores no. of projections performed in the last call to admm!		 			 

	const xu_u_admm = [MVector{nxnu}(randn(nxnu)) for _ in 1:N-1]
	const xu_y_admm = [MVector{nxnu}(randn(nxnu)) for _ in 1:N-1]

	function proj_xu1!(xu,zv,t)
		# Projection on segment 1 of x
		BLAS.blascopy!(nx1,zv,1,κ2_x1,1)
		proj.soc!(κ1_x1,κ2_x1,α_gs,nx1,γ1)
		BLAS.blascopy!(nx1,κ1_x1,1,xu,1)

		# Projection on segment 2 of x
		for j in ix2
			κ2_x2[j-nx1] = zv[j]
		end
		proj.ball!(κ1_x2,κ2_x2,r_vmax,nx2)
		for j in ix2
			xu[j] = κ1_x2[j-nx1]
		end

		# No projection on segment 1 of u 
		for j=nx+1:nx+nu1
			xu[j] = zv[j]
		end

		# Projection on segment 3 of x and segment 2 of u
		κ1_zσ[1] = zv[ix3]
		κ1_zσ[2] = zv[nx+iu2]
		proj.hlfspace2!(κ2_zσ,κ1_zσ,u1_zσ[t+1],u2_zσ[t+1],ζ1_zσ[t+1],ζ2_zσ[t+1],2,γ1,γ2,γ3)
		xu[ix3] = κ2_zσ[1]
		xu[nx+iu2] = κ2_zσ[2]  		
	end

	function proj_xu2!(xu,zv,t)
		# No projection on segments 1 and 2 of x
		BLAS.blascopy!(6,zv,1,xu,1)

		# Projection on segment 3 of x
		xu[ix3] = min(u_z[t+1],max(zv[ix3],l_z[t+1]))

		# Projection on u (involves both segments 1 and 2)
		for j in nx+1:nxnu
			κ2_u[j-nx] = zv[j]
		end
		proj.soc!(κ1_u,κ2_u,α_Tσ,nu,γ1)
		for j in nx+1:nxnu
			xu[j] = κ1_u[j-nx]
		end
	end	

	f_xu1_temp = Array{Function}(undef,N)
	f_xu2_temp = Array{Function}(undef,N)
	for t=1:N
		f_xu1_temp[t] = (x,z) -> proj_xu1!(x,z,t)
		f_xu2_temp[t] = (x,z) -> proj_xu2!(x,z,t)
	end
	const f_xu1 = tuple(f_xu1_temp...)
	const f_xu2 = tuple(f_xu2_temp...)

	function project_xu!(x,u,z,v,t)
		BLAS.blascopy!(nx,z,1,κ2_xu,1)
		for j=nx+1:nxnu
			κ2_xu[j] = v[j-nx]
		end
		proj.admm!(κ1_xu,xu_y_admm[t],xu_u_admm[t],κ2_xu,f_xu1[t],f_xu2[t],ρ_admm,ν_admm,ϵ_admm,nxnu,γ1,γ2,κ3_xu,κ4_xu)
		BLAS.blascopy!(nx,κ1_xu,1,x,1)
		for j=nx+1:nxnu
			u[j-nx] = κ1_xu[j]
		end		
	end

	# Diagnostics projection functions
	function project_xu_diagnostic!(x,u,z,v,t)
		project_xu!(x,u,z,v,t)
		proj_counter[1] = Int64(κ4_xu[1])
		proj_counter[2] = Int64(κ4_xu[1])
	end

	# Specify constraints in a JuMP problem
	function set_constr_JuMP!(model,x,u)			

		# Thrust upper and lower bound after lcvx
		@constraint(model,[t=1:N-1],u1_zσ[t+1]⋅cat(x[7,t+1],u[4,t],dims=1) ≤ ζ1_zσ[t+1])
		@constraint(model,[t=1:N-1],u2_zσ[t+1]⋅cat(x[7,t+1],u[4,t],dims=1) ≤ ζ2_zσ[t+1])

		@constraint(model,[t=1:N-1],[α_Tσ*u[4,t],u[1:3,t]...] in SecondOrderCone())

		# Mass upper and lower bound
		@constraint(model,[t=1:N-1],l_z[t+1] ≤ x[7,t+1] ≤ u_z[t+1])

		# Glide-slope constraint 
		@constraint(model,[t=1:N-1],[α_gs*x[3,t+1],x[1:2,t+1]...] in SecondOrderCone())

		# Velocity magnitude upper bound
		@constraint(model,[t=1:N-1],[r_vmax,x[4:6,t+1]...] in SecondOrderCone())

		# Thrust pointing constraint
		# @constraint(model,[t=1:N-1],[α_tl*u[3,t],u[1:2,t]...] in SecondOrderCone())
	end
end	
