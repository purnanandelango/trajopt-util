#  minimize   0.5 zᵀHz + hᵀz
# subject to  Gz = g 
#   		  z ∈ 𝔻 
module utils		
# This module contains: 
# 1) variables updated during power iteration
# 2) function for computing ground truth estimate of maximum singular value of matrix G from the linear equality constraint: Gz=g
# 3) function for computing the equality constraint error (Gz-g) and the relative distance between two solutions: |z1-z2|/|z2|, which is same as the relative 
#    distance to optimum (rd2o) when z2 is the optimal solution

# Module eX must be loaded before attempting to load this module
using LinearAlgebra, StaticArrays, JuMP	      				# Required packages
using ECOS, Gurobi, MosekTools, COSMO, SCS 			  	    # Optimization solvers
using ..eX	  							      				# Module with problem definition

	# Power iteration
		# Container variables
		const xL = [MVector{eX.nx}(randn(eX.nx)) for _ in 1:eX.N]	   # 0:N-1
		const uL = [MVector{eX.nu}(randn(eX.nu)) for _ in 1:(eX.N-1)]  # 0:N-2
		const vL = [MVector{eX.nx}(randn(eX.nx)) for _ in 1:eX.N]	   # 1:N
		# Special initialization for power iterations
		xL[1] .= zeros(eX.nx)
		vL[eX.N] .= zeros(eX.nx) 

		function reset_powiter!(xLL,uLL,vLL)
			# Reset the containers for power iterations
			for t = 1:eX.N-1
				xLL[t+1] .= randn(eX.nx)
				uLL[t] .= randn(eX.nu)
				vLL[t] .= randn(eX.nx)
			end
			# Special initialization for power
			xLL[1] .= zeros(eX.nx)
			vLL[eX.N] .= zeros(eX.nx)
		end


	# Estimate max singular value σ of G, matrix in the linear equality constraint of the stacked problem
	function compute_σ()
	    nunx = eX.nu+eX.nx   
	    
	    # Constuct G and H
	    G1 = zeros(Float64,eX.nx,nunx*(eX.N-1))
	    G2 = zeros(Float64,(eX.N-2)*eX.nx,nunx*(eX.N-1))

	    G1[1:eX.nx,1:nunx] = [-eX.Bd Diagonal(ones(eX.nx))]
	       
	    for j = 1:eX.N-2
	         G2[(j-1)*eX.nx+1:j*eX.nx,eX.nu+(j-1)*nunx+1:eX.nu+j*nunx+eX.nx] = hcat(-eX.Ad,-eX.Bd,Diagonal(ones(eX.nx)))
	    end
	    G = vcat(G1,G2);
	    
	    # Computation of σ, μ and λ
	    sGTG = svd(transpose(G)*G)
	    σ = sGTG.S[1]
	    
	    return σ
	end

	const σ = compute_σ() # Maximum singular value of G

	function set_scale!(x,u,flg::Symbol)
		if flg == :scale
			x[1] .= eX.scl_x_mat*x[1]
			for t = 1:eX.N-1
				x[t+1] .= eX.scl_x_mat*x[t+1]
				u[t] .= eX.scl_u_mat*u[t]
			end
		elseif flg == :unscale
			x[1] .= eX.scl_x_imat*x[1]
			for t = 1:eX.N-1
				x[t+1] .= eX.scl_x_imat*x[t+1]
				u[t] .= eX.scl_u_imat*u[t]
			end		
		else
			error("Invalid flag for scaling.")
		end 
	end

	function compute_error(x1,u1,x2,u2,print_flag=true,scl=:scaled,kmax=eX.kmax_pipg)
	# Reports the distance to optimum and dye for unscales solution variables
		θ = 0.0
		ϕ = 0.0 
		ψ = 0.0
		if eX.err_type == :linf   # ∞-norm
			if scl == :scaled
				ψ = max(ψ,norm(eX.scl_x_mat*x2[1],Inf))
			elseif scl == :unscaled 
				ψ = max(ψ,norm(x2[1],Inf))
			else
				error("Invalid scaling flag.")
			end
			for t = 1:eX.N-1
				if scl == :scaled
					θ = max(θ,norm(eX.scl_x_mat*(x1[t+1] .- x2[t+1]),Inf))
					θ = max(θ,norm(eX.scl_u_mat*(u1[t] .- u2[t]),Inf))
					ψ = max(ψ,norm(eX.scl_x_mat*x2[t+1],Inf))
					ψ = max(ψ,norm(eX.scl_u_mat*u2[t],Inf))
					ϕ = max(ϕ,norm(eX.scl_x_mat*(x1[t+1] .- eX.Ad*x1[t] .- eX.Bd*u1[t] .- eX.gd),Inf))
				elseif scl == :unscaled
					θ = max(θ,norm(x1[t+1] .- x2[t+1],Inf))
					θ = max(θ,norm(u1[t] .- u2[t],Inf))
					ψ = max(ψ,norm(x2[t+1],Inf))
					ψ = max(ψ,norm(u2[t],Inf))
					ϕ = max(ϕ,norm(x1[t+1] .- eX.Ad_unscl*x1[t] .- eX.Bd_unscl*u1[t] .- eX.gd_unscl,Inf))
				end
			end
		elseif eX.err_type == :l2 # 2-norm
			if scl == :scaled
				ψ += norm(eX.scl_x_mat*x2[1],2)^2
			elseif scl == :unscaled 
				ψ += norm(x2[1],2)^2
			else
				error("Invalid scaling flag.")
			end
			for t = 1:eX.N-1
				if scl == :scaled
					θ += norm(eX.scl_x_mat*(x1[t+1] .- x2[t+1]),2)^2
					θ += norm(eX.scl_u_mat*(u1[t] .- u2[t]),2)^2
					ψ += norm(eX.scl_x_mat*x2[t+1],2)^2
					ψ += norm(eX.scl_u_mat*u2[t],2)^2
					ϕ += norm(eX.scl_x_mat*(x1[t+1] .- eX.Ad*x1[t] .- eX.Bd*u1[t] .- eX.gd),2)^2
				elseif scl == :unscaled
					θ += norm(x1[t+1] .- x2[t+1],2)^2
					θ += norm(u1[t] .- u2[t],2)^2
					ψ += norm(x2[t+1],2)^2
					ψ += norm(u2[t],2)^2
					ϕ += norm(x1[t+1] .- eX.Ad_unscl*x1[t] .- eX.Bd_unscl*u1[t] .- eX.gd_unscl,2)^2
				end
			end
			θ = sqrt(θ)
			ψ = sqrt(ψ)
			ϕ = sqrt(ϕ)
		end	
		if print_flag
			print("Total PIPG iterations          : $kmax\n")
			print("Norm of optimum                : $ψ\n")
			print("Distance to optimum (absolute) : $θ\nDistance to optimum (relative) : $(θ/ψ) \nDynamics error                 : $ϕ")
		end
		return (θ/ψ),ϕ
	end

	# Solve problem via JuMP
	solver_JuMP_choice = [:gurobi]
	function solve_JuMP!(xopt,uopt,slvr=eX.solver_JuMP,verbosity=true,is_mpc=false,xinit=eX.x0,yref=eX.y,ϵ_pd=eX.ϵ_pd_JuMP,ϵ_gap=eX.ϵ_gap_JuMP)
		if slvr == :ecos
			model= Model(ECOS.Optimizer)
			set_optimizer_attribute(model,"printlevel",0)
			set_optimizer_attribute(model,"feastol",ϵ_pd)
			set_optimizer_attribute(model,"abstol",ϵ_gap)
			set_optimizer_attribute(model,"reltol",ϵ_gap)
		elseif slvr == :gurobi
			model= Model(Gurobi.Optimizer)
			set_optimizer_attribute(model,"Presolve",0)
			set_optimizer_attribute(model,"FeasibilityTol",ϵ_pd)
			set_optimizer_attribute(model,"OptimalityTol",ϵ_pd)
			set_optimizer_attribute(model,"BarConvTol",ϵ_gap)
			set_optimizer_attribute(model,"BarQCPConvTol",ϵ_gap)
		elseif slvr == :mosek
			model= Model(Mosek.Optimizer)
			set_optimizer_attribute(model,"INTPNT_CO_TOL_DFEAS",ϵ_pd)
			set_optimizer_attribute(model,"INTPNT_CO_TOL_PFEAS",ϵ_pd)
			set_optimizer_attribute(model,"INTPNT_CO_TOL_REL_GAP",ϵ_gap)
		elseif slvr == :cosmo
			model= Model(COSMO.Optimizer)
			set_optimizer_attribute(model,"eps_abs",ϵ_pd)
			set_optimizer_attribute(model,"eps_rel",ϵ_pd)
			set_optimizer_attribute(model,"eps_prim_inf",ϵ_pd)
			set_optimizer_attribute(model,"eps_dual_inf",ϵ_pd)
		elseif slvr == :scs
			model = Model(SCS.Optimizer)
			set_optimizer_attribute(model,"eps_abs",ϵ_pd)
			set_optimizer_attribute(model,"eps_rel",ϵ_pd)
			set_optimizer_attribute(model,"eps_infeas",ϵ_pd)
		elseif slvr == :osqp
			model = Model(OSQP.Optimizer)
			set_optimizer_attribute(model,"eps_abs",ϵ_pd)
			set_optimizer_attribute(model,"eps_rel",ϵ_pd)
			set_optimizer_attribute(model,"eps_prim_inf",ϵ_pd)
			set_optimizer_attribute(model,"eps_dual_inf",ϵ_pd)
		elseif slvr == :ipopt
			model = Model(Ipopt.Optimizer)
		else
			error("Invalid solver choice for JuMP.")
		end
		if ~verbosity
			set_silent(model)
		end

		@variable(model,x[1:eX.nx,1:eX.N])
		@variable(model,u[1:eX.nu,1:eX.N-1])
		@variable(model,β1[1:eX.N-1] ≥ 0)        			# Aux variable for defining objective
		@variable(model,β2[1:eX.N-1] ≥ 0)					# Aux variable for defining objective

		# Dynamics constraint
		Inx = Array(Diagonal(ones(eX.nx)))
		@constraint(model,[t=1:eX.N-1],[Inx -eX.Ad -eX.Bd -Inx]*cat(x[:,t+1],x[:,t],u[:,t],eX.gd,dims=1) .== 0)

		# Initial condition
		@constraint(model,x[:,1] .== xinit)

		# Set problem constraints
		eX.set_constr_JuMP!(model,x,u)

		# Define objective terms
		if eX.uref == nothing
			@constraint(model,[t=1:eX.N-1],[β1[t],(sqrt(eX.R)*u[:,t])...] in SecondOrderCone())
		else
			@constraint(model,[t=1:eX.N-1],[β1[t],(sqrt(eX.R)*(u[:,t].-eX.uref[t]))...] in SecondOrderCone())
		end
		@constraint(model,[t=1:eX.N-2],[β2[t],(sqrt(eX.Q)*(x[:,t+1] .- yref[t+1]))...] in SecondOrderCone())
		@constraint(model,[β2[eX.N-1],(sqrt(eX.Qf)*(x[:,eX.N] .- yref[eX.N]))...] in SecondOrderCone())

		@objective(model,Min,0.5*sum(β1 .^ 2) + 0.5*sum(β2 .^ 2))

		optimize!(model)

		exit_status = string(termination_status(model))
		@assert (exit_status ∈ ("OPTIMAL","SLOW_PROGRESS","ALMOST_OPTIMAL")) "Problem not solved correctly."

		compute_time = solve_time(model)

		for t=1:eX.N-1
			xopt[t] .= value.(x)[:,t]
			uopt[t] .= value.(u)[:,t]
		end
		xopt[eX.N] .= value.(x)[:,eX.N]

		solver_JuMP_choice[1] = slvr
		return exit_status,compute_time
	end

end

module common_plotter
using ..eX	  							      # Module with problem definition
using Plots, LaTeXStrings
pgfplotsx()
# This ensures that a legend entry is not created by default
default(lab="",markersize=2,markerstrokewidth=0.1,xtickfontsize=12,ytickfontsize=12,
	ztickfontsize=12,legendfontsize=9)

	# Plot relative distance to optimum and dynamics error for the unscaled solution
	function solution_quality(rd2o,dye,kmax=eX.kmax_pipg)

		p1 = plot(1:kmax-2,rd2o[1:kmax-2],line=:solid,color=:blue,lw=1.1,lab="",yaxis=:log,xlabel="Iteration",title="Distance to Optimum (relative)")
		p2 = plot(1:kmax-2,dye[1:kmax-2],line=:solid,color=:red,lw=1.1,lab="",yaxis=:log,xlabel="Iteration",title="Dynamics Violation")

		λ = @layout [a b]

		display(plot(p1,p2,layout=λ,size=[800,400]))

	end

	# Plot average number of projections per iteration
	function projection_count(px,pu,px_limit,pu_limit,kmax=eX.kmax_pipg)

		# y-axis limit
		ylimset = [0.001+0.5*min(minimum(px[1:kmax-2]),minimum(pu[1:kmax-2]),pu_limit,px_limit),1.5*max(10/15,maximum(px[1:kmax-2]),maximum(pu[1:kmax-2]))]

		p1 = plot(1:kmax-2,[px[k] .+ 0.001 for k=1:kmax-2],line=:solid,color=:blue,lw=:1.1,lab="",xlabel="Iteration",title="State Constraint Projections",yaxis=:log)
		plot!(ylim=ylimset)
		plot!(1:kmax-2,(px_limit+0.001) .* ones(kmax-2),line=:dash,color=:cyan,lw=1,lab="",yaxis=:log)
		plot!(ylabel="Projection/Horizon")
		p2 = plot(1:kmax-2,[pu[k] .+ 0.001 for k=1:kmax-2],line=:solid,color=:red,lw=:1.1,lab="",xlabel="Iteration",title="Input Constraint Projections",yaxis=:log)	
		plot!(ylim=ylimset)
		plot!(1:kmax-2,(pu_limit+0.001) .* ones(kmax-2),line=:dash,color=:cyan,lw=1,lab="Expected asymptote",legend=:bottomright,yaxis=:log)
		plot!(ylabel="Projection/Horizon")

		λ =@layout [a b]

		display(plot(p1,p2,layout=λ,size=[800,400]))

	end

end


module pipg 
# Provides implementation of the proportional-integral projected gradient (pipg) algorithm
# for solving trajectory optimization problems with strongly convex objective function for e.g. LQR, MPC sub-problems
# See https://arxiv.org/abs/2009.06980 for details
# ---
# Before loading pipg, module eX with problem definition, and module utils with utility functions must be loaded

using LinearAlgebra, StaticArrays	      # Required packages
using ..eX	  							  # Module with problem definition
using ..utils							  # Module with utility functions

	function powiter!(x,u,v,αk,βk,σ1,σ2,θx1,θx2)
		# Power iteration for estimating σ

		# Main iteration
		αk = abs(σ1-σ2)
		@fastmath @inbounds while αk ≥ eX.ϵ_powiter
			σ2 = copy(σ1)
			
			βk = 0.0
			@fastmath @inbounds @simd for t = 1:eX.N-1
				for j = 1:eX.nx
					σ1 = x[t+1][j]
					βk += σ1*σ1
				end
				for j = 1:eX.nu
					σ1 = u[t][j]
					βk += σ1*σ1
				end				
			end
			σ1 = sqrt(βk)

			βk = 1/σ1
			@fastmath @inbounds @simd for t=1:eX.N-1
				lmul!(βk,x[t+1])
				lmul!(βk,u[t])
			end

			@fastmath @inbounds @simd for t = 1:eX.N-1
				mul!(θx1,eX.Ad,x[t])
				mul!(θx2,eX.Bd,u[t])
				v[t] .= x[t+1] .- θx1 .- θx2
			end

			@fastmath @inbounds @simd for t = 1:eX.N-1
				mul!(θx1,eX.mAdT,v[t+1])
				x[t+1] .= v[t] .+ θx1
				mul!(u[t],eX.mBdT,v[t])
			end

			αk = abs(σ1-σ2)
		end

		return σ1
	end # powiter!

	# PIPG solver 
	function solver!(x,u,v,w,q,k,αk,βk,σ1,σ2,σ3,θx1,θx2,θu1,θu2,yref=eX.y,ϵp=eX.ϵ_primal,ϵd=eX.ϵ_dual)

		# σ1 = powiter!(utils.xL,utils.uL,utils.vL,αk,βk,σ1,σ2,θx1,θx2)

		σ1 = eX.μ/(2.2*σ1) 						# μ/(2σ) (make the σ estimate an over approximation)
		σ2 = eX.μ + 2*eX.λ 						# μ + 2λ

		@inbounds @simd for t = 1:eX.N-1	
			θx1 .= q[t]
			lmul!(2.0*σ1,θx1)
			v[t] .= w[t] .+ θx1
		end

		# PIPG iteration
		k = 1.0 # Iteration counter

		if eX.uref == nothing

				@fastmath @inbounds while k ≤ eX.kmax_pipg
					αk = 2/(k*eX.μ + σ2)				# PIPG step size

					σ3 = 0.0							# For primal stop crit.
					βk = 0.0 							# For primal stop crit.
					@fastmath @inbounds @simd for t = 1:eX.N-1
						mul!(θu1,eX.mBdT,v[t])
						mul!(θu2,eX.R,u[t])
						θu1 .= θu1 .+ θu2
						lmul!(-αk,θu1)
						θu2 .= u[t] .+ θu1
						θu1 .= u[t] 					# For primal stop crit.
						eX.project_u!(u[t],θu2,t)		# Projection onto input constraint set
						θu1 .= u[t] .- θu1				# For primal stop crit.

						θx1 .= x[t+1] .- yref[t+1]
						if t < eX.N-1
							mul!(θx2,eX.Q,θx1)
						else
							mul!(θx2,eX.Qf,θx1)
						end
						mul!(θx1,eX.mAdT,v[t+1])
						θx1 .= θx1 .+ v[t] .+ θx2
						lmul!(-αk,θx1)
						θx2 .= x[t+1] .+ θx1
						θx1 .= x[t+1] 					# For primal stop crit.
						eX.project_x!(x[t+1],θx2,t)		# Projection onto state constraint set
						θx1 .= x[t+1] .- θx1			# For primal stop crit.

						# Computation of primal stopping criteria
						# Relative change of primal variable: |zkp1-zk|/|zkp1| 
						for j = 1:eX.nx
							σ3 = max(σ3,abs(θx1[j]))
							βk = max(βk,abs(x[t+1][j]))
						end
						for j = 1:eX.nu
							σ3 = max(σ3,abs(θu1[j]))
							βk = max(βk,abs(u[t][j]))
						end
					end
					σ3 = σ3/βk							# Primal stop crit. value

					βk = 0.0							# For dual stop crit.
					αk = 0.0							# For dual stop crit.
					@fastmath @inbounds @simd for t = 1:eX.N-1
						mul!(θx1,eX.Ad,x[t])
						mul!(θx2,eX.Bd,u[t])		
						q[t] .= x[t+1] .- θx1 .- θx2 .- eX.gd
						axpy!((k+1.0)*σ1,q[t],w[t])		# PIPG step size: (k+1)μ/(2σ)
						θx1 .= q[t]
						lmul!((k+2.0)*σ1,θx1)			# PIPG step size: (k+2)μ/(2σ)
						v[t] .= w[t] .+ θx1	

						# Computation of dual stopping criteria
						# |G*zkp1-g|/|wkp1|
						for j=1:eX.nx
							βk = max(βk,abs(q[t][j]))
							αk = max(αk,abs(w[t][j]))
						end
					end

					# Primal-Dual stopping criteria
					if σ3 ≤ ϵp 									# Primal test
						if βk ≤ max(eX.ϵ_abs,ϵd*αk)				# Dual test

							# Current iterate is acceptable
							θx1[1] = k+1 # Store current iteration index
							break

						end
					end

					θx1[1] = k # Store current iteration index
					k += 1
				end			

		else

				@fastmath @inbounds while k ≤ eX.kmax_pipg
					αk = 2/(k*eX.μ + σ2)				# PIPG step size

					σ3 = 0.0							# For primal stop crit.
					βk = 0.0 							# For primal stop crit.
					@fastmath @inbounds @simd for t = 1:eX.N-1
						θu2 .= u[t] .- eX.uref[t]
						mul!(θu1,eX.R,θu2)
						mul!(θu2,eX.mBdT,v[t])
						θu1 .= θu1 .+ θu2
						lmul!(-αk,θu1)
						θu2 .= u[t] .+ θu1
						θu1 .= u[t] 					# For primal stop crit.
						eX.project_u!(u[t],θu2,t)		# Projection onto input constraint set
						θu1 .= u[t] .- θu1				# For primal stop crit.

						θx1 .= x[t+1] .- yref[t+1]
						if t < eX.N-1
							mul!(θx2,eX.Q,θx1)
						else
							mul!(θx2,eX.Qf,θx1)
						end
						mul!(θx1,eX.mAdT,v[t+1])
						θx1 .= θx1 .+ v[t] .+ θx2
						lmul!(-αk,θx1)
						θx2 .= x[t+1] .+ θx1
						θx1 .= x[t+1] 					# For primal stop crit.
						eX.project_x!(x[t+1],θx2,t)		# Projection onto state constraint set
						θx1 .= x[t+1] .- θx1			# For primal stop crit.

						# Computation of primal stopping criteria
						# Relative change of primal variable: |zkp1-zk|/|zkp1| 
						for j = 1:eX.nx
							σ3 = max(σ3,abs(θx1[j]))
							βk = max(βk,abs(x[t+1][j]))
						end
						for j = 1:eX.nu
							σ3 = max(σ3,abs(θu1[j]))
							βk = max(βk,abs(u[t][j]))
						end
					end
					σ3 = σ3/βk							# Primal stop crit. value

					βk = 0.0							# For dual stop crit.
					αk = 0.0							# For dual stop crit.
					@fastmath @inbounds @simd for t = 1:eX.N-1
						mul!(θx1,eX.Ad,x[t])
						mul!(θx2,eX.Bd,u[t])		
						q[t] .= x[t+1] .- θx1 .- θx2 .- eX.gd
						axpy!((k+1.0)*σ1,q[t],w[t])		# PIPG step size: (k+1)μ/(2σ)
						θx1 .= q[t]
						lmul!((k+2.0)*σ1,θx1)			# PIPG step size: (k+2)μ/(2σ)
						v[t] .= w[t] .+ θx1	

						# Computation of dual stopping criteria
						# |G*zkp1-g|/|wkp1|
						for j=1:eX.nx
							βk = max(βk,abs(q[t][j]))
							αk = max(αk,abs(w[t][j]))
						end
					end

					# Primal dual stopping criteria
					if σ3 ≤ eX.ϵp							# Primal test
						if βk ≤ max(eX.ϵ_abs,ϵd*αk)			# Dual test

							# Current iterate is acceptable
							θx1[1] = k+1 # Store current iteration index
							break

						end
					end

					θx1[1] = k # Store current iteration index
					k += 1
				end

		end

	end # solver!

	# Same as solver! but records diagnostic information
	function solver_diagnostic!(x,u,v,w,q,k,αk,βk,σ1,σ2,σ3,θx1,θx2,θu1,θu2,xopt,uopt,rd2o,dye,proj_count_x,proj_count_u,term_crit=true,yref=eX.y,powiter_text=true,ϵp=eX.ϵ_primal,ϵd=eX.ϵ_dual)

		# σ1 = powiter!(utils.xL,utils.uL,utils.vL,αk,βk,σ1,σ2,θx1,θx2)

		if powiter_text
			print("Power iteration estimate of σ  : $σ1\nSVD estimate of σ              : $(utils.σ)")
		end
			
		σ1 = eX.μ/(2.2*σ1) 						# μ/(2σ) (make the σ estimate an over approximation)
		σ2 = eX.μ + 2*eX.λ 						# μ + 2λ

		for t = 1:eX.N-1	
			θx1 .= q[t]
			lmul!(2.0*σ1,θx1)
			v[t] .= w[t] .+ θx1
		end

		# PIPG iteration
		k = 1.0 # Iteration counter
		if eX.uref == nothing

				while k ≤ eX.kmax_pipg
					αk = 2/(k*eX.μ + σ2)				# PIPG step size

					σ3 = 0.0							# For primal stop crit.
					βk = 0.0 							# For primal stop crit.
					for t = 1:eX.N-1
						mul!(θu1,eX.mBdT,v[t])
						mul!(θu2,eX.R,u[t])
						θu1 .= θu1 .+ θu2
						lmul!(-αk,θu1)
						θu2 .= u[t] .+ θu1
						θu1 .= u[t] 								# For primal stop crit.
						eX.project_u_diagnostic!(u[t],θu2,t)		# Projection onto input constraint set
						θu1 .= u[t] .- θu1							# For primal stop crit.

						# Count no. of projections on input constraints set
						proj_count_u[Int64(k)] += eX.proj_counter[2]

						θx1 .= x[t+1] .- yref[t+1]
						if t < eX.N-1
							mul!(θx2,eX.Q,θx1)
						else
							mul!(θx2,eX.Qf,θx1)
						end
						mul!(θx1,eX.mAdT,v[t+1])
						θx1 .= θx1 .+ v[t] .+ θx2
						lmul!(-αk,θx1)
						θx2 .= x[t+1] .+ θx1
						θx1 .= x[t+1] 								# For primal stop crit.
						eX.project_x_diagnostic!(x[t+1],θx2,t)		# Projection onto state constraint set
						θx1 .= x[t+1] .- θx1						# For primal stop crit.

						# Count no. of projections on input constraints set
						proj_count_x[Int64(k)] += eX.proj_counter[1]

						# Computation of primal stopping criteria
						# Relative change of primal variable: |zkp1-zk|/|zkp1| 
						for j = 1:eX.nx
							σ3 = max(σ3,abs(θx1[j]))
							βk = max(βk,abs(x[t+1][j]))
						end
						for j = 1:eX.nu
							σ3 = max(σ3,abs(θu1[j]))
							βk = max(βk,abs(u[t][j]))
						end
					end
					σ3 = σ3/βk							# Primal stop crit. value

					# Normalize projection count by horizon length
					proj_count_u[Int64(k)] = proj_count_u[Int64(k)]/(eX.N-1)
					proj_count_x[Int64(k)] = proj_count_x[Int64(k)]/(eX.N-1)

					βk = 0.0							# For dual stop crit.
					αk = 0.0							# For dual stop crit.
					for t = 1:eX.N-1
						mul!(θx1,eX.Ad,x[t])
						mul!(θx2,eX.Bd,u[t])		
						q[t] .= x[t+1] .- θx1 .- θx2 .- eX.gd
						axpy!((k+1.0)*σ1,q[t],w[t])		# PIPG step size: (k+1)μ/(2σ)
						θx1 .= q[t]
						lmul!((k+2.0)*σ1,θx1)			# PIPG step size: (k+2)μ/(2σ)
						v[t] .= w[t] .+ θx1	

						# Computation of dual stopping criteria
						# |G*zkp1-g|/|wkp1|
						for j=1:eX.nx
							βk = max(βk,abs(q[t][j]))
							αk = max(αk,abs(w[t][j]))
						end
					end

					if term_crit
						# Primal dual stopping criteria
						if σ3 ≤ ϵp								# Primal test
							if βk ≤ max(eX.ϵ_abs,ϵd*αk)			# Dual test

								# Current iterate is acceptable
								θx1[1] = k+1 # Store current iteration index
								# print("\nTotal PIPG iterations        : $(k+1)")
								break

							end
						end
					end

					rd2o[Int64(k)],dye[Int64(k)] = utils.compute_error(x,u,xopt,uopt,false)

					k += 1

					## Diagnostics
					# print("\nIter:\n")
					# print(Array(x))
					# print("\n")
					# print(Array(u))
					# print("\n")
				end

		else

				while k ≤ eX.kmax_pipg
					αk = 2/(k*eX.μ + σ2)				# PIPG step size

					σ3 = 0.0							# For primal stop crit.
					βk = 0.0 							# For primal stop crit.
					for t = 1:eX.N-1
						θu2 .= u[t] .- eX.uref[t]
						mul!(θu1,eX.R,θu2)
						mul!(θu2,eX.mBdT,v[t])
						θu1 .= θu1 .+ θu2
						lmul!(-αk,θu1)
						θu2 .= u[t] .+ θu1
						θu1 .= u[t] 								# For primal stop crit.
						eX.project_u_diagnostic!(u[t],θu2,t)		# Projection onto input constraint set
						θu1 .= u[t] .- θu1							# For primal stop crit.

						# Count no. of projections on input constraints set
						proj_count_u[Int64(k)] += eX.proj_counter[2]

						θx1 .= x[t+1] .- yref[t+1]
						if t < eX.N-1
							mul!(θx2,eX.Q,θx1)
						else
							mul!(θx2,eX.Qf,θx1)
						end
						mul!(θx1,eX.mAdT,v[t+1])
						θx1 .= θx1 .+ v[t] .+ θx2
						lmul!(-αk,θx1)
						θx2 .= x[t+1] .+ θx1
						θx1 .= x[t+1] 								# For primal stop crit.
						eX.project_x_diagnostic!(x[t+1],θx2,t)		# Projection onto state constraint set
						θx1 .= x[t+1] .- θx1						# For primal stop crit.

						# Count no. of projections on input constraints set
						proj_count_x[Int64(k)] += eX.proj_counter[1]

						# Computation of primal stopping criteria
						# Relative change of primal variable: |zkp1-zk|/|zkp1| 
						for j = 1:eX.nx
							σ3 = max(σ3,abs(θx1[j]))
							βk = max(βk,abs(x[t+1][j]))
						end
						for j = 1:eX.nu
							σ3 = max(σ3,abs(θu1[j]))
							βk = max(βk,abs(u[t][j]))
						end
					end
					σ3 = σ3/βk							# Primal stop crit. value

					# Normalize projection count by horizon length
					proj_count_u[Int64(k)] = proj_count_u[Int64(k)]/(eX.N-1)
					proj_count_x[Int64(k)] = proj_count_x[Int64(k)]/(eX.N-1)

					βk = 0.0							# For dual stop crit.
					αk = 0.0							# For dual stop crit.
					for t = 1:eX.N-1
						mul!(θx1,eX.Ad,x[t])
						mul!(θx2,eX.Bd,u[t])		
						q[t] .= x[t+1] .- θx1 .- θx2 .- eX.gd
						axpy!((k+1.0)*σ1,q[t],w[t])		# PIPG step size: (k+1)μ/(2σ)
						θx1 .= q[t]
						lmul!((k+2.0)*σ1,θx1)			# PIPG step size: (k+2)μ/(2σ)
						v[t] .= w[t] .+ θx1	

						# Computation of dual stopping criteria
						# |G*zkp1-g|/|wkp1|
						for j=1:eX.nx
							βk = max(βk,abs(q[t][j]))
							αk = max(αk,abs(w[t][j]))
						end
					end

					if term_crit
						# Primal dual stopping criteria
						if σ3 ≤ ϵp									# Primal test
							if βk ≤ max(eX.ϵ_abs,ϵd*αk)				# Dual test

								# Current iterate is acceptable
								θx1[1] = k+1 # Store current iteration index
								# print("\nTotal PIPG iterations        : $(k+1)")
								break

							end
						end
					end

					rd2o[Int64(k)],dye[Int64(k)] = utils.compute_error(x,u,xopt,uopt,false)

					k += 1

					## Diagnostics
					# print("\nIter:\n")
					# print(Array(x))
					# print("\n")
					# print(Array(u))
					# print("\n")
				end

		end
		
		if k > eX.kmax_pipg-1 
			θx1[1] = k-1 # = eX.kmax_pipg
			# print("\nTotal PIPG iterations        : $(eX.kmax_pipg)")
		end
	end # solver_diagnostic!

	# Same as solver! but explicitly passes the variables updated in the calls to admm 
	function solver_v2!(x,u,v,w,q,k,αk,βk,σ1,σ2,σ3,θx1,θx2,θu1,θu2,yref,x_y_admm,x_u_admm,u_y_admm,u_u_admm,ϵp=eX.ϵ_primal,ϵd=eX.ϵ_dual)

		# σ1 = powiter!(utils.xL,utils.uL,utils.vL,αk,βk,σ1,σ2,θx1,θx2)

		σ1 = eX.μ/(2.2*σ1) 						# μ/(2σ) (make the σ estimate an over approximation)
		σ2 = eX.μ + 2*eX.λ 						# μ + 2λ

		@inbounds @simd for t = 1:eX.N-1	
			θx1 .= q[t]
			lmul!(2.0*σ1,θx1)
			v[t] .= w[t] .+ θx1
		end

		# PIPG iteration
		k = 1.0 # Iteration counter
		if eX.uref == nothing

				@fastmath @inbounds while k ≤ eX.kmax_pipg
					αk = 2/(k*eX.μ + σ2)				# Pipg step size

					σ3 = 0.0							# For primal stop crit.
					βk = 0.0 							# For primal stop crit.
					@fastmath @inbounds @simd for t = 1:eX.N-1
						mul!(θu1,eX.mBdT,v[t])
						mul!(θu2,eX.R,u[t])
						θu1 .= θu1 .+ θu2
						lmul!(-αk,θu1)
						θu2 .= u[t] .+ θu1
						θu1 .= u[t] 												# For primal stop crit.
						eX.project_u_v2!(u[t],θu2,t,u_y_admm[t],u_u_admm[t])		# Projection onto input constraint set
						θu1 .= u[t] .- θu1											# For primal stop crit.

						θx1 .= x[t+1] .- yref[t+1]
						if t < eX.N-1
							mul!(θx2,eX.Q,θx1)
						else
							mul!(θx2,eX.Qf,θx1)
						end
						mul!(θx1,eX.mAdT,v[t+1])
						θx1 .= θx1 .+ v[t] .+ θx2
						lmul!(-αk,θx1)
						θx2 .= x[t+1] .+ θx1
						θx1 .= x[t+1] 												# For primal stop crit.
						eX.project_x_v2!(x[t+1],θx2,t,x_y_admm[t],x_u_admm[t])		# Projection onto state constraint set
						θx1 .= x[t+1] .- θx1										# For primal stop crit.

						# Computation of primal stopping criteria
						# Relative change of primal variable: |zkp1-zk|/|zkp1| 
						for j = 1:eX.nx
							σ3 = max(σ3,abs(θx1[j]))
							βk = max(βk,abs(x[t+1][j]))
						end
						for j = 1:eX.nu
							σ3 = max(σ3,abs(θu1[j]))
							βk = max(βk,abs(u[t][j]))
						end
					end
					σ3 = σ3/βk							# Primal stop crit. value

					βk = 0.0							# For dual stop crit.
					αk = 0.0							# For dual stop crit.
					@fastmath @inbounds @simd for t = 1:eX.N-1
						mul!(θx1,eX.Ad,x[t])
						mul!(θx2,eX.Bd,u[t])		
						q[t] .= x[t+1] .- θx1 .- θx2 .- eX.gd
						axpy!((k+1.0)*σ1,q[t],w[t])		# PIPG step size: (k+1)μ/(2σ)
						θx1 .= q[t]
						lmul!((k+2.0)*σ1,θx1)			# PIPG step size: (k+2)μ/(2σ)
						v[t] .= w[t] .+ θx1	

						# Computation of dual stopping criteria
						# |G*zkp1-g|/|wkp1|
						for j=1:eX.nx
							βk = max(βk,abs(q[t][j]))
							αk = max(αk,abs(w[t][j]))
						end
					end

					# Primal dual stopping criteria
					if σ3 ≤ ϵp									# Primal test
						if βk ≤ max(eX.ϵ_abs,ϵd*αk)				# Dual test

							# Current iterate is acceptable
							θx1[1] = k+1 # Store current iteration index
							break

						end
					end

					θx1[1] = k # Store current iteration index
					k += 1
				end

		else

				@fastmath @inbounds while k ≤ eX.kmax_pipg
					αk = 2/(k*eX.μ + σ2)				# PIPG step size

					σ3 = 0.0							# For primal stop crit.
					βk = 0.0 							# For primal stop crit.
					@fastmath @inbounds @simd for t = 1:eX.N-1
						θu2 .= u[t] .- eX.uref[t]
						mul!(θu1,eX.R,θu2)
						mul!(θu2,eX.mBdT,v[t])
						θu1 .= θu1 .+ θu2
						lmul!(-αk,θu1)
						θu2 .= u[t] .+ θu1
						θu1 .= u[t] 												# For primal stop crit.
						eX.project_u_v2!(u[t],θu2,t,u_y_admm[t],u_u_admm[t])		# Projection onto input constraint set
						θu1 .= u[t] .- θu1											# For primal stop crit.

						θx1 .= x[t+1] .- yref[t+1]
						if t < eX.N-1
							mul!(θx2,eX.Q,θx1)
						else
							mul!(θx2,eX.Qf,θx1)
						end
						mul!(θx1,eX.mAdT,v[t+1])
						θx1 .= θx1 .+ v[t] .+ θx2
						lmul!(-αk,θx1)
						θx2 .= x[t+1] .+ θx1
						θx1 .= x[t+1] 												# For primal stop crit.
						eX.project_x_v2!(x[t+1],θx2,t,x_y_admm[t],x_u_admm[t])		# Projection onto state constraint set
						θx1 .= x[t+1] .- θx1										# For primal stop crit.

						# Computation of primal stopping criteria
						# Relative change of primal variable: |zkp1-zk|/|zkp1| 
						for j = 1:eX.nx
							σ3 = max(σ3,abs(θx1[j]))
							βk = max(βk,abs(x[t+1][j]))
						end
						for j = 1:eX.nu
							σ3 = max(σ3,abs(θu1[j]))
							βk = max(βk,abs(u[t][j]))
						end
					end
					σ3 = σ3/βk							# Primal stop crit. value

					βk = 0.0							# For dual stop crit.
					αk = 0.0							# For dual stop crit.
					@fastmath @inbounds @simd for t = 1:eX.N-1
						mul!(θx1,eX.Ad,x[t])
						mul!(θx2,eX.Bd,u[t])		
						q[t] .= x[t+1] .- θx1 .- θx2 .- eX.gd
						axpy!((k+1.0)*σ1,q[t],w[t])		# PIPG step size: (k+1)μ/(2σ)
						θx1 .= q[t]
						lmul!((k+2.0)*σ1,θx1)			# PIPG step size: (k+2)μ/(2σ)
						v[t] .= w[t] .+ θx1	

						# Computation of dual stopping criteria
						# |G*zkp1-g|/|wkp1|
						for j=1:eX.nx
							βk = max(βk,abs(q[t][j]))
							αk = max(αk,abs(w[t][j]))
						end
					end

					# Primal dual stopping criteria
					if σ3 ≤ ϵp 	 								# primal test
						if βk ≤ max(eX.ϵ_abs,ϵd*αk)				# dual test

							# Current iterate is acceptable
							θx1[1] = k+1 # Store current iteration index
							break

						end
					end

					θx1[1] = k # Store current iteration index
					k += 1
				end

		end
	end # solver!

	# Same as solver_v2! but records diagnostic information
	function solver_diagnostic_v2!(x,u,v,w,q,k,αk,βk,σ1,σ2,σ3,θx1,θx2,θu1,θu2,xopt,uopt,rd2o,dye,proj_count_x,proj_count_u,yref,x_y_admm,x_u_admm,u_y_admm,u_u_admm,term_crit=true,powiter_text=true,ϵp=eX.ϵ_primal,ϵd=eX.ϵ_dual)

		# σ1 = powiter!(utils.xL,utils.uL,utils.vL,αk,βk,σ1,σ2,θx1,θx2)

		if powiter_text
			print("Power iteration estimate of σ  : $σ1\nSVD estimate of σ              : $(utils.σ)")
		end
			
		σ1 = eX.μ/(2.2*σ1) 						# μ/(2σ) (make the σ estimate an over approximation)
		σ2 = eX.μ + 2*eX.λ 						# μ + 2λ

		for t = 1:eX.N-1	
			θx1 .= q[t]
			lmul!(2.0*σ1,θx1)
			v[t] .= w[t] .+ θx1
		end

		# PIPG iteration
		k = 1.0 # Iteration counter

		if eX.uref == nothing

				while k ≤ eX.kmax_pipg
					αk = 2/(k*eX.μ + σ2)				# PIPG step size

					σ3 = 0.0							# For primal stop crit.
					βk = 0.0 							# For primal stop crit.
					for t = 1:eX.N-1
						mul!(θu1,eX.mBdT,v[t])
						mul!(θu2,eX.R,u[t])
						θu1 .= θu1 .+ θu2
						lmul!(-αk,θu1)
						θu2 .= u[t] .+ θu1
						θu1 .= u[t] 														# For primal stop crit.
						eX.project_u_diagnostic_v2!(u[t],θu2,t,u_y_admm[t],u_u_admm[t])		# Projection onto input constraint set
						θu1 .= u[t] .- θu1													# For primal stop crit.

						# Count no. of projections on input constraints set
						proj_count_u[Int64(k)] += eX.proj_counter[2]

						θx1 .= x[t+1] .- yref[t+1]
						if t < eX.N-1
							mul!(θx2,eX.Q,θx1)
						else
							mul!(θx2,eX.Qf,θx1)
						end
						mul!(θx1,eX.mAdT,v[t+1])
						θx1 .= θx1 .+ v[t] .+ θx2
						lmul!(-αk,θx1)
						θx2 .= x[t+1] .+ θx1
						θx1 .= x[t+1] 														# For primal stop crit.
						eX.project_x_diagnostic_v2!(x[t+1],θx2,t,x_y_admm[t],x_u_admm[t])	# Projection onto state constraint set
						θx1 .= x[t+1] .- θx1												# For primal stop crit.

						# Count no. of projections on input constraints set
						proj_count_x[Int64(k)] += eX.proj_counter[1]

						# Computation of primal stopping criteria
						# Relative change of primal variable: |zkp1-zk|/|zkp1| 
						for j = 1:eX.nx
							σ3 = max(σ3,abs(θx1[j]))
							βk = max(βk,abs(x[t+1][j]))
						end
						for j = 1:eX.nu
							σ3 = max(σ3,abs(θu1[j]))
							βk = max(βk,abs(u[t][j]))
						end
					end
					σ3 = σ3/βk							# Primal stop crit. value

					# Normalize projection count by horizon length
					proj_count_u[Int64(k)] = proj_count_u[Int64(k)]/(eX.N-1)
					proj_count_x[Int64(k)] = proj_count_x[Int64(k)]/(eX.N-1)

					βk = 0.0							# For dual stop crit.
					αk = 0.0							# For dual stop crit.
					for t = 1:eX.N-1
						mul!(θx1,eX.Ad,x[t])
						mul!(θx2,eX.Bd,u[t])		
						q[t] .= x[t+1] .- θx1 .- θx2 .- eX.gd
						axpy!((k+1.0)*σ1,q[t],w[t])		# PIPG step size: (k+1)μ/(2σ)
						θx1 .= q[t]
						lmul!((k+2.0)*σ1,θx1)			# PIPG step size: (k+2)μ/(2σ)
						v[t] .= w[t] .+ θx1	

						# Computation of dual stopping criteria
						# |G*zkp1-g|/|wkp1|
						for j=1:eX.nx
							βk = max(βk,abs(q[t][j]))
							αk = max(αk,abs(w[t][j]))
						end
					end

					if term_crit
						# Primal dual stopping criteria
						if σ3 ≤ ϵp 									# Primal test
							if βk ≤ max(eX.ϵ_abs,ϵd*αk)				# Dual test

								# Current iterate is acceptable
								θx1[1] = k+1 # Store current iteration index
								# print("\nTotal PIPG iterations        : $(k+1)")
								break

							end
						end
					end

					rd2o[Int64(k)],dye[Int64(k)] = utils.compute_error(x,u,xopt,uopt,false)

					k += 1

					## Diagnostics
					# print("\nIter:\n")
					# print(Array(x))
					# print("\n")
					# print(Array(u))
					# print("\n")
				end

		else

				while k ≤ eX.kmax_pipg
					αk = 2/(k*eX.μ + σ2)				# PIPG step size

					σ3 = 0.0							# For primal stop crit.
					βk = 0.0 							# For primal stop crit.
					for t = 1:eX.N-1
						θu2 .= u[t] .- eX.uref[t]
						mul!(θu1,eX.R,θu2)
						mul!(θu2,eX.mBdT,v[t])
						θu1 .= θu1 .+ θu2
						lmul!(-αk,θu1)
						θu2 .= u[t] .+ θu1
						θu1 .= u[t] 														# For primal stop crit.
						eX.project_u_diagnostic_v2!(u[t],θu2,t,u_y_admm[t],u_u_admm[t])		# Projection onto input constraint set
						θu1 .= u[t] .- θu1													# For primal stop crit.

						# Count no. of projections on input constraints set
						proj_count_u[Int64(k)] += eX.proj_counter[2]

						θx1 .= x[t+1] .- yref[t+1]
						if t < eX.N-1
							mul!(θx2,eX.Q,θx1)
						else
							mul!(θx2,eX.Qf,θx1)
						end
						mul!(θx1,eX.mAdT,v[t+1])
						θx1 .= θx1 .+ v[t] .+ θx2
						lmul!(-αk,θx1)
						θx2 .= x[t+1] .+ θx1
						θx1 .= x[t+1] 														# For primal stop crit.
						eX.project_x_diagnostic_v2!(x[t+1],θx2,t,x_y_admm[t],x_u_admm[t])	# Projection onto state constraint set
						θx1 .= x[t+1] .- θx1												# For primal stop crit.

						# Count no. of projections on input constraints set
						proj_count_x[Int64(k)] += eX.proj_counter[1]

						# Computation of primal stopping criteria
						# Relative change of primal variable: |zkp1-zk|/|zkp1| 
						for j = 1:eX.nx
							σ3 = max(σ3,abs(θx1[j]))
							βk = max(βk,abs(x[t+1][j]))
						end
						for j = 1:eX.nu
							σ3 = max(σ3,abs(θu1[j]))
							βk = max(βk,abs(u[t][j]))
						end
					end
					σ3 = σ3/βk							# Primal stop crit. value

					# Normalize projection count by horizon length
					proj_count_u[Int64(k)] = proj_count_u[Int64(k)]/(eX.N-1)
					proj_count_x[Int64(k)] = proj_count_x[Int64(k)]/(eX.N-1)

					βk = 0.0							# For dual stop crit.
					αk = 0.0							# For dual stop crit.
					for t = 1:eX.N-1
						mul!(θx1,eX.Ad,x[t])
						mul!(θx2,eX.Bd,u[t])		
						q[t] .= x[t+1] .- θx1 .- θx2 .- eX.gd
						axpy!((k+1.0)*σ1,q[t],w[t])		# PIPG step size: (k+1)μ/(2σ)
						θx1 .= q[t]
						lmul!((k+2.0)*σ1,θx1)			# PIPG step size: (k+2)μ/(2σ)
						v[t] .= w[t] .+ θx1	

						# Computation of dual stopping criteria
						# |G*zkp1-g|/|wkp1|
						for j=1:eX.nx
							βk = max(βk,abs(q[t][j]))
							αk = max(αk,abs(w[t][j]))
						end
					end

					if term_crit
						# Primal dual stopping criteria
						if σ3 ≤ ϵp									# Primal test
							if βk ≤ max(eX.ϵ_abs,ϵd*αk)	# Dual test

								# Current iterate is acceptable
								θx1[1] = k+1 # Store current iteration index
								# print("\nTotal PIPG iterations        : $(k+1)")
								break

							end
						end
					end

					rd2o[Int64(k)],dye[Int64(k)] = utils.compute_error(x,u,xopt,uopt,false)

					k += 1

					## Diagnostics
					# print("\nIter:\n")
					# print(Array(x))
					# print("\n")
					# print(Array(u))
					# print("\n")
				end

		end

		if k > eX.kmax_pipg-1 
			θx1[1] = k-1 # = eX.kmax_pipg
			# print("\nTotal PIPG iterations        : $(eX.kmax_pipg)")
		end
	end # solver_diagnostic!

	# Handles coupled state and input constraint
	function solver_v3!(x,u,v,w,q,k,αk,βk,σ1,σ2,σ3,θx1,θx2,θu1,θu2,yref=eX.y,ϵp=eX.ϵ_primal,ϵd=eX.ϵ_dual)

		# σ1 = powiter!(utils.xL,utils.uL,utils.vL,αk,βk,σ1,σ2,θx1,θx2)

		σ1 = eX.μ/(2.2*σ1) 						# μ/(2σ) (make the σ estimate an over approximation)
		σ2 = eX.μ + 2*eX.λ 						# μ + 2λ

		@inbounds @simd for t = 1:eX.N-1	
			θx1 .= q[t]
			lmul!(2.0*σ1,θx1)
			v[t] .= w[t] .+ θx1
		end

		# PIPG iteration
		k = 1.0 # iteration counter

		if eX.uref == nothing

				@fastmath @inbounds while k ≤ eX.kmax_pipg
					αk = 2/(k*eX.μ + σ2)				# PIPG step size

					σ3 = 0.0							# For primal stop crit.
					βk = 0.0 							# For primal stop crit.
					@fastmath @inbounds @simd for t = 1:eX.N-1
						mul!(θu1,eX.mBdT,v[t])
						mul!(θu2,eX.R,u[t])
						θu1 .= θu1 .+ θu2
						lmul!(-αk,θu1)
						θu2 .= u[t] .+ θu1
						θu1 .= u[t] 					# For primal stop crit.
						# eX.project_u!(u[t],θu2,t)		# Projection onto input constraint set

						θx1 .= x[t+1] .- yref[t+1]
						if t < eX.N-1
							mul!(θx2,eX.Q,θx1)
						else
							mul!(θx2,eX.Qf,θx1)
						end
						mul!(θx1,eX.mAdT,v[t+1])
						θx1 .= θx1 .+ v[t] .+ θx2
						lmul!(-αk,θx1)
						θx2 .= x[t+1] .+ θx1
						θx1 .= x[t+1] 					# For primal stop crit.
						# eX.project_x!(x[t+1],θx2,t)	# projection onto state constraint set

						eX.project_xu!(x[t+1],u[t],θx2,θu2,t)

						θu1 .= u[t] .- θu1				# For primal stop crit.
						θx1 .= x[t+1] .- θx1			# For primal stop crit.

						# Computation of primal stopping criteria
						# Relative change of primal variable: |zkp1-zk|/|zkp1| 
						for j = 1:eX.nx
							σ3 = max(σ3,abs(θx1[j]))
							βk = max(βk,abs(x[t+1][j]))
						end
						for j = 1:eX.nu
							σ3 = max(σ3,abs(θu1[j]))
							βk = max(βk,abs(u[t][j]))
						end
					end
					σ3 = σ3/βk							# Primal stop crit. value

					βk = 0.0							# For dual stop crit.
					αk = 0.0							# For dual stop crit.
					@fastmath @inbounds @simd for t = 1:eX.N-1
						mul!(θx1,eX.Ad,x[t])
						mul!(θx2,eX.Bd,u[t])		
						q[t] .= x[t+1] .- θx1 .- θx2 .- eX.gd
						axpy!((k+1.0)*σ1,q[t],w[t])		# PIPG step size: (k+1)μ/(2σ)
						θx1 .= q[t]
						lmul!((k+2.0)*σ1,θx1)			# PIPG step size: (k+2)μ/(2σ)
						v[t] .= w[t] .+ θx1	

						# Computation of dual stopping criteria
						# |G*zkp1-g|/|wkp1|
						for j=1:eX.nx
							βk = max(βk,abs(q[t][j]))
							αk = max(αk,abs(w[t][j]))
						end
					end

					# Primal dual stopping criteria
					if σ3 ≤ ϵp 									# Primal test
						if βk ≤ max(eX.ϵ_abs,ϵd*αk)				# Dual test

							# Current iterate is acceptable
							θx1[1] = k+1 # Store current iteration index
							break

						end
					end

					θx1[1] = k # Store current iteration index
					k += 1
				end			

		else

				@fastmath @inbounds while k ≤ eX.kmax_pipg
					αk = 2/(k*eX.μ + σ2)				# PIPG step size

					σ3 = 0.0							# For primal stop crit.
					βk = 0.0 							# For primal stop crit.
					@fastmath @inbounds @simd for t = 1:eX.N-1
						θu2 .= u[t] .- eX.uref[t]
						mul!(θu1,eX.R,θu2)
						mul!(θu2,eX.mBdT,v[t])
						θu1 .= θu1 .+ θu2
						lmul!(-αk,θu1)
						θu2 .= u[t] .+ θu1
						θu1 .= u[t] 					# For primal stop crit.
						# eX.project_u!(u[t],θu2,t)		# Projection onto input constraint set
						
						θx1 .= x[t+1] .- yref[t+1]
						if t < eX.N-1
							mul!(θx2,eX.Q,θx1)
						else
							mul!(θx2,eX.Qf,θx1)
						end
						mul!(θx1,eX.mAdT,v[t+1])
						θx1 .= θx1 .+ v[t] .+ θx2
						lmul!(-αk,θx1)
						θx2 .= x[t+1] .+ θx1
						θx1 .= x[t+1] 					# For primal stop crit.
						# eX.project_x!(x[t+1],θx2,t)	# Projection onto state constraint set

						eX.project_xu!(x[t+1],u[t],θx2,θu2,t)

						θx1 .= x[t+1] .- θx1			# For primal stop crit.
						θu1 .= u[t] .- θu1				# For primal stop crit.

						# Computation of primal stopping criteria
						# Relative change of primal variable: |zkp1-zk|/|zkp1| 
						for j = 1:eX.nx
							σ3 = max(σ3,abs(θx1[j]))
							βk = max(βk,abs(x[t+1][j]))
						end
						for j = 1:eX.nu
							σ3 = max(σ3,abs(θu1[j]))
							βk = max(βk,abs(u[t][j]))
						end
					end
					σ3 = σ3/βk							# Primal stop crit. value

					βk = 0.0							# For dual stop crit.
					αk = 0.0							# For dual stop crit.
					@fastmath @inbounds @simd for t = 1:eX.N-1
						mul!(θx1,eX.Ad,x[t])
						mul!(θx2,eX.Bd,u[t])		
						q[t] .= x[t+1] .- θx1 .- θx2 .- eX.gd
						axpy!((k+1.0)*σ1,q[t],w[t])		# PIPG step size: (k+1)μ/(2σ)
						θx1 .= q[t]
						lmul!((k+2.0)*σ1,θx1)			# PIPG step size: (k+2)μ/(2σ)
						v[t] .= w[t] .+ θx1	

						# Computation of dual stopping criteria
						# |G*zkp1-g|/|wkp1|
						for j=1:eX.nx
							βk = max(βk,abs(q[t][j]))
							αk = max(αk,abs(w[t][j]))
						end
					end

					# Primal dual stopping criteria
					if σ3 ≤ ϵp  								# Primal test
						if βk ≤ max(eX.ϵ_abs,ϵd*αk)				# Dual test

							# Current iterate is acceptable
							θx1[1] = k+1 # store current iteration index
							break

						end
					end

					θx1[1] = k # Store current iteration index
					k += 1
				end

		end

	end # solver_v3!

	# Same as solver_v3! but records diagnostic information
	function solver_diagnostic_v3!(x,u,v,w,q,k,αk,βk,σ1,σ2,σ3,θx1,θx2,θu1,θu2,xopt,uopt,rd2o,dye,proj_count_x,proj_count_u,term_crit=true,yref=eX.y,powiter_text=true,ϵp=eX.ϵ_primal,ϵd=eX.ϵ_dual)

		# σ1 = powiter!(utils.xL,utils.uL,utils.vL,αk,βk,σ1,σ2,θx1,θx2)

		if powiter_text
			print("Power iteration estimate of σ  : $σ1\nSVD estimate of σ              : $(utils.σ)")
		end
			
		σ1 = eX.μ/(2.2*σ1) 						# μ/(2σ) (make the σ estimate an over approximation)
		σ2 = eX.μ + 2*eX.λ 						# μ + 2λ

		for t = 1:eX.N-1	
			θx1 .= q[t]
			lmul!(2.0*σ1,θx1)
			v[t] .= w[t] .+ θx1
		end

		# PIPG iteration
		k = 1.0 # Iteration counter
		if eX.uref == nothing

				while k ≤ eX.kmax_pipg
					αk = 2/(k*eX.μ + σ2)				# PIPG step size

					σ3 = 0.0							# For primal stop crit.
					βk = 0.0 							# For primal stop crit.
					for t = 1:eX.N-1
						mul!(θu1,eX.mBdT,v[t])
						mul!(θu2,eX.R,u[t])
						θu1 .= θu1 .+ θu2
						lmul!(-αk,θu1)
						θu2 .= u[t] .+ θu1
						θu1 .= u[t] 								# For primal stop crit.
						# eX.project_u_diagnostic!(u[t],θu2,t)		# Projection onto input constraint set

						θx1 .= x[t+1] .- yref[t+1]
						if t < eX.N-1
							mul!(θx2,eX.Q,θx1)
						else
							mul!(θx2,eX.Qf,θx1)
						end
						mul!(θx1,eX.mAdT,v[t+1])
						θx1 .= θx1 .+ v[t] .+ θx2
						lmul!(-αk,θx1)
						θx2 .= x[t+1] .+ θx1
						θx1 .= x[t+1] 								# For primal stop crit.
						# eX.project_x_diagnostic!(x[t+1],θx2,t)		# Projection onto state constraint set

						eX.project_xu_diagnostic!(x[t+1],u[t],θx2,θu2,t)

						θx1 .= x[t+1] .- θx1						# For primal stop crit.
						θu1 .= u[t] .- θu1							# For primal stop crit.

						# Count no. of projections on input constraints set
						proj_count_x[Int64(k)] += eX.proj_counter[1]

						# Count no. of projections on input constraints set
						proj_count_u[Int64(k)] += eX.proj_counter[2]

						# Computation of primal stopping criteria
						# Relative change of primal variable: |zkp1-zk|/|zkp1| 
						for j = 1:eX.nx
							σ3 = max(σ3,abs(θx1[j]))
							βk = max(βk,abs(x[t+1][j]))
						end
						for j = 1:eX.nu
							σ3 = max(σ3,abs(θu1[j]))
							βk = max(βk,abs(u[t][j]))
						end
					end
					σ3 = σ3/βk							# Primal stop crit. value

					# Normalize projection count by horizon length
					proj_count_u[Int64(k)] = proj_count_u[Int64(k)]/(eX.N-1)
					proj_count_x[Int64(k)] = proj_count_x[Int64(k)]/(eX.N-1)

					βk = 0.0							# For dual stop crit.
					αk = 0.0							# For dual stop crit.
					for t = 1:eX.N-1
						mul!(θx1,eX.Ad,x[t])
						mul!(θx2,eX.Bd,u[t])		
						q[t] .= x[t+1] .- θx1 .- θx2 .- eX.gd
						axpy!((k+1.0)*σ1,q[t],w[t])		# PIPG step size: (k+1)μ/(2σ)
						θx1 .= q[t]
						lmul!((k+2.0)*σ1,θx1)			# PIPG step size: (k+2)μ/(2σ)
						v[t] .= w[t] .+ θx1	

						# computation of dual stopping criteria
						# |G*zkp1-g|/|wkp1|
						for j=1:eX.nx
							βk = max(βk,abs(q[t][j]))
							αk = max(αk,abs(w[t][j]))
						end
					end

					if term_crit
						# Primal dual stopping criteria
						if σ3 ≤ ϵp		  							# Primal test
							if βk ≤ max(eX.ϵ_abs,ϵd*αk)				# Dual test

								# Current iterate is acceptable
								θx1[1] = k+1 # Store current iteration index
								# print("\nTotal PIPG iterations        : $(k+1)")
								break

							end
						end
					end

					rd2o[Int64(k)],dye[Int64(k)] = utils.compute_error(x,u,xopt,uopt,false)

					k += 1

					## Diagnostics
					# print("\nIter:\n")
					# print(Array(x))
					# print("\n")
					# print(Array(u))
					# print("\n")
				end

		else

				while k ≤ eX.kmax_pipg
					αk = 2/(k*eX.μ + σ2)				# PIPG step size

					σ3 = 0.0							# For primal stop crit.
					βk = 0.0 							# For primal stop crit.
					for t = 1:eX.N-1
						θu2 .= u[t] .- eX.uref[t]
						mul!(θu1,eX.R,θu2)
						mul!(θu2,eX.mBdT,v[t])
						θu1 .= θu1 .+ θu2
						lmul!(-αk,θu1)
						θu2 .= u[t] .+ θu1
						θu1 .= u[t] 								# For primal stop crit.
						# eX.project_u_diagnostic!(u[t],θu2,t)		# Projection onto input constraint set

						θx1 .= x[t+1] .- yref[t+1]
						if t < eX.N-1
							mul!(θx2,eX.Q,θx1)
						else
							mul!(θx2,eX.Qf,θx1)
						end
						mul!(θx1,eX.mAdT,v[t+1])
						θx1 .= θx1 .+ v[t] .+ θx2
						lmul!(-αk,θx1)
						θx2 .= x[t+1] .+ θx1
						θx1 .= x[t+1] 								# For primal stop crit.
						# eX.project_x_diagnostic!(x[t+1],θx2,t)		# Projection onto state constraint set

						eX.project_xu_diagnostic!(x[t+1],u[t],θx2,θu2,t)

						θx1 .= x[t+1] .- θx1						# For primal stop crit.
						θu1 .= u[t] .- θu1							# For primal stop crit.

						# Count no. of projections on input constraints set
						proj_count_x[Int64(k)] += eX.proj_counter[1]

						# Count no. of projections on input constraints set
						proj_count_u[Int64(k)] += eX.proj_counter[2]		

						# Computation of primal stopping criteria
						# Relative change of primal variable: |zkp1-zk|/|zkp1| 
						for j = 1:eX.nx
							σ3 = max(σ3,abs(θx1[j]))
							βk = max(βk,abs(x[t+1][j]))
						end
						for j = 1:eX.nu
							σ3 = max(σ3,abs(θu1[j]))
							βk = max(βk,abs(u[t][j]))
						end
					end
					σ3 = σ3/βk							# Primal stop crit. value

					# Normalize projection count by horizon length
					proj_count_u[Int64(k)] = proj_count_u[Int64(k)]/(eX.N-1)
					proj_count_x[Int64(k)] = proj_count_x[Int64(k)]/(eX.N-1)

					βk = 0.0							# For dual stop crit.
					αk = 0.0							# For dual stop crit.
					for t = 1:eX.N-1
						mul!(θx1,eX.Ad,x[t])
						mul!(θx2,eX.Bd,u[t])		
						q[t] .= x[t+1] .- θx1 .- θx2 .- eX.gd
						axpy!((k+1.0)*σ1,q[t],w[t])		# PIPG step size: (k+1)μ/(2σ)
						θx1 .= q[t]
						lmul!((k+2.0)*σ1,θx1)			# PIPG step size: (k+2)μ/(2σ)
						v[t] .= w[t] .+ θx1	

						# Computation of dual stopping criteria
						# |G*zkp1-g|/|wkp1|
						for j=1:eX.nx
							βk = max(βk,abs(q[t][j]))
							αk = max(αk,abs(w[t][j]))
						end
					end

					if term_crit
						# Primal dual stopping criteria
						if σ3 ≤ ϵp										# Primal test
							if βk ≤ max(eX.ϵ_abs,ϵd*αk)	 				# Dual test

								# Current iterate is acceptable
								θx1[1] = k+1 # Store current iteration index
								# print("\nTotal PIPG iterations        : $(k+1)")
								break

							end
						end
					end

					rd2o[Int64(k)],dye[Int64(k)] = utils.compute_error(x,u,xopt,uopt,false)

					k += 1

					## Diagnostics
					# print("\nIter:\n")
					# print(Array(x))
					# print("\n")
					# print(Array(u))
					# print("\n")
				end

		end
		
		if k > eX.kmax_pipg-1 
			θx1[1] = k-1 # = eX.kmax_pipg
			# print("\nTotal PIPG iterations        : $(eX.kmax_pipg)")
		end
	end # solver_diagnostic_v3!

end


# ..:: DIAGNOSTICS ::..
# ^^^^^^^^^^^^^^^^^^^^^
#     @code_warntype pipg.solver!(iX.x,iX.u,iX.v,iX.w,iX.q,iX.γk,iX.γ4,iX.γ5,iX.γ1,iX.γ2,iX.γ3,
#                 iX.κ1x,iX.κ2x,iX.κ1u,iX.κ2u)
#     @code_warntype pipg.solver_diagnostic!(iX.x,iX.u,iX.v,iX.w,iX.q,iX.γk,iX.γ4,iX.γ5,iX.γ1,iX.γ2,iX.γ3,
#                             iX.κ1x,iX.κ2x,iX.κ1u,iX.κ2u,iX.xopt,iX.uopt,iX.rd2o,iX.dye)
#     @benchmark eX.project_u!($(iX.u[1]),$(iX.u[2]),$1)
#     @benchmark eX.project_x!($(iX.x[1]),$(iX.x[2]),$1)
#     @code_warntype eX.project_x!(iX.x[1],iX.x[2],1)
#     @code_warntype eX.project_u!(iX.u[1],iX.u[2],1)
