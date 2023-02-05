module iX
using LinearAlgebra, StaticArrays	      # Required packages
using ..eX	  							  # Module with problem definition
using ..utils

# PIPG
	const x = [MVector{eX.nx}(randn(eX.nx)) for _ in 1:eX.N]	   # 0:N-1
	const u = [MVector{eX.nu}(randn(eX.nu)) for _ in 1:(eX.N-1)]   # 0:N-2
	const v = [MVector{eX.nx}(randn(eX.nx)) for _ in 1:eX.N]       # 1:N
	const w = [MVector{eX.nx}(randn(eX.nx)) for _ in 1:(eX.N-1)]   # 1:N-1
	const q = [MVector{eX.nx}(randn(eX.nx)) for _ in 1:(eX.N-1)]   # 1:N-1

	# Special initialization for PIPG
	x[1] .= eX.x0
	v[eX.N] .= zeros(eX.nx)

	# Temp. storage variables
	const γ1 = 100.0		
	const γ2 = 200.0								# Must ensure that γ1 and γ2 are different
	const γ3 = randn()		
	const γ4 = randn()		
	const γ5 = randn()
	const γk = randn() 

	const κ1x = MVector{eX.nx}(randn(eX.nx)) 		# Upon pipg exit, κ1x[1] holds the index of the final iteration 
	const κ2x = MVector{eX.nx}(randn(eX.nx))
	const κ1u = MVector{eX.nu}(randn(eX.nu))
	const κ2u = MVector{eX.nu}(randn(eX.nu))

	function reset_pipg!(xx,uu,vv,ww,qq)
		# Reset the container for PIPG
		for t = 1:eX.N-1
			xx[t+1] .= randn(eX.nx)
			uu[t] .= randn(eX.nu)
			vv[t] .= randn(eX.nx)
			ww[t] .= randn(eX.nx)
		end
		xx[1] .= eX.x0
		vv[eX.N] .= zeros(eX.nx)
		for t = 1:eX.N-1	
			qq[t] .= xx[t+1] .- eX.Ad*xx[t] .- eX.Bd*uu[t] .- eX.gd
		end

		# Reset the temp storage vectors
		# κ1x .= randn(eX.nx)
		# κ2x .= randn(eX.nx)
		# κ1u .= randn(eX.nu)
		# κ2u .= randn(eX.nu)
	end

	function reset_pipg_powiter!(xx,uu,vv,ww,qq,xLL,uLL,vLL)
		reset_pipg!(xx,uu,vv,ww,qq)
		utils.reset_powiter!(xLL,uLL,vLL)
	end

	function get_pipg_itercount()
		# Return the number of PIPG iterations in the last run
		return Int64(κ1x[1])
	end

# Optimal value and solution performance containers
	const xopt = [MVector{eX.nx}(randn(eX.nx)) for _ in 1:eX.N]
	const uopt = [MVector{eX.nu}(randn(eX.nu)) for _ in 1:(eX.N-1)]
	
	rd2o = ones(eX.kmax_pipg)											# Distance to optimality
	dye = ones(eX.kmax_pipg) 											# Dynamics constraint violation
	proj_count_u = zeros(eX.kmax_pipg)									# Counts no. of input constraints set projections at each iteration
	proj_count_x = zeros(eX.kmax_pipg)									# Counts no. of state constraints set projections at each iteration

	function reset_optvar!(xxopt,uuopt,rrd2o,ddye,pproj_count_x,pproj_count_u)
		for t = 1:eX.N-1
			xxopt[t+1] .= randn(eX.nx)
			uuopt[t] .= randn(eX.nu)
		end
		xxopt[1] .= randn(eX.nx)

		rrd2o .= ones(eX.kmax_pipg)
		ddye .= ones(eX.kmax_pipg)
		pproj_count_u .= zeros(eX.kmax_pipg)
		pproj_count_x .= zeros(eX.kmax_pipg)
	end

	const xopt2 = [MVector{eX.nx}(randn(eX.nx)) for _ in 1:eX.N]
	const uopt2 = [MVector{eX.nu}(randn(eX.nu)) for _ in 1:(eX.N-1)]

end