module exe_utils
# Module defining convenient wrapper functions that call various helper functions in the correct sequence for solving 
# the problem using PIPG and JuMP
# Wrapper functions for benchmarking the performance of PIPG using BenchmarkTools and TimerOutputs are also provided
# Order in which these modules are loaded is important

using ..eX
using ..utils
using ..common_plotter
using ..pipg
using ..iX
using BenchmarkTools, TimerOutputs

function execute_JuMP_proxy!(xopt,uopt,rd2o,dye,proj_count_x,proj_count_u,slvr,ϵ_pd=eX.ϵ_pd_JuMP,ϵ_gap=eX.ϵ_gap_JuMP)
# Obtain optimal solution using a solver in JuMP defined in eX

	# Initialize/reset solution containers
	if ~(rd2o == nothing)
		iX.reset_optvar!(xopt,uopt,rd2o,dye,proj_count_x,proj_count_u)
	end

	# Call JuMP (parse+solve)
	_,compute_time = utils.solve_JuMP!(xopt,uopt,slvr,false,false,eX.x0,eX.y,ϵ_pd,ϵ_gap)

	# Bring the solution to original scale
	utils.set_scale!(xopt,uopt,:scale)

	return compute_time

end

execute_JuMP!(slvr,ϵ_pd=eX.ϵ_pd_JuMP,ϵ_gap=eX.ϵ_gap_JuMP) = execute_JuMP_proxy!(iX.xopt,iX.uopt,iX.rd2o,iX.dye,iX.proj_count_x,iX.proj_count_u,slvr,ϵ_pd,ϵ_gap)

execute_JuMP_benchmark(slvr,ϵ_pd=eX.ϵ_pd_JuMP,ϵ_gap=eX.ϵ_gap_JuMP) = begin
	bench_var = @benchmarkable execute_JuMP_proxy!($(iX.xopt),$(iX.uopt),$(iX.rd2o),$(iX.dye),$(iX.proj_count_x),$(iX.proj_count_u),$(slvr),$(ϵ_pd),$(ϵ_gap))
	run(bench_var)
end

function execute_PIPG_proxy!(x,u,v,w,q,xL,uL,vL,xopt,uopt,rd2o,dye,proj_count_x,proj_count_u,diagnose_flag=false,σ1=0.0,ϵp=eX.ϵ_primal,ϵd=eX.ϵ_dual)
# Executes PIPG
# If max singular value of G ( in the linear equality constraint Gx=g ) is provided (i.e. σ1 == 0) then power iteration is not called
# Before using the diagnostics mode, make sure that solution obtained using JuMP is available
# Remember that the JuMP solution is in original scale

	if σ1 ≈ 0.0
		# Initialize/reset PIPG and power iteration solution variables
		iX.reset_pipg_powiter!(x,u,v,w,q,xL,uL,vL)
		# Call power iteration
		σ1 = pipg.powiter!(xL,uL,vL,iX.γ4,iX.γ5,iX.γ1,iX.γ2,iX.κ1x,iX.κ2x)
	else
		iX.reset_pipg!(x,u,v,w,q)
	end

	# Call PIPG
	if diagnose_flag
		# Scale the JuMP solution
		utils.set_scale!(xopt,uopt,:unscale)
		pipg.solver_diagnostic!(x,u,v,w,q,iX.γk,iX.γ4,iX.γ5,σ1,iX.γ2,iX.γ3,iX.κ1x,iX.κ2x,iX.κ1u,iX.κ2u,xopt,uopt,rd2o,dye,proj_count_x,proj_count_u,false,eX.y,true,ϵp,ϵd)
		print("\n\n")

		# Bring both JuMP and PIPG solutions to original scale
		utils.set_scale!(x,u,:scale)
		utils.set_scale!(xopt,uopt,:scale);	
	else
		pipg.solver!(x,u,v,w,q,iX.γk,iX.γ4,iX.γ5,σ1,iX.γ2,iX.γ3,iX.κ1x,iX.κ2x,iX.κ1u,iX.κ2u,eX.y,ϵp,ϵd)

		# Bring PIPG solution to original scale
		utils.set_scale!(x,u,:scale)
	end

	# Report solution statistics
	rd2o_final,dye_final = utils.compute_error(x,u,xopt,uopt,true,:unscaled,iX.get_pipg_itercount())

	# If diagnostics is enabled, plot relative d2o and dye
	if diagnose_flag
		common_plotter.solution_quality(rd2o,dye,iX.get_pipg_itercount())
		common_plotter.projection_count(proj_count_x,proj_count_u,eX.proj_count_limit...,iX.get_pipg_itercount())
	end

	return rd2o_final,dye_final,iX.get_pipg_itercount()

end

execute_PIPG!(diagnose_flag=false,σ1=0.0,ϵp=eX.ϵ_primal,ϵd=eX.ϵ_dual) = execute_PIPG_proxy!(iX.x,iX.u,iX.v,iX.w,iX.q,utils.xL,utils.uL,utils.vL,iX.xopt,iX.uopt,iX.rd2o,iX.dye,iX.proj_count_x,iX.proj_count_u,diagnose_flag,σ1,ϵp,ϵd)


execute_PIPG_benchmark(ϵp=eX.ϵ_primal,ϵd=eX.ϵ_dual) = begin
	σ1 = pipg.powiter!(utils.xL,utils.uL,utils.vL,iX.γ4,iX.γ5,iX.γ1,iX.γ2,iX.κ1x,iX.κ2x)
	bench_var = @benchmarkable pipg.solver!($(iX.x),$(iX.u),$(iX.v),$(iX.w),$(iX.q),$(iX.γk),$(iX.γ4),$(iX.γ5),$σ1,$(iX.γ2),$(iX.γ3),$(iX.κ1x),$(iX.κ2x),$(iX.κ1u),$(iX.κ2u),$(eX.y),$(ϵp),$(ϵd)) setup=(iX.reset_pipg!(iX.x,iX.u,iX.v,iX.w,iX.q))
	run(bench_var)
end

execute_PIPG_timeit(ϵp=eX.ϵ_primal,ϵd=eX.ϵ_dual) = begin
		σ1 = pipg.powiter!(utils.xL,utils.uL,utils.vL,iX.γ4,iX.γ5,iX.γ1,iX.γ2,iX.κ1x,iX.κ2x) 
		to_var = TimerOutput()
		enable_timer!(to_var)
		reset_timer!(to_var)
		for _ in 1:20
		    iX.reset_pipg!(iX.x,iX.u,iX.v,iX.w,iX.q)
		   @timeit to_var "PIPG" pipg.solver!(iX.x,iX.u,iX.v,iX.w,iX.q,iX.γk,iX.γ4,iX.γ5,σ1,iX.γ2,iX.γ3,iX.κ1x,iX.κ2x,iX.κ1u,iX.κ2u,eX.y,ϵp,ϵd)
		end
		disable_timer!(to_var)
		show(to_var)
end

# When coupled state and control constraint exist, use PIPG v3
function execute_PIPG_v3_proxy!(x,u,v,w,q,xL,uL,vL,xopt,uopt,rd2o,dye,proj_count_x,proj_count_u,diagnose_flag=false,σ1=0.0,ϵp=eX.ϵ_primal,ϵd=eX.ϵ_dual)
# Executes PIPG : solver_v3
# If max singular value of G ( in the linear equality constraint Gx=g ) is provided (i.e. σ1 == 0) then power iteration is not called
# Before using the diagnostics mode, make sure that solution obtained using JuMP is available
# Remember that the JuMP solution is in original scale

	if σ1 ≈ 0.0
		# Initialize/reset PIPG and power iteration solution variables
		iX.reset_pipg_powiter!(x,u,v,w,q,xL,uL,vL)
		# Call power iteration
		σ1 = pipg.powiter!(xL,uL,vL,iX.γ4,iX.γ5,iX.γ1,iX.γ2,iX.κ1x,iX.κ2x)
	else
		iX.reset_pipg!(x,u,v,w,q)
	end

	# Call PIPG
	if diagnose_flag
		# Scale the JuMP solution
		utils.set_scale!(xopt,uopt,:unscale)
		pipg.solver_diagnostic_v3!(x,u,v,w,q,iX.γk,iX.γ4,iX.γ5,σ1,iX.γ2,iX.γ3,iX.κ1x,iX.κ2x,iX.κ1u,iX.κ2u,xopt,uopt,rd2o,dye,proj_count_x,proj_count_u,false,eX.y,true,ϵp,ϵd)
		print("\n\n")

		# Bring both JuMP and PIPG solutions to original scale
		utils.set_scale!(x,u,:scale)
		utils.set_scale!(xopt,uopt,:scale);	
	else
		pipg.solver_v3!(x,u,v,w,q,iX.γk,iX.γ4,iX.γ5,σ1,iX.γ2,iX.γ3,iX.κ1x,iX.κ2x,iX.κ1u,iX.κ2u,eX.y,ϵp,ϵd)

		# Bring PIPG solution to original scale
		utils.set_scale!(x,u,:scale)
	end

	# Report solution statistics
	rd2o_final,dye_final = utils.compute_error(x,u,xopt,uopt,true,:unscaled,iX.get_pipg_itercount())

	# If diagnostics is enabled, plot relative d2o and dye
	if diagnose_flag
		common_plotter.solution_quality(rd2o,dye,iX.get_pipg_itercount())
		common_plotter.projection_count(proj_count_x,proj_count_u,eX.proj_count_limit...,iX.get_pipg_itercount())
	end

	return rd2o_final,dye_final,iX.get_pipg_itercount()

end

execute_PIPG_v3!(diagnose_flag=false,σ1=0.0,ϵp=eX.ϵ_primal,ϵd=eX.ϵ_dual) = execute_PIPG_v3_proxy!(iX.x,iX.u,iX.v,iX.w,iX.q,utils.xL,utils.uL,utils.vL,iX.xopt,iX.uopt,iX.rd2o,iX.dye,iX.proj_count_x,iX.proj_count_u,diagnose_flag,σ1,ϵp,ϵd)

execute_PIPG_v3_benchmark(ϵp=eX.ϵ_primal,ϵd=eX.ϵ_dual) = begin
	σ1 = pipg.powiter!(utils.xL,utils.uL,utils.vL,iX.γ4,iX.γ5,iX.γ1,iX.γ2,iX.κ1x,iX.κ2x)
	bench_var = @benchmarkable pipg.solver_v3!($(iX.x),$(iX.u),$(iX.v),$(iX.w),$(iX.q),$(iX.γk),$(iX.γ4),$(iX.γ5),$σ1,$(iX.γ2),$(iX.γ3),$(iX.κ1x),$(iX.κ2x),$(iX.κ1u),$(iX.κ2u),$(eX.y),$(ϵp),$(ϵd)) setup=(iX.reset_pipg!(iX.x,iX.u,iX.v,iX.w,iX.q))
	run(bench_var)
end

execute_PIPG_v3_timeit(ϵp=eX.ϵ_primal,ϵd=eX.ϵ_dual) = begin
		σ1 = pipg.powiter!(utils.xL,utils.uL,utils.vL,iX.γ4,iX.γ5,iX.γ1,iX.γ2,iX.κ1x,iX.κ2x) 
		to_var = TimerOutput()
		enable_timer!(to_var)
		reset_timer!(to_var)
		for _ in 1:20
		    iX.reset_pipg!(iX.x,iX.u,iX.v,iX.w,iX.q)
		   @timeit to_var "PIPG" pipg.solver_v3!(iX.x,iX.u,iX.v,iX.w,iX.q,iX.γk,iX.γ4,iX.γ5,σ1,iX.γ2,iX.γ3,iX.κ1x,iX.κ2x,iX.κ1u,iX.κ2u,eX.y,ϵp,ϵd)
		end
		disable_timer!(to_var)
		show(to_var)
end

end