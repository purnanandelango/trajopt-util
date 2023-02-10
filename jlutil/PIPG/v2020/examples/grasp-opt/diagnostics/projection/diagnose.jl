
include("../../../../toolkit/proj_funcs.jl")
include("../../../../toolkit/linalg_utils.jl")

module eXd
using LinearAlgebra, StaticArrays
using JuMP, ECOS, Gurobi, MosekTools
using ..proj
using ..la_utils

# const slvr = :gurobi
const slvr = :ecos

const nu = 9 
const nu1 = 3
const nu2 = 3
const nu3 = 3
const a_blk = 0.1					  							# half the width of the block (m)
const ϵ_admm = 1e-8 				  							# admm    	exit tolerance
const ρ_admm = 2.0					  							# admm 		step-size
const ν_admm = 1/(ρ_admm+1)
const ϵ_pd_JuMP = 1e-9				  							# JuMP solver primal and dual feasibility tolerance	
const ϵ_gap_JuMP = 1e-9 			  							# JuMP solver primal and dual objective gap tolerance


const s1 = SVector{3}([-a_blk,0.0,0.0]) 							# position vector first contact point
const s2 = SVector{3}([-a_blk,0.0,0.0]) 						#                 second 
const s3 = SVector{3}([a_blk,0.0,0.0])						 	#                 third	

const μ1 = 0.5					      							# coefficient of friction associated with first contact point
const μ2 = 0.5					      							# 										  second
const μ3 = 0.75				      								# 										  third 

const F1max = 4.0												# upper bound on magnitude of first contact force
const F2max = 16.0												#							  second	
const F3max = 22.0												#							  third	
	
const Γrot = [la_utils.skew_mat(s1) la_utils.skew_mat(s2) la_utils.skew_mat(s3)]  
# projection matrix for ker Γrot
const proj_ker_Γrot = la_utils.proj_ker(Γrot,nu)

# temp storage variables
const γ1 = randn()
const γ2 = randn()
const κ1_u1 = MVector{nu1}(randn(nu1))
const κ2_u1 = MVector{nu1}(randn(nu1))
const κ1_u2 = MVector{nu2}(randn(nu2))
const κ2_u2 = MVector{nu2}(randn(nu2))
const κ1_u3 = MVector{nu3}(randn(nu3))
const κ2_u3 = MVector{nu3}(randn(nu3))	
const κ1_u = MVector{nu}(randn(nu))
const κ2_u = MVector{nu}(randn(nu))					# first element stores no. of projections performed in the last call to admm!
const y_u = MVector{nu}(randn(nu))					# container for admm warm-start
const u_u = MVector{nu}(randn(nu))					# container for admm warm-start

function proj_friction_soc_ball!(u,w)
	κ2_u1[1] = w[2]
	κ2_u1[2] = w[3]
	κ2_u1[3] = w[1]
	proj.soc_ball!(κ1_u1,κ2_u1,μ1,F1max,nu1,γ1)
	u[1] = κ1_u1[3]
	u[2] = κ1_u1[1]
	u[3] = κ1_u1[2]

	κ2_u2[1] = w[5]
	κ2_u2[2] = w[6]
	κ2_u2[3] = w[4]
	proj.soc_ball!(κ1_u2,κ2_u2,μ2,F2max,nu2,γ1)
	u[4] = κ1_u2[3]
	u[5] = κ1_u2[1]
	u[6] = κ1_u2[2]

	κ2_u3[1] = w[8]
	κ2_u3[2] = w[9]
	κ2_u3[3] = -w[7]
	proj.soc_ball!(κ1_u3,κ2_u3,μ3,F3max,nu3,γ1)
	u[7] = -κ1_u3[3]
	u[8] = κ1_u3[1]
	u[9] = κ1_u3[2]
end

const fu1(u,z) = proj_friction_soc_ball!(u,z)
const fu2(u,z) = mul!(u,proj_ker_Γrot,z)
const project_u!(u,z) = proj.admm!(u,y_u,u_u,z,fu1,fu2,ρ_admm,ν_admm,ϵ_admm,nu,γ1,γ2,κ1_u,κ2_u)

function projsolv_JuMP!(x,z)
# determine projection of z onto a set by solving 
# an optimization problem defined using JuMP
# ---
# x : projected result
# z : vector to be projected onto desired set

	if slvr == :ecos
		model= Model(ECOS.Optimizer)
		set_optimizer_attribute(model,"printlevel",0)
		set_optimizer_attribute(model,"feastol",ϵ_pd_JuMP)
		set_optimizer_attribute(model,"abstol",ϵ_gap_JuMP)
		set_optimizer_attribute(model,"reltol",ϵ_gap_JuMP)
	elseif slvr == :gurobi
		model= Model(Gurobi.Optimizer)
		set_optimizer_attribute(model,"Presolve",2)
		set_optimizer_attribute(model,"FeasibilityTol",ϵ_pd_JuMP)
		set_optimizer_attribute(model,"OptimalityTol",ϵ_pd_JuMP)
		set_optimizer_attribute(model,"BarConvTol",ϵ_gap_JuMP)
	elseif slvr == :mosek
		model= Model(Mosek.Optimizer)
		set_optimizer_attribute(model,"INTPNT_CO_TOL_DFEAS",ϵ_pd_JuMP)
		set_optimizer_attribute(model,"INTPNT_CO_TOL_PFEAS",ϵ_pd_JuMP)
		set_optimizer_attribute(model,"INTPNT_CO_TOL_REL_GAP",ϵ_gap_JuMP)
	else
		error("Invalid solver choice for JuMP.")
	end
	set_silent(model)
	@variable(model,c[1:nu])
	@variable(model,d ≥ 0) 
	@objective(model,Min,d^2)
	@constraint(model,[d,(z .- c)...] in SecondOrderCone())

	@constraint(model,[μ1*c[1],c[2:3]...] in SecondOrderCone())

	@constraint(model,[μ2*c[4],c[5:6]...] in SecondOrderCone())

	@constraint(model,[-μ3*c[7],c[8:9]...] in SecondOrderCone()) 

	@constraint(model,[F1max,c[1:3]...] in SecondOrderCone())

	@constraint(model,[F2max,c[4:6]...] in SecondOrderCone())

	@constraint(model,[F3max,c[7:9]...] in SecondOrderCone())

	@constraint(model,Γrot*c .== 0)

	optimize!(model)

	@assert (termination_status(model) == MOI.OPTIMAL) || (termination_status(model) == MOI.ALMOST_OPTIMAL) "Problem not solved correctly."


	x .= value.(c)[1:nu]

end

end