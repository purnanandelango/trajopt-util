module plotter
using Plots, LaTeXStrings
pgfplotsx()
# This ensures that a legend entry is not created by default
default(lab="",markersize=2,markerstrokewidth=0.1,xtickfontsize=12,ytickfontsize=12,
	ztickfontsize=12,legendfontsize=9)
using ..eX
using ..utils
using ..iX

function trajectory2D(x1,u1,pltopt::Int64=0,y_unscl=eX.y_unscl)
# pltopt (plot option)
# plot feasible boundary for all plot options (pltopt)
# 0 : plot ref
# 1 : plot ref and x1
# 2 : plot ref, x1 and iX.xopt

	# Position

	# obstacle regions
	# x = eX.p1[1] .+ collect(-eX.r1:(2*eX.r1/(100-1)):eX.r1)
	# y = eX.p1[2] .+ sqrt.( abs.(eX.r1^2 .- (x .- eX.p1[1]) .^ 2) )
	# p1 = plot(x,y,color=:black,line=:solid)
	# y = eX.p1[2] .- sqrt.( abs.(eX.r1^2 .- (x .- eX.p1[1]) .^ 2) )
	# plot!(x,y,color=:black,line=:solid)

	x = eX.p2[1] .+ collect(-eX.r2:(2*eX.r2/(100-1)):eX.r2)
	y = eX.p2[2] .+ sqrt.( abs.(eX.r2^2 .- (x .- eX.p2[1]) .^ 2) )
	p1 = plot(x,y,color=:black,line=:solid)
	y = eX.p2[2] .- sqrt.( abs.(eX.r2^2 .- (x .- eX.p2[1]) .^ 2) )
	plot!(x,y,color=:black,line=:solid)

	x = eX.p3[1] .+ collect(-eX.r3:(2*eX.r3/(100-1)):eX.r3)
	y = eX.p3[2] .+ sqrt.( abs.(eX.r3^2 .- (x .- eX.p3[1]) .^ 2) )
	plot!(x,y,color=:black,line=:solid)
	y = eX.p3[2] .- sqrt.( abs.(eX.r3^2 .- (x .- eX.p3[1]) .^ 2) )
	plot!(x,y,color=:black,line=:solid)

	plot!([eX.x0_unscl[1]],[eX.x0_unscl[2]],color=:orange,marker=:circle,linealpha=0.01,lab=L"x_0[1:2]")
	plot!([eX.xf_unscl[1]],[eX.xf_unscl[2]],color=:orange,marker=:square,linealpha=0.01,lab=L"x_0[1:2]")

	plot!([y_unscl[t][1] for t=1:eX.N],[y_unscl[t][2] for t=1:eX.N],color=:red,line=:solid,lw=0.7,lab="Cost Ref.")
	plot!([eX.ȳ_unobs[t][1]*eX.scl_x1 for t=1:eX.N],[eX.ȳ_unobs[t][2]*eX.scl_x1 for t=1:eX.N],color=:magenta,line=:solid,lw=1.3,lab="Obstcl. Ref.")
	plot!(xlim=[-2,2],ylim=[-2,2])

	if pltopt == 1
		plot!([x1[t][1] for t=1:eX.N],[x1[t][2] for t=1:eX.N],color=:blue,line=:solid,lab="PIPG")
	end

	if pltopt == 2
		plot!([x1[t][1] for t=1:eX.N],[x1[t][2] for t=1:eX.N],color=:blue,line=:solid,lab="PIPG")
		
		plot!([iX.xopt[t][1] for t=1:eX.N],[iX.xopt[t][2] for t=1:eX.N],color=:cyan,line=:dash,lab=string("JuMP - ",String(utils.solver_JuMP_choice[1])))
	end 	

	plot!(size=[500,500],legend=:topleft,title="Position")

	vmax = eX.scl_x2*eX.vmax
	umax = eX.scl_u*eX.umax

	# Velocity
	x = collect(-vmax:(2*vmax/(100-1)):vmax)
	y = sqrt.( abs.(vmax^2 .- x .^ 2) )
	p2 = plot(x,y,color=:black,line=:solid)
	y = -1 .* sqrt.( abs.(vmax^2 .- x .^ 2) )
	plot!(x,y,color=:black,line=:solid)

	if pltopt==1
		plot!([x1[t][3] for t=1:eX.N],[x1[t][4] for t=1:eX.N],color=:blue,line=:solid)
		plot!([x1[1][3]],[x1[1][4]],color=:blue,marker=:circle,linealpha=0.01,lab=L"x_0[3:4]")
		plot!([x1[eX.N][3]],[x1[eX.N][4]],color=:blue,marker=:square,linealpha=0.01,lab=L"x_{N}[3:4]")		
	end

	if pltopt==2
		plot!([x1[t][3] for t=1:eX.N],[x1[t][4] for t=1:eX.N],color=:blue,line=:solid)
		plot!([x1[1][3]],[x1[1][4]],color=:blue,marker=:circle,linealpha=0.01,lab=L"x_0[3:4]")
		plot!([x1[eX.N][3]],[x1[eX.N][4]],color=:blue,marker=:square,linealpha=0.01,lab=L"x_{N}[3:4]")		
		
		plot!([iX.xopt[t][3] for t=1:eX.N],[iX.xopt[t][4] for t=1:eX.N],color=:cyan,line=:dash)
		plot!([iX.xopt[1][3]],[iX.xopt[1][4]],color=:cyan,marker=:circle,linealpha=0.01,lab=L"x^{\star}_0[3:4]")
		plot!([iX.xopt[eX.N][3]],[iX.xopt[eX.N][4]],color=:cyan,marker=:square,linealpha=0.01,lab=L"x^{\star}_N[3:4]")		
	end

	plot!(size=[500,500],legend=:topleft,title="Velocity")

	# Acceleration
	x = collect(-umax:(2*umax/(100-1)):umax)
	y = sqrt.( abs.(umax^2 .- x .^ 2) )
	p3 = plot(x,y,color=:black,line=:solid)
	y = -1 .* sqrt.( abs.(umax^2 .- x .^ 2) )
	plot!(x,y,color=:black,line=:solid)

	if pltopt==1
		plot!([u1[t][1] for t=1:eX.N-1],[u1[t][2] for t=1:eX.N-1],color=:blue,line=:solid)		
		plot!([u1[1][1]],[u1[1][2]],color=:blue,marker=:circle,linealpha=0.01,lab=L"u_{0}")
		plot!([u1[eX.N-1][1]],[u1[eX.N-1][2]],color=:blue,marker=:square,linealpha=0.01,lab=L"u_{N-1}")
	end

	if pltopt==2
		plot!([u1[t][1] for t=1:eX.N-1],[u1[t][2] for t=1:eX.N-1],color=:blue,line=:solid)		
		plot!([u1[1][1]],[u1[1][2]],color=:blue,marker=:circle,linealpha=0.01,lab=L"u_{0}")
		plot!([u1[eX.N-1][1]],[u1[eX.N-1][2]],color=:blue,marker=:square,linealpha=0.01,lab=L"u_{N-1}")

		plot!([iX.uopt[t][1] for t=1:eX.N-1],[iX.uopt[t][2] for t=1:eX.N-1],color=:cyan,line=:dash)		
		plot!([iX.uopt[1][1]],[iX.uopt[1][2]],color=:cyan,marker=:circle,linealpha=0.01,lab=L"u^{\star}_{0}")
		plot!([iX.uopt[eX.N-1][1]],[iX.uopt[eX.N-1][2]],color=:cyan,marker=:square,linealpha=0.01,lab=L"u^{\star}_{N-1}")
	end

	plot!(size=[500,500],legend=:topleft,title="Acceleration")

	λ = @layout [a b c]
	return plot(p1,p2,p3,layout=λ,size=[999,333])


end


end