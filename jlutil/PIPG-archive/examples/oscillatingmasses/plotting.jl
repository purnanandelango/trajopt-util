module plotter
using LinearAlgebra
using Plots, LaTeXStrings
pgfplotsx()
# This ensures that a legend entry is not created by default
default(lab="",markersize=2,markerstrokewidth=0.1,xtickfontsize=12,ytickfontsize=12,
	ztickfontsize=12,legendfontsize=9)
using ..eX
using ..utils
using ..iX

function pos_vel_force(x,u)
	Nx = length(x)
	Nu = length(u)
	nrm_pos = [norm(x[t][eX.ix1],Inf) for t in 1:Nx]
	nrm_vel = [norm(x[t][eX.ix2],Inf) for t in 1:Nx]
	nrm_force = [norm(u[t],Inf) for t in 1:Nu]
	p1 = plot(0:(eX.Δ):((Nx-1)*eX.Δ),nrm_pos,line=:solid,color=:red,title=L"\|p\|_{\infty}",ylim=[0,1.1*eX.pmax])
	plot!(0:(eX.Δ):((Nx-1)*eX.Δ),eX.pmax .* ones(Nx),line=:dash,lw=2,color=:black,lab="")
	p2 = plot(0:(eX.Δ):((Nx-1)*eX.Δ),nrm_vel,line=:solid,color=:green,title=L"\|v\|_{\infty}",ylim=[0,1.1*eX.vmax])
	plot!(0:(eX.Δ):((Nx-1)*eX.Δ),eX.vmax .* ones(Nx),line=:dash,lw=2,color=:black,lab="")
	p3 = plot(0:(eX.Δ):((Nu-1)*eX.Δ),nrm_force,line=:solid,color=:blue,title=L"\|u\|_{\infty}",ylim=[0,1.1*eX.umax])
	plot!(0:(eX.Δ):((Nu-1)*eX.Δ),eX.umax .* ones(Nu),line=:dash,lw=2,color=:black,lab="")

	λ = @layout [a b c]
	plot(p1,p2,p3,size=0.87*[999,333],layout=λ)
end


end