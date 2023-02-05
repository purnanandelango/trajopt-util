module plotter
using Plots, LaTeXStrings
pgfplotsx()
# This ensures that a legend entry is not created by default
default(lab="",markersize=2,markerstrokewidth=0.1,xtickfontsize=12,ytickfontsize=12,
	ztickfontsize=12,legendfontsize=9)
using ..eX

function trajectory3D(x,u)

	rx = [x[t][1] for t = 1:eX.N]
	ry = [x[t][2] for t = 1:eX.N]
	rz = [x[t][3] for t = 1:eX.N]

	u_scl = 0.05
	u1x = (u_scl .* [u[t][1] for t =1:eX.N-1]) .+ rx[1:eX.N-1]
	u1y = (u_scl .* [u[t][2] for t =1:eX.N-1]) .+ ry[1:eX.N-1]
	u1z = (u_scl .* [u[t][3] for t =1:eX.N-1]) .+ rz[1:eX.N-1]

	u2x = (u_scl .* [u[t][4] for t =1:eX.N-1]) .+ rx[1:eX.N-1]
	u2y = (u_scl .* [u[t][5] for t =1:eX.N-1]) .+ ry[1:eX.N-1]
	u2z = (u_scl .* [u[t][6] for t =1:eX.N-1]) .+ rz[1:eX.N-1]

	rylim = max([abs(ry[t]) for t=1:eX.N]...)
	rxlim = max([abs(rx[t]) for t=1:eX.N]...)
	ryxlim = 1.1*max(rylim,rxlim)
	rzlim = 1.1*max([abs(rz[t]) for t=1:eX.N]...)

	p1 = plot(rx,ry,rz,title="Trajectory",lab="",linealpha=0.01,marker=:circle,color=:red,camera=[6,9])
	plot!(xlabel = L"x\ \mathrm{[m]}",ylabel = L"y\ \mathrm{[m]}",zlabel = L"z\ \mathrm{[m]}")
	for t=1:2:eX.N-1
		plot!([rx[t],u1x[t]],[ry[t],u1y[t]],[rz[t],u1z[t]],line=:solid,lw=1.1,color=:red)
		plot!([rx[t],u2x[t]],[ry[t],u2y[t]],[rz[t],u2z[t]],line=:solid,lw=1.1,color=:blue)
	end
	plot!(xlim=[-ryxlim,ryxlim],ylim=[-ryxlim,ryxlim],zlim=[-0.08*rzlim,rzlim])

	display(p1)

end

end