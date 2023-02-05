module plotter
using LinearAlgebra
using Plots, LaTeXStrings
using MATLAB
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
	x = eX.p1[1] .+ collect(-eX.r1:(2*eX.r1/(100-1)):eX.r1)
	y = eX.p1[2] .+ sqrt.( abs.(eX.r1^2 .- (x .- eX.p1[1]) .^ 2) )
	p1 = plot(x,y,color=:black,line=:solid)
	y = eX.p1[2] .- sqrt.( abs.(eX.r1^2 .- (x .- eX.p1[1]) .^ 2) )
	plot!(x,y,color=:black,line=:solid)

	x = eX.p2[1] .+ collect(-eX.r2:(2*eX.r2/(100-1)):eX.r2)
	y = eX.p2[2] .+ sqrt.( abs.(eX.r2^2 .- (x .- eX.p2[1]) .^ 2) )
	plot!(x,y,color=:black,line=:solid)
	y = eX.p2[2] .- sqrt.( abs.(eX.r2^2 .- (x .- eX.p2[1]) .^ 2) )
	plot!(x,y,color=:black,line=:solid)

	x = eX.p3[1] .+ collect(-eX.r3:(2*eX.r3/(100-1)):eX.r3)
	y = eX.p3[2] .+ sqrt.( abs.(eX.r3^2 .- (x .- eX.p3[1]) .^ 2) )
	plot!(x,y,color=:black,line=:solid)
	y = eX.p3[2] .- sqrt.( abs.(eX.r3^2 .- (x .- eX.p3[1]) .^ 2) )
	plot!(x,y,color=:black,line=:solid)

	plot!([eX.x0_unscl[1]],[eX.x0_unscl[2]],color=:orange,marker=:circle,linealpha=0.01,lab=L"x_0[1:3]")
	plot!([eX.xf_unscl[1]],[eX.xf_unscl[2]],color=:orange,marker=:square,linealpha=0.01,lab=L"x_0[1:3]")

	plot!([y_unscl[t][1] for t=1:eX.N],[y_unscl[t][2] for t=1:eX.N],color=:red,line=:solid,lw=0.7,lab="Cost Ref.")
	plot!([eX.ȳ_unobs[t][1]*eX.scl_x1 for t=1:eX.N],[eX.ȳ_unobs[t][2]*eX.scl_x1 for t=1:eX.N],color=:magenta,line=:solid,lw=1.3,lab="Obstcl. Ref.")
	# lim_set = [-max(eX.r1,eX.r2,eX.r3),max(eX.r1,eX.r2,eX.r3)]
	lim_set = [-20,20]
	plot!(xlim=lim_set,ylim=lim_set)

	if pltopt == 1
		plot!([x1[t][1] for t=1:eX.N],[x1[t][2] for t=1:eX.N],color=:blue,line=:solid,lab="PIPG")
		plot!(size=[500,500],legend=:topleft,title="Position")

		p0 = plot(eX.tvec,[x1[t][3] for t=1:eX.N],color=:orange,line=:solid,lab="PIPG",title="Altitude",size=[400,400])
	end

	if pltopt == 2
		plot!([x1[t][1] for t=1:eX.N],[x1[t][2] for t=1:eX.N],color=:blue,line=:solid,lab="PIPG")
		plot!([iX.xopt[t][1] for t=1:eX.N],[iX.xopt[t][2] for t=1:eX.N],color=:cyan,line=:dash,lab=string("JuMP - ",String(utils.solver_JuMP_choice[1])))
		plot!(size=[500,500],legend=:topleft,title="Position")

		p0 = plot(eX.tvec,[x1[t][3] for t=1:eX.N],color=:orange,line=:solid,lab="PIPG",title="Altitude",size=[400,400])
		plot!(eX.tvec,[iX.xopt[t][3] for t=1:eX.N],color=:green,line=:dash,lab=string("JuMP - ",String(utils.solver_JuMP_choice[1])))
	end 	

	vmax = eX.vmax
	umax = eX.umax

	# Velocity
	p2 = plot(eX.tvec,vmax .* ones(eX.N),color=:black,line=:solid)

	if pltopt==1
		plot!(eX.tvec,[norm(x1[t][4:6]) for t=1:eX.N],color=:blue,line=:solid,lab="PIPG")
	end

	if pltopt==2
		plot!(eX.tvec,[norm(x1[t][4:6]) for t=1:eX.N],color=:blue,line=:solid,lab="PIPG")
		plot!(eX.tvec,[norm(iX.xopt[t][4:6]) for t=1:eX.N],color=:cyan,line=:dash,lab=string("JuMP - ",String(utils.solver_JuMP_choice[1])))
	end

	plot!(size=[500,500],legend=:none,title="Velocity")

	# Acceleration
	p3 = plot(eX.tvec[1:eX.N-1],umax .* ones(eX.N-1),color=:black,line=:solid)

	if pltopt==1
		plot!(eX.tvec[1:eX.N-1],[norm(u1[t]) for t=1:eX.N-1],color=:blue,line=:solid,lab="PIPG")
	end

	if pltopt==2
		plot!(eX.tvec[1:eX.N-1],[norm(u1[t]) for t=1:eX.N-1],color=:blue,line=:solid,lab="PIPG")		
		plot!(eX.tvec[1:eX.N-1],[norm(iX.uopt[t]) for t=1:eX.N-1],color=:cyan,line=:dash,lab=string("JuMP - ",String(utils.solver_JuMP_choice[1])))
	end

	plot!(size=[500,500],legend=:bottomleft,title="Acceleration")

	λ = @layout [a b c]
	display(plot(p1,p2,p3,layout=λ,size=[999,333]))

	display(p0)
end

function trajectory3D(x,u,y)
	Nx = length(x)
	Nu = length(u)
	Ny = length(y)

	xm = [zeros(eX.nx) for _ in 1:Nx]
	um = [zeros(eX.nu) for _ in 1:Nu]
	ym = [zeros(eX.nx) for _ in 1:Ny]
	xm .= x
	um .= u
	ym .= y

	p1m = zeros(3);
	p2m = zeros(3);
	p3m = zeros(3);
	p1m .= eX.p1
	p2m .= eX.p2
	p3m .= eX.p3

	mat"""
	setfig;

	Nx = $(Nx);
	Nu = $(Nu);
	Ny = $(Ny);
	nx = $(eX.nx);
	nu = $(eX.nu);
	x = zeros(nx,Nx);
	u = zeros(nu,Nu);
	y = zeros(nx,Ny);
	for j=1:Nx
		x(:,j) = $(xm){j};
	end
	for j=1:Nu
		u(:,j) = $(um){j};
	end
	for j=1:Ny
		y(:,j) = $(ym){j};
	end

	r1 = $(eX.r1);
	r2 = $(eX.r2);
	r3 = $(eX.r3);
	p1 = $(p1m);
	p2 = $(p2m);
	p3 = $(p3m);

	figure
	[x1,y1,z1] = cylinder();
	z1 = 5*(z1*2-1);

	surf(r1*x1+p1(1),r1*y1+p1(2),z1,'EdgeAlpha',0.1,'FaceAlpha',0.1);
	hold on
	surf(r3*x1+p3(1),r3*y1+p3(2),z1,'EdgeAlpha',0.1,'FaceAlpha',0.1);
	surf(r2*x1+p2(1),r2*y1+p2(2),z1,'EdgeAlpha',0.1,'FaceAlpha',0.1);

	plt1 = plot3(x(1,:),x(2,:),x(3,:),'.-b','LineWidth',2,'DisplayName','Solution');
	plt2 = plot3(y(1,:),y(2,:),y(3,:),'.-r','LineWidth',2,'DisplayName','Reference');
	xlabel('\$x\$');
	ylabel('\$y\$');
	zlabel('\$z\$');
	plt3 = plot(x(1,1),x(2,1),'ob');
	view(-88,46);

	scl_u = -0.05;
	for j=1:Nu
		if j==1
			plt4 = quiver3(x(1,j),x(2,j),x(3,j),scl_u*u(1,j),scl_u*u(2,j),scl_u*u(3,j),'-m','LineWidth',2);
		else
			quiver3(x(1,j),x(2,j),x(3,j),scl_u*u(1,j),scl_u*u(2,j),scl_u*u(3,j),'-m','LineWidth',2);
		end
	end
	legend([plt1,plt2,plt3,plt4],{'Solution','Reference','\$t=0\$','Thrust'});

	title('Quadrotor Planning');

	axis equal;

	"""
end


end