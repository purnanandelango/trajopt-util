module plotter
using Plots, LaTeXStrings, MATLAB
pgfplotsx()
# This ensures that a legend entry is not created by default
default(lab="",markersize=2,markerstrokewidth=0.1,xtickfontsize=12,ytickfontsize=12,
	ztickfontsize=12,legendfontsize=9)
using ..eX

function trajectory3D(x,u)

	rx = [x[t][1] for t = 1:eX.N]
	ry = [x[t][2] for t = 1:eX.N]
	rz = [x[t][3] for t = 1:eX.N]

	u_scl = -0.02
	ux = (u_scl .* [u[t][1] for t =1:eX.N-1]) .+ rx[1:eX.N-1]
	uy = (u_scl .* [u[t][2] for t =1:eX.N-1]) .+ ry[1:eX.N-1]
	uz = (u_scl .* [u[t][3] for t =1:eX.N-1]) .+ rz[1:eX.N-1]

	rylim = max([abs(ry[t]) for t=1:eX.N]...)
	rxlim = max([abs(rx[t]) for t=1:eX.N]...)
	ryxlim = 1.1*max(rylim,rxlim)
	rzlim = 1.1*max([abs(rz[t]) for t=1:eX.N]...)

	p1 = plot(rx,ry,rz,title="Trajectory",lab="",linealpha=0.01,marker=:circle,color=:red,camera=[-36,18])
	plot!(xlabel = L"x\ \mathrm{[m]}",ylabel = L"y\ \mathrm{[m]}",zlabel = L"z\ \mathrm{[m]}")
	for t=1:eX.N-1
		plot!([rx[t],ux[t]],[ry[t],uy[t]],[rz[t],uz[t]],line=:solid,lw=1.1,color=:blue)
	end
	plot!(xlim=[-ryxlim,ryxlim],ylim=[-ryxlim,ryxlim],zlim=[-0.08*rzlim,rzlim])

	display(p1)

end

function trajectory3D_mpc(x,u)

	M = eX.NN-eX.N+1

	rx = [x[t][1] for t = 1:M]
	ry = [x[t][2] for t = 1:M]
	rz = [x[t][3] for t = 1:M]

	u_scl = -0.02
	ux = (u_scl .* [u[t][1] for t =1:M]) .+ rx[1:M]
	uy = (u_scl .* [u[t][2] for t =1:M]) .+ ry[1:M]
	uz = (u_scl .* [u[t][3] for t =1:M]) .+ rz[1:M]

	rylim = max([abs(ry[t]) for t=1:M]...)
	rxlim = max([abs(rx[t]) for t=1:M]...)
	ryxlim = 1.1*max(rylim,rxlim)
	rzlim = 1.1*max([abs(rz[t]) for t=1:M]...)

	p1 = plot(rx,ry,rz,title="Trajectory",lab="",linealpha=0.01,marker=:circle,color=:red,camera=[-36,18])
	plot!(xlabel = L"x\ \mathrm{[m]}",ylabel = L"y\ \mathrm{[m]}",zlabel = L"z\ \mathrm{[m]}")
	for t=1:M
		plot!([rx[t],ux[t]],[ry[t],uy[t]],[rz[t],uz[t]],line=:solid,lw=1.1,color=:blue)
	end
	plot!(xlim=[-ryxlim,ryxlim],ylim=[-ryxlim,ryxlim],zlim=[-0.08*rzlim,rzlim])

	display(p1)

end

function visualization_matlab(X,U)

	N = length(X)
	xtmp = [randn(eX.nx) for _ in 1:N]
	xtmp .= X 
	utmp = [randn(eX.nu) for _ in 1:length(U)]
	utmp .= U

	xm = mxarray(xtmp)
	um = mxarray(utmp)

	mat"""
	setfig
	scl_val = -0.02; % thrust vector scaling

	N = double($(N));
	Del = double($(eX.Δ));

	tvec = 0:Del:Del*(N-1);

	x = zeros($(eX.nx),$(N));
	u = zeros($(eX.nu),$(N)-1);
	x(:,1) = $(xm){1};
	for j=1:$(N)-1;
		x(:,j+1) = $(xm){j+1};
		u(:,j) = $(um){j};
	end

	%%%%%%%%%%%%%%%%%%%%%%%%% NOT NECESSARY

	figure
	traj = plot3(x(1,:),x(2,:),x(3,:),'o-b','MarkerSize',7);
	hold on
	for i = 1:2:N-1
	    thrust = quiver3(x(1,i),x(2,i),x(3,i),scl_val*u(1,i),scl_val*u(2,i),scl_val*u(3,i),'-r','LineWidth',2);
	end
	xlabel('\$x\$');
	ylabel('\$y\$');
	zlabel('\$z\$');
	axis equal

	% plot glideslope cone
	r = linspace(0,30,100);
	th = linspace(0,2*pi,100);
	[R,T] = meshgrid(r,th);
	X = R.*cos(T);
	Y = R.*sin(T);
	Z = 0.5*R/tan($(eX.γ_gs));
	cone = surf(X,Y,Z);
	set(cone,'EdgeAlpha',0,'FaceColor',[0.7,0.5,1],'FaceAlpha',0.5)

	% axis square
	% ax = gca;
	% ax.DataAspectRatio = [1,1,1];
	% ax.PlotBoxAspectRatio = [1,1,1];
	axis equal;
	view(-71,14);

	%%%%%%%%%%%%%%%%%%%%%%%%%


	figure
	subplot(2,3,[1,2,4,5])
	traj = plot3(x(1,:),x(2,:),x(3,:),'o-b','MarkerSize',7);
	hold on
	for i = 1:N-1
	    thrust = quiver3(x(1,i),x(2,i),x(3,i),scl_val*u(1,i),scl_val*u(2,i),scl_val*u(3,i),'-r','LineWidth',2);
	end
	xlabel('\$x\$');
	ylabel('\$y\$');
	zlabel('\$z\$');
	title('3D point-mass soft landing')

	% plot glideslope cone
	r = linspace(0,30,100);
	th = linspace(0,2*pi,100);
	[R,T] = meshgrid(r,th);
	X = R.*cos(T);
	Y = R.*sin(T);
	Z = R/tan($(eX.γ_gs));
	cone = surf(X,Y,Z);
	set(cone,'EdgeAlpha',0,'FaceColor',[0.7,0.5,1],'FaceAlpha',0.5)


	% plot gimbal cone
	r = linspace(0,1,100);
	th = linspace(0,2*pi,100);
	[R,T] = meshgrid(r,th);
	for i = 1:3
	    X = R.*cos(T) + x(1,i);
	    Y = R.*sin(T) + x(2,i);
	    Z = -R/tan($(eX.θ_tl)) + x(3,i);
	    cone2 = surf(X,Y,Z);
	    set(cone2,'EdgeAlpha',0,'FaceColor',[0.5,0.7,1],'FaceAlpha',0.5)
	end

	legend([traj,thrust,cone,cone2],{'Trajectory','Thrust vector','Glide slope cone','Gimbal cone'},'Location',[0.3,0.3,0.03,0.02])

	axis square
	ax = gca;
	ax.DataAspectRatio = [1,1,1];
	ax.PlotBoxAspectRatio = [1,1,1];
	view(-71,14)

	Vmax = $(eX.vmax*eX.scl_x2);

	subplot(2,3,3)
	plot(tvec,Vmax*ones(1,N),'--k','LineWidth',2);
	hold on
	plot(tvec,sqrt(x(4,:).^2 + x(5,:).^2 + x(6,:).^2),'o-r')
	ylim([0,1.1*Vmax])
	xlabel('\$t\$');
	title('Velocity magnitude');
	legend('\$V_{\\max}\$')

	umax = $(eX.umax*eX.scl_u);
	umin = $(eX.umin*eX.scl_u);

	Nu = length(u(1,:));

	subplot(2,3,6)
	hold on
	p_umax = plot(tvec,umax*ones(1,N),'--k','LineWidth',2);
	p_umin = plot(tvec,umin*ones(1,N),'--r','LineWidth',2);
	stairs(tvec(1:Nu),sqrt(u(1,:).^2 + u(2,:).^2 + u(3,:).^2),'o-b','MarkerSize',10)
	xlabel('\$t\$');
	title('Thrust magnitude \$\\|u\\|_2\$');
	legend([p_umax,p_umin],{'\$u_{\\max}\$','\$u_{\\min}\$'})

	"""

end

end