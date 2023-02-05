module plotter
using Plots, LaTeXStrings, MATLAB
pgfplotsx()
# This ensures that a legend entry is not created by default
default(lab="",markersize=2,markerstrokewidth=0.1,xtickfontsize=12,ytickfontsize=12,
	ztickfontsize=12,legendfontsize=9)
using ..eX

function trajectory3D_matlab(x,u,y,N)
	# make sure that x and u are unscaled
	# MATLAB tries to convert array of object to cell array of objects
	# converting array of static array to cell array is not great
	# first convert array of static array to regular array of Float64 array
	# then convert this array of Float64 array to cell array of Float64 array

	xtmp = [randn(eX.nx) for _ in 1:N]
	xtmp .= x 
	ytmp = [randn(eX.nx) for _ in 1:N]
	ytmp .= y
	utmp = [randn(eX.nu) for _ in 1:N-1]
	utmp .= u

	xm = mxarray(xtmp)
	ym = mxarray(ytmp)
	um = mxarray(utmp)
	s1 = mxarray(Array(eX.s1))
	s2 = mxarray(Array(eX.s2))
	s3 = mxarray(Array(eX.s3))

	mat"""
	setfig

	x = zeros($(eX.nx),$(N));
	yy = zeros($(eX.nx),$(N));
	u = zeros($(eX.nu),$(N)-1);
	x(:,1) = $(xm){1};
	yy(:,1) = $(ym){1};
	for j=1:$(N)-1;
		x(:,j+1) = $(xm){j+1};
		u(:,j) = $(um){j};
		yy(:,j+1) = $(ym){j+1};
	end
	s1 = $(s1);
	s2 = $(s2);
	s3 = $(s3);

	figure
	plot3(x(1,:),x(2,:),x(3,:),'.-k','DisplayName','Optimal','MarkerSize',20);
	hold on
	plot3(yy(1,:),yy(2,:),yy(3,:),'-.k','DisplayName','Reference');
	legend('AutoUpdate','off','location','northeast')
	axis equal
	axis('manual')

	lw = 2;
	scl = 0.05;
	scl_box = 2;
	for j=1:2:$(N)-1
	    start_pos = x(1:3,j) + scl_box*s1 - scl*u(1:3,j); 
	    end_pos = x(1:3,j) + scl_box*s1;
	    arrow(start_pos,end_pos,'tipangle',30,'Color',[0,0,1],'Width',1,'Length',5)

	    start_pos = x(1:3,j) + scl_box*s2 - scl*u(4:6,j); 
	    end_pos = x(1:3,j) + scl_box*s2;
	    arrow(start_pos,end_pos,'tipangle',30,'Color',[0,1,0],'Width',1,'Length',5) 
	    
	    start_pos = x(1:3,j) + scl_box*s3 - scl*u(7:9,j); 
	    end_pos = x(1:3,j) + scl_box*s3;
	    arrow(start_pos,end_pos,'tipangle',30,'Color',[1,0,0],'Width',1,'Length',5)
	    
	    cube_origin = x(1:3,j)' - scl_box*$(eX.a_blk)*ones(1,3);
	    plotcube(scl_box*2*$(eX.a_blk)*ones(1,3),cube_origin,0.3,[1,0.8,0])
	end
	xlabel('\$x\$');
	ylabel('\$y\$');
	zlabel('\$z\$');
	view(24,31)

	axis equal
	axis('manual')
	"""
end

function trajectory3D(x,u)

	N = length(x)

	rx = [x[t][1] for t = 1:N]
	ry = [x[t][2] for t = 1:N]
	rz = [x[t][3] for t = 1:N]

	u_scl = -0.05
	u1x = (u_scl .* [u[t][1] for t =1:N-1]) .+ rx[1:N-1]
	u1y = (u_scl .* [u[t][2] for t =1:N-1]) .+ ry[1:N-1]
	u1z = (u_scl .* [u[t][3] for t =1:N-1]) .+ rz[1:N-1]

	u2x = (u_scl .* [u[t][4] for t =1:N-1]) .+ rx[1:N-1]
	u2y = (u_scl .* [u[t][5] for t =1:N-1]) .+ ry[1:N-1]
	u2z = (u_scl .* [u[t][6] for t =1:N-1]) .+ rz[1:N-1]

	u3x = (u_scl .* [u[t][7] for t =1:N-1]) .+ rx[1:N-1]
	u3y = (u_scl .* [u[t][8] for t =1:N-1]) .+ ry[1:N-1]
	u3z = (u_scl .* [u[t][9] for t =1:N-1]) .+ rz[1:N-1]

	rylim = max([abs(ry[t]) for t=1:N]...)
	rxlim = max([abs(rx[t]) for t=1:N]...)
	ryxlim = 1.1*max(rylim,rxlim)
	rzlim = 1.1*max([abs(rz[t]) for t=1:N]...)

	p1 = plot(rx,ry,rz,title="Trajectory",lab="",linealpha=0.01,marker=:circle,color=:red,camera=[6,9])
	plot!(xlabel = L"x\ \mathrm{[m]}",ylabel = L"y\ \mathrm{[m]}",zlabel = L"z\ \mathrm{[m]}")
	for t=1:2:N-1
		plot!([rx[t],u1x[t]],[ry[t],u1y[t]],[rz[t],u1z[t]],line=:solid,lw=1.1,color=:red)
		plot!([rx[t],u2x[t]],[ry[t],u2y[t]],[rz[t],u2z[t]],line=:solid,lw=1.1,color=:blue)
		plot!([rx[t],u3x[t]],[ry[t],u3y[t]],[rz[t],u3z[t]],line=:solid,lw=1.1,color=:green)
	end
	plot!(xlim=[-ryxlim,ryxlim],ylim=[-ryxlim,ryxlim],zlim=[-0.08*rzlim,rzlim])

	display(p1)

end

end