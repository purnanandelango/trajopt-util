clearvars
clc
close all

load solfile.mat

a_blk = 0.1;

x = zeros(nx,N);
yy = zeros(nx,N);
u = zeros(nu,N-1);
x(:,1) = xm{1};
yy(:,1) = ym{1};
for j=1:N-1
	x(:,j+1) = xm{j+1};
	u(:,j) = um{j};
	yy(:,j+1) = ym{j+1};
end

figure
plot3(x(1,:),x(2,:),x(3,:),'-k','DisplayName','Optimal','LineWidth',1);
hold on
% plot3(yy(1,:),yy(2,:),yy(3,:),'-.k','DisplayName','Reference');
% legend('AutoUpdate','off','location','northeast')
% axis equal
axis('manual')

lw = 2;
scl = 0.6;
scl_box = 2;
F((N-1)/4) = struct('cdata',[],'colormap',[]);
cntr = 0;
for j=1:1:N-1
    start_pos = x(1:3,j) + scl_box*s1 - scl*u(1:3,j); 
    end_pos = x(1:3,j) + scl_box*s1;
    h1 = arrow(start_pos,end_pos,'tipangle',30,'Color',[0,0,1],'Width',2,'Length',5);

    start_pos = x(1:3,j) + scl_box*s2 - scl*u(4:6,j); 
    end_pos = x(1:3,j) + scl_box*s2;
    h2 = arrow(start_pos,end_pos,'tipangle',30,'Color',[0,1,0],'Width',2,'Length',5); 
    
    start_pos = x(1:3,j) + scl_box*s3 - scl*u(7:9,j); 
    end_pos = x(1:3,j) + scl_box*s3;
    h3 = arrow(start_pos,end_pos,'tipangle',30,'Color',[1,0,0],'Width',2,'Length',5);

    start_pos = x(1:3,j) + scl_box*s4 - scl*u(10:12,j); 
    end_pos = x(1:3,j) + scl_box*s4;
    h4 = arrow(start_pos,end_pos,'tipangle',30,'Color',[1,0,1],'Width',2,'Length',5);
    
    cube_origin = x(1:3,j)' - scl_box*a_blk*ones(1,3);
    c1 = plotcube(scl_box*2*a_blk*ones(1,3),cube_origin,0.3,[1,0.8,0]);

    xlabel('$x$ [m]','FontSize',25);
    ylabel('$y$ [m]','FontSize',25);
    zlabel('$z$ [m]','FontSize',25);
    view(34,15)
    ax = gca;
    ax.DataAspectRatio = [1,1,1];
    ax.PlotBoxAspectRatio = [1,1,1];
    ax.XLim = [-2.5,2.5];
    ax.YLim = [-5,4.5];
    ax.ZLim = [-0.5,3];
    
    grid on

    cntr = cntr + 1;
    F(cntr+1) = getframe(gcf);
    
    if j < N-1
        delete([h1;h2;h3;h4;c1]);
    end
end
F = F(2:end);
v = VideoWriter('sol.mp4','MPEG-4');
v.FrameRate = 10;
open(v);
writeVideo(v,F)
close(v);