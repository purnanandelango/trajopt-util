function plot_soln(x,u,yy,N,s1,s2,s3,s4,a)
    figure
    plot3(x(1,:),x(2,:),x(3,:),'.-k','DisplayName','Optimal','MarkerSize',20);
    hold on
    plot3(yy(1,:),yy(2,:),yy(3,:),'-.k','DisplayName','Reference');
    legend('AutoUpdate','off','location','northeast')
    axis equal
    axis('manual')

    lw = 1;
    scl = 0.5;
    scl_box = 2;
    for j=1:2:N-1
        start_pos = x(1:3,j) + scl_box*s1 - scl*u(1:3,j); 
        end_pos = x(1:3,j) + scl_box*s1;
        arrow(start_pos,end_pos,'tipangle',30,'Color',[0,0,1],'Width',lw,'Length',5)

        start_pos = x(1:3,j) + scl_box*s2 - scl*u(4:6,j); 
        end_pos = x(1:3,j) + scl_box*s2;
        arrow(start_pos,end_pos,'tipangle',30,'Color',[0,1,0],'Width',lw,'Length',5) 

        start_pos = x(1:3,j) + scl_box*s3 - scl*u(7:9,j); 
        end_pos = x(1:3,j) + scl_box*s3;
        arrow(start_pos,end_pos,'tipangle',30,'Color',[1,0,0],'Width',lw,'Length',5)
        
        start_pos = x(1:3,j) + scl_box*s4 - scl*u(10:12,j); 
        end_pos = x(1:3,j) + scl_box*s4;
        arrow(start_pos,end_pos,'tipangle',30,'Color',[1,0,1],'Width',lw,'Length',5)
        
        cube_origin = x(1:3,j)' - scl_box*a*ones(1,3);
        plotcube(scl_box*2*a*ones(1,3),cube_origin,0.3,[1,0.8,0])
    end
    xlabel('$x$');
    ylabel('$y$');
    zlabel('$z$');
    view(6,9)

    axis equal
    axis('manual')
end


