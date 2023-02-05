clear variables 
clc
close all

problem_data
load soln_dat


figure
plot3(x(1,:),x(2,:),x(3,:),'.-k','DisplayName','Optimal','MarkerSize',20);
hold on
plot3(yy(1,:),yy(2,:),yy(3,:),'-.k','DisplayName','Reference');
legend('AutoUpdate','off','location','northeast')
axis equal
axis('manual')

lw = 2;
scl = 0.02;
scl_box = 2;
for j=1:2:N-1
    start_pos = x(1:3,j) + scl_box*s1 - scl*u(1:3,j); 
    end_pos = x(1:3,j) + scl_box*s1;
    arrow(start_pos,end_pos,'tipangle',30,'Color',[0,0,1],'Width',1,'Length',5)

    start_pos = x(1:3,j) + scl_box*s2 - scl*u(4:6,j); 
    end_pos = x(1:3,j) + scl_box*s2;
    arrow(start_pos,end_pos,'tipangle',30,'Color',[0,1,0],'Width',1,'Length',5) 
    
    cube_origin = x(1:3,j)' - scl_box*a*ones(1,3);
    plotcube(scl_box*2*a*ones(1,3),cube_origin,0.3,[1,0.8,0])
end
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
view(6,9)

axis equal
axis('manual')


