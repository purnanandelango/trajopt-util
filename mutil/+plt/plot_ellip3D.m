function [] = plot_ellip3D(A_grid)
% Plot sequence of 3D ellipsoids given by A_grid 

    % Ambient dimension
    n = size(A_grid{1},1);

    % Number of ellipsoids
    M = length(A_grid);

    % Check if dimension is 3
    assert(n==3);

    % Grid of points on unit sphere
    NN = 100;
    [thet_grid,phi_grid] = ndgrid(linspace(0,pi,NN),linspace(0,2*pi,NN));
    x = sin(thet_grid).*cos(phi_grid);
    y = sin(thet_grid).*sin(phi_grid);
    z = cos(thet_grid);

     for j = 1:M
        ux = x;
        uy = y;
        uz = z;
        for k1 = 1:NN
            for k2 = 1:NN
                u = A_grid{j}*[x(k1,k2);y(k1,k2);z(k1,k2)];
                ux(k1,k2) = u(1);
                uy(k1,k2) = u(2);
                uz(k1,k2) = u(3);
            end
        end
        surf(ux,uy,uz,'FaceColor',[0,0,1],'EdgeColor',[0,0,0.7],'FaceAlpha',0.3,...
            'EdgeAlpha',0.6,...
            'EdgeLighting','gouraud','FaceLighting','gouraud','BackFaceLighting','lit',...
            'AmbientStrength',0.8);
        axis equal        
        hold on
     end
     

end
