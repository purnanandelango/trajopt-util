function [] = plot_ellip2D(A_grid)
% Plot sequence of 2D ellipsoids given by A_grid

    % Ambient dimension
    n = size(A_grid{1},1);

    % Number of ellipses
    M = length(A_grid);

    % Check if dimension is 2
    assert(n==2);

    % Uniformly spaced points on unit circle
    u_grid = [cos(linspace(0,2*pi,100));
              sin(linspace(0,2*pi,100))];

    ellip = cell(1,M);
    for j = 1:M
        ellip{j} = A_grid{j}*u_grid;
        % Plot ellipse
        plot(ellip{j}(1,:),ellip{j}(2,:),'-b');
        hold on
    end
    axis equal

end
