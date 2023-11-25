function [] = plot_polyhed(H_grid,g_grid)
% Create 2D or 3D plot of sequence of polyhedron using MPT3
% { z | H_grid{j}*z <= g_grid{j} } for j = 1,...,M

    % Ambient dimension
    n = size(H_grid{1},2);

    % Check if dimension is 2D or 3D
    assert(ismember(n,[2,3]),'Only 2D or 3D polyhedron can be plotted.');

    % Number of polyhedrons 
    M = length(H_grid);

    for j = 1:M

        P = Polyhedron('A',H_grid{j},'b',g_grid{j});
        P.minVRep; % Compute minimal vertices
        P.plot('alpha',0.3,'color','blue','EdgeAlpha',0.4)
        hold on

    end

end
