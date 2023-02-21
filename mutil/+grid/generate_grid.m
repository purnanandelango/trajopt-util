function x = generate_grid(xmin,xmax,N,flg)
% Wrapper function for generating 1D grid with different spacing
    switch flg
        case 'uniform'
            x = linspace(xmin,xmax,N);
        case 'cosine'
            x = grid.cosinespace(xmin,xmax,N);
        case 'sine-plus'
            x = grid.sinespace(xmin,xmax,N,'plus');
        case 'sine-minus'
            x = grid.sinespace(xmin,xmax,N,'minus');
        otherwise
            error("Incorrect flag for generating grid.");
    end
end