function x = generate_grid(xmin,xmax,N,flg)
% Wrapper function for generating 1D grid with different spacing
    switch flg
        case 'uniform'
            x = linspace(xmin,xmax,N);
        case 'cosine'
            x = cosinespace(xmin,xmax,N);
        case 'sine-plus'
            x = sinespace(xmin,xmax,N,'plus');
        case 'sine-minus'
            x = sinespace(xmin,xmax,N,'minus');
        otherwise
            error("Incorrect flag for generating grid.");
    end
end