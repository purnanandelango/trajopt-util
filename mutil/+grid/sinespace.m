function x = sinespace(xmin,xmax,N,flg)
% Sine spaced grid with N points
    tau = sin((0:(N-1))*0.5*pi/(N-1)); 
    switch flg
        case 'plus' % high density at xmax
            x = tau*(xmax-xmin) +xmin;            
        case 'minus' % high density at xmin
            x = fliplr((1-tau)*(xmax-xmin) + xmin);            
    end
end