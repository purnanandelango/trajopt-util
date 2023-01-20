function x = cosinespace(xmin,xmax,N)
% Cosine spaced grid with N points
    tau = cos((0:(N-1))*pi/(N-1)); 
    x = 0.5*(1-tau)*(xmax-xmin) + xmin;
end