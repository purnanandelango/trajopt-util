function [H_proj,g_proj] = project_polyhed2dims(H,g,dims_lower)
% Compute shadow (projection) of a polytope { z | H*z <= g } onto lower dimensions 
% specified by dims_lower
    P = Polyhedron('A',H,'b',g);
    P_proj = P.projection(dims_lower,'fourier');
    P_proj.minHRep; % Computer minimal H-Rep
    H_proj = P_proj.A;
    g_proj = P_proj.b;
end
