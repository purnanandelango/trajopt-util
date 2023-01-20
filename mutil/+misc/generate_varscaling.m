function [Sz,cz] = generate_varscaling(zbnd,sclbnd)
% Construct parameters of affine transformation which scale the dimensional primal variables
% z = Sz*z_scl + cz
% Scaling parameters vary across node points

    n = length(zbnd);
    N = size(zbnd{1},3);
    Sz = cell(N,n);
    cz = cell(N,n);

    sclmin = sclbnd(1);
    sclmax = sclbnd(2);    

    for j = 1:n
        nz = size(zbnd{j},1);        
        A = [eye(nz)*sclmin eye(nz);
             eye(nz)*sclmax eye(nz)];    
        for l = 1:N            
            zmin = zbnd{j}(:,1,l);
            zmax = zbnd{j}(:,2,l); 
    
            for k = 1:nz
                assert(zmin(k) < zmax(k),'Dimensional lower bound must be **strictly smaller** than the dimensional upper bound.');
            end
   
            b = A\[zmin;zmax];
            Sz{l,j} = diag(b(1:nz));
            cz{l,j} = b(nz+1:2*nz);
        end
    end    

end