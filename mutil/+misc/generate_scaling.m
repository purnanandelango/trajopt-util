function [Sz,cz] = generate_scaling(zbnd,sclbnd)
% Construct parameters of affine transformation which scale the dimensional primal variables
% z = Sz*z_scl + cz

    n = length(zbnd);
    Sz = cell(1,n);
    cz = cell(1,n);

    sclmin = sclbnd(1);
    sclmax = sclbnd(2);    

    for j = 1:n
        zmin = zbnd{j}(:,1);
        zmax = zbnd{j}(:,2);
    
        nz = size(zbnd{j},1); 

        for k = 1:nz
            assert(zmin(k) < zmax(k),'Dimensional lower bound must be **strictly smaller** than the dimensional upper bound.');
        end

        A = [eye(nz)*sclmin eye(nz);
             eye(nz)*sclmax eye(nz)];
        b = A\[zmin;zmax];
        Sz{j} = diag(b(1:nz));
        cz{j} = b(nz+1:2*nz);
    end    

end