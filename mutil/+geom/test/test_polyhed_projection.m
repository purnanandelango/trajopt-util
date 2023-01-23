% Test computation of projection onto polyhedron
  
clc
clearvars
close all

% Ambient dimension
n = 3;
dims_lower = 1:2;

In = eye(n);

% Random point to project onto polytope
x = randn(n,1);

% Construct polytope

%     % From ellipsoid { A*u | \|u\|_2 <= 1 }
%     V = qr(randn(n,n));
%     A = V*diag(randi([1,3],[1,n]) .^ 2);
%     % Alternate ellipsoid representation { z | z'*inv(Q)*z <= 1 }
%     Q = (In/A')/A;
%     % Create SReachEllipsoid object
%     E = SReachEllipsoid(zeros(n,1),Q);
%     % Create outer approximation
%     P = E.outer_approx(2*n+2^n,true);
%     if n <= 3
%         P.plot('color',[0,0,1],'alpha',0.5,'EdgeColor',[0,0,0.7]);
%         hold on
%         if n == 2
%             E.plot('color',[0,1,0]);
%             plot(x(1),x(2),'.k','MarkerSize',12);
%         else
%             plot3(x(1),x(2),x(3),'.k','MarkerSize',12);
%         end
%         axis equal
%     end
%     H = P.A;
%     g = P.b;

    % From a box { z | -z_max <= z <= z_max}
    V = qr(randn(n));
    z_max = randi([1,5],[n,1]);
    H = [In;-In]*V;
    g = [z_max;z_max];

    if n <= 3
        P = Polyhedron('A',H,'b',g);
        P.plot('color',[0,0,1],'alpha',0.3);
        hold on
        if n == 2
            plot(x(1),x(2),'.k','MarkerSize',12);
        else
            plot3(x(1),x(2),x(3),'.k','MarkerSize',12);
        end        
    end

% Compute projection
tic
y = geom.sign_dist_polyhed(x,H,g,true);
toc

if ismember(n,[2,3])
    if n == 2
        plot(y(1),y(2),'.r','MarkerSize',12);
    elseif n == 3
        plot3(y(1),y(2),y(3),'.r','MarkerSize',12);
    end
    axis equal
    ax = gca;
    ax.DataAspectRatioMode = 'manual';
    ax.DataAspectRatio = [1,1,1];
    ax.PlotBoxAspectRatioMode = 'manual';
    ax.PlotBoxAspectRatio = [1,1,1];
end
figure
% Project polytope to lower dimensions
if n >= 3
    [Hl,gl] = geom.project_polyhed2dims(H,g,dims_lower);
    Pl = Polyhedron('A',Hl,'b',gl);
    Pl.plot('color',[0,0.7,0.3],'alpha',0.3);    
    axis equal
    ax = gca;
    ax.DataAspectRatioMode = 'manual';
    ax.DataAspectRatio = [1,1,1];
    ax.PlotBoxAspectRatioMode = 'manual';
    ax.PlotBoxAspectRatio = [1,1,1];    
end
