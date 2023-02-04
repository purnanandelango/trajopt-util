%% Oscillating masses example via FiOrdOs
% state and input constrained MPC

clear all
clc

%% problem definition

nx = 20;
nu = 10;

A = [zeros(nx/2) eye(nx/2);
      (diag(ones(nx/2-1,1),-1)+diag(-2*ones(nx/2,1))+diag(ones(nx/2-1,1),1)) zeros(nx/2)];
B = [zeros(nx/2);
     eye(nx/2)];
 
N=10;
Q=eye(nx);
R=eye(nu);
P=eye(nx);

umin=-2*ones(nu,1);
umax= -umin;
xmin = -[2*ones(nx/2,1);0.8*ones(nx/2,1)];
xmax = -xmin;

% u_ref = (0:(N-1))'/10;
% x_ref(nx*(N+1),1) = 0;
% for j=1:N+1
%     x_ref((j-1)*nx+1:nx*j) = [j;2*j]/10;
% end

u_ref = randn(nu*N,1);
x_ref = repmat([ones(nx/2,1);zeros(nx/2,1)],[1+N,1]);

z_ref = [u_ref;x_ref(nx+1:end)];

noise_std = 0.001;
 
ZZ=SimpleSet(2*N);
ZZ.addSet(1:N,EssBox(nu,'l',umin,'u',umax));
ZZ.addSet(N+1:2*N,EssBox(nx,'l',xmin,'u',xmax));
 
H =blkdiag(kron(eye(N),R),  kron(eye(N-1),Q),eye(size(P)));
 
Aeu=kron(eye(N),-B);
Aex=blkdiag(eye((N-1)*nx),P^(-1/2)) - kron(diag(ones(N-1,1),-1),A);
Ae =[Aeu,Aex];
 
op=OptProb('H',H, 'g','param', 'c','param', 'X',ZZ, 'Ae',Ae, 'be','param');     
 

%% solver definition and code generation

s=Solver(op,'approach','dual','algoOuter','fgm','algoInner','fgm');
% maximum number of iterations
s.setSettings('algoOuter', 'maxit',5000);
s.setSettings('algoInner', 'maxit',5000);
% gradient-map stopping criterion
s.setSettings('algoOuter', 'stopg',true, 'stopgEps',1e-3); 
s.setSettings('algoInner', 'stopg',true, 'stopgEps',1e-4);
s.generateCode('prefix','ex12_','forceOverwrite',true);
rehash;
ex12_mex_make;

%%
Nmpc=30;
 
trajU=nan(nu,Nmpc);
trajX=nan(nx,Nmpc+1);
solve_time(Nmpc,1) = 0;
rel_err(Nmpc,1) = 0;

noise_val = noise_std*randn(nx,Nmpc);

trajX(:,1)=1.5*[ones(nx/2,1);zeros(nx/2,1)]; % initial state
 
mparams=struct();
msetgs=struct();

cvx_solver_name = 'mosek';
cvx_solver(cvx_solver_name)

clc
fprintf('                            ..:: MPC simulation ::..                                     \n');
fprintf('|-------------|-----------------------------------------------|-------------------------|\n');
fprintf(horzcat('|  CVX Status |  Relative error between FiOrdOs and cvx-',cvx_solver_name,' | FiOrdOs solve time (ms) |\n'));
fprintf('|-------------|-----------------------------------------------|-------------------------|\n');
for k=1:Nmpc
    x0=trajX(:,k);
 
    % parametric data
    mparams.g=-H*z_ref;
    mparams.c=0.5*x0'*Q*x0;
    mparams.be =[A;zeros((N-1)*nx,nx)]*x0;
    
    if k>1 % warmstarting
        % % if approach is primal-dual
        % msetgs.approach.apprInitX  = mres.x;                  
        % msetgs.approach.apprInitLa = mres.la;
        
        % % if approach is dual
        msetgs.algoInner.init=[mres.x(nu+1:N*nu);zeros(nu,1);mres.x(N*nu+nx+1:end);zeros(nx,1)];   % shifted
        msetgs.algoOuter.init=mres.la;                                                             % unshifted
    end
 
    % call solver
    tic
    mres=ex12_mex(mparams,msetgs);
    solve_time(k) = 1000*toc;
    
    % solve problem via cvx
    cvx_begin quiet
        variables x_cvx(nx,N+1) u_cvx(nu,N)
        expression obj_val(N+1)
        for j=1:N
            x_cvx(:,j+1) == A*x_cvx(:,j) + B*u_cvx(:,j);
            umin <= u_cvx(:,j) <= umax;
            xmin <= x_cvx(:,j) <= xmax;
            obj_val(j) = quad_form(x_cvx(:,j)-x_ref((j-1)*nx+1:nx*j),Q) + quad_form(u_cvx(:,j)-u_ref((j-1)*nu+1:j*nu),R); 
        end
        x_cvx(:,1) == x0;
        xmin <= x_cvx(:,N+1) <= xmax;
        obj_val(N+1) = quad_form(x_cvx(:,N+1)-x_ref(N*nx+1:(N+1)*nx),Q);
        minimize sum(obj_val)
    cvx_end
    if ~strcmp('Solved',cvx_status)
        cvx_status = 'Unslvd';
    end
    u_cvx = reshape(u_cvx,[nu*N,1]);
    x_cvx = reshape(x_cvx(nx+1:end),[nx*N,1]);
    
    rel_err(k) = norm([u_cvx;x_cvx] - mres.x,Inf)/norm([u_cvx;x_cvx],Inf);

    if isnan(rel_err(k))
        fprintf('|    %s   |                      %.3e                      |         %.3e       |\n',cvx_status,rel_err(k),solve_time(k));    
    else
        fprintf('|    %s   |                   %.3e                   |         %.3e       |\n',cvx_status,rel_err(k),solve_time(k));        
    end
    
    % apply input u_0 to system
    trajU(:,k) = mres.x(1:nu);
    trajX(:,k+1) = A*trajX(:,k) + B*trajU(:,k) + noise_val(:,k);
end

fprintf('\n\nAverage FiOrdOs solve time: %.3f     ms\nAverage relative error:     %.3e\n',mean(solve_time(2:end)),mean(rel_err));
