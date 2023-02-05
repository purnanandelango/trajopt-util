%% MPC example 1 for testing FiOrdOs
% input constrained MPC
% no constraints on state

clear all
clc

%% problem definition

A=[0.8 1; 0 0.9];
B=[-1;2];
[nx,nu]=size(B);
 
N=10;
Q=eye(nx);
R=eye(nu);
P=eye(nx);


umin=-7;
umax= 7;

% uref = (0:(N-1))'/10;
uref = randn(nu*N,1);

noise_std = 0.001;

% P=dlyap(A',Q); % terminal weights ensuring closed loop stability
 
UU=SimpleSet(N);
UU.addSet(1:N,EssBox(nu,'l',umin,'u',umax));
 
Ai=A;
AA=Ai;
for i=2:N
    Ai=A*Ai;
    AA=[AA;Ai];
end
 
AiB=B;
BB=kron(eye(N),AiB);
for i=1:N-1
    AiB=A*AiB;
    BB=BB+kron(diag(ones(N-i,1),-i),AiB);
end
 
QQ=blkdiag(kron(eye(N-1),Q),P);
RR=kron(eye(N),R);
H =(BB'*QQ*BB + RR);
 
op=OptProb('H',H, 'g','param', 'c','param', 'X',UU);

%% solver definition and code generation

s=Solver(op,'approach','primal','algo','fgm');
s.setSettings('approach','inlineH',false); % write the matrix-vector multiplication involving $H$ element by element in C.
s.setPrecond('algo',eye(size(H)));
s.setSettings('algo', 'stopg',true, 'stopgEps', 1e-3); % stop the algorithm when the function value is epsilon-suboptimal
s.setSettings('algo', 'maxit',5000) % stop the algorithm when max iterations is reached
s.generateCode('prefix','ex1_','forceOverwrite',true);
rehash;
ex1_mex_make;

%%
Nmpc=30;
 
trajU=nan(nu,Nmpc);
trajX=nan(nx,Nmpc+1);
solve_time(Nmpc,1) = 0;
rel_err(Nmpc,1) = 0;

noise_val = noise_std*randn(nx,Nmpc);

trajX(:,1)=[1;2];  % initial state
 
mparams=struct();
msetgs=struct();
% msetgs.algo.maxit=5000;

cvx_solver_name = 'ecos';
cvx_solver(cvx_solver_name);

clc
fprintf('                            ..:: MPC simulation ::..                                     \n');
fprintf('|-------------|-----------------------------------------------|-------------------------|\n');
fprintf(horzcat('|  CVX Status |  Relative error between FiOrdOs and cvx-',cvx_solver_name,'  | FiOrdOs solve time (ms) |\n'));
fprintf('|-------------|-----------------------------------------------|-------------------------|\n');
for k=1:Nmpc
    x0=trajX(:,k);
 
    % parametric data
    mparams.g=BB'*QQ*AA*x0 - RR*uref;
    mparams.c=0.5*x0'*(AA'*QQ*AA + Q)*x0;
    
    % warmstarting
    if k>1
        % msetgs.algo.init=mres.x;                         % unshifted
        msetgs.algo.init=[mres.x(nu+1:end);zeros(nu,1)];   % shifted
    end
 
    % call solver
    tic
    mres=ex1_mex(mparams,msetgs);
    solve_time(k) = 1000*toc;

    % solve problem via cvx
    cvx_begin quiet
        variables x_cvx(nx,N+1) u_cvx(nu,N)
        expression obj_val(N+1)
        for j=1:N
            x_cvx(:,j+1) == A*x_cvx(:,j) + B*u_cvx(:,j);
            umin <= u_cvx(:,j) <= umax;
            obj_val(j) = quad_form(x_cvx(:,j),Q) + quad_form(u_cvx(:,j)-uref((j-1)*nu+1:j*nu),R); 
            % obj_val(j) = quad_form(x(:,j)-x_ref((j-1)*nx+1:j*nx),Q) + quad_form(u(:,j)-u_ref((j-1)*nu+1:j*nu),R); 
        end
        x_cvx(:,1) == x0;
        % obj_val(N+1) = quad_form(x(:,N+1)-x_ref(N*nx+1:(N+1)*nx),Q);
        obj_val(N+1) = quad_form(x_cvx(:,N+1),Q);
        minimize sum(obj_val)
    cvx_end
    if ~strcmp('Solved',cvx_status)
        cvx_status = 'Unslvd';
    end
    u_cvx = reshape(u_cvx,[nu*N,1]);
    
    rel_err(k) = norm(u_cvx - mres.x,Inf)/norm(u_cvx,Inf);
    
    if isnan(rel_err(k))
        fprintf('|    %s   |                      %.3e                      |         %.3e       |\n',cvx_status,rel_err(k),solve_time(k));    
    else
        fprintf('|    %s   |                   %.3e                   |         %.3e       |\n',cvx_status,rel_err(k),solve_time(k));        
    end
    
    % apply input u_0 to system
    trajU(:,k) = mres.x(1:nu);
    trajX(:,k+1) = A*trajX(:,k) + B*trajU(:,k) +noise_val(:,k);
end
 
fprintf('\n\nAverage FiOrdOs solve time: %.3f ms\n',mean(solve_time(2:end)));

% plot results
% figure(1); clf; stairs(0:Nsteps-1,trajU'); title('inputs'); xlabel('steps');
% figure(2); clf; plot(trajX(1,:),trajX(2,:),'*:'); title('state trajectory'); xlabel('$x_1$'); ylabel('$x_2$');
