clc
clearvars
close all

nx = 10;
nu = nx/2;
N = 50;

tic
P_2D = parse_via_optimizer_2D(nx,nu,N);
fprintf("Parse time by optimizer with 2D arrays: %.3f ms\n",toc*1000);

tic
P_3D = parse_via_optimizer_3D(nx,nu,N);
fprintf("Parse time by optimizer with 3D arrays: %.3f ms\n",toc*1000);

M = 5;
% Store execution times
s1 = zeros(1,M);
s2 = zeros(1,M);
s3 = zeros(1,M);
% Store error in optimal values
err_obj_12 = zeros(1,M);
err_obj_13 = zeros(1,M);
for k = 1:M
    % Get random problem data
    [A,B,Q,R,xmax,umax,x0] = generate_random_data(nx,nu,N);  

    tic
    objval = solve_via_optimize(nx,nu,N,A,B,Q,R,xmax,umax,x0);
    s1(k) = toc*1000;

    sqrtQ = sqrt(Q);
    sqrtR = sqrt(R);    

    tic
    sol2 = P_2D(A,B,sqrtQ,sqrtR,xmax,umax,x0);
    s2(k) = toc*1000;

    % Convert parameter matrices from 2D to 3D
    [A,B,Q,R] = reshape_data(nx,nu,N,A,B,Q,R);
    sqrtQ = sqrt(Q);
    sqrtR = sqrt(R);

    tic
    sol3 = P_3D(A,B,sqrtQ,sqrtR,xmax,umax,x0);
    s3(k) = toc*1000;    

    err_obj_12(k) = norm(objval-sol2{1})/sol2{1};
    err_obj_13(k) = norm(objval-sol3{1})/sol3{1};
end

figure
subplot(1,2,1)
plot(1:M,s1,'ob','DisplayName','\texttt{optimize}');
hold on
plot(1:M,s2,'or','DisplayName','\texttt{optimizer} 2D');
plot(1:M,s3,'og','DisplayName','\texttt{optimizer} 3D');
legend('Location','best');
ylabel('ms')
xlabel('Samples')
xlim([1,M])

subplot(1,2,2)
semilogy(1:M,err_obj_12*100,'ob');
hold on
semilogy(1:M,err_obj_13*100,'+r');
title('\% Error in optimal value','FontSize',16);
xlabel('Samples')
xlim([1,M]);

function [A,B,Q,R,xmax,umax,x0] = generate_random_data(nx,nu,N)
    A = randn(nx^2,N-1);
    B = randn(nx*nu,N-1);
    Q = zeros(nx^2,N);
    R = zeros(nu^2,N-1);
    for k = 1:N-1
        A(:,k) = reshape(A(:,k)/norm(A(:,k)),[nx^2,1]);
        Q(:,k) = reshape(diag(randn(nx,1) .^2),[nx^2,1]);
        R(:,k) = reshape(diag(randn(nu,1) .^2),[nu^2,1]);
    end
    Q(:,N) = reshape(diag(randn(nx,1) .^2),[nx^2,1]);
    umax = 0.5; 
    xmax = 1;    
    x0 = 0.1*ones(nx,1);    
end

% Convert parameters matrices from 2D to 3D
function [A2,B2,Q2,R2] = reshape_data(nx,nu,N,A,B,Q,R)
    A2 = zeros(nx,nx,N-1);
    B2 = zeros(nx,nu,N-1);
    Q2 = zeros(nx,nx,N);
    R2 = zeros(nu,nu,N-1);
    for k = 1:N-1
        A2(:,:,k) = reshape(A(:,k),[nx,nx]);
        B2(:,:,k) = reshape(B(:,k),[nx,nu]);
        Q2(:,:,k) = reshape(Q(:,k),[nx,nx]);
        R2(:,:,k) = reshape(R(:,k),[nu,nu]);
    end
    Q2(:,:,N) = reshape(Q(:,N),[nx,nx]);
end

function [objval,x,u] = solve_via_optimize(nx,nu,N,A,B,Q,R,xmax,umax,x0)
    % Caution: this command will mess up the optimizer object definition.
    % Refrain from using.
    % yalmip clear
    x = sdpvar(nx,N);
    u = sdpvar(nu,N-1);
    cnstr = x(:,1) == x0;
    objval = 0;
    for k = 1:N-1
        cnstr = [cnstr;
                 x(:,k+1) == reshape(A(:,k),[nx,nx])*x(:,k) + reshape(B(:,k),[nx,nu])*u(:,k);
                 norm(x(:,k),2) <= xmax;
                 norm(u(:,k),inf) <= umax];
        objval = objval + x(:,k)'*reshape(Q(:,k),[nx,nx])*x(:,k) + u(:,k)'*reshape(R(:,k),[nu,nu])*u(:,k);
    end
    objval = objval + x(:,N)'*reshape(Q(:,N),[nx,nx])*x(:,N);
    optimize(cnstr,objval,sdpsettings('solver','ecos','verbose',0));
    objval = value(objval);
    x = value(x);
    u = value(u);
end

% Uses 2D arrays to store matrices
function P = parse_via_optimizer_2D(nx,nu,N)
    xvar = sdpvar(nx,N);
    uvar = sdpvar(nu,N-1);
    Avar = sdpvar(nx^2,N-1);
    Bvar = sdpvar(nx*nu,N-1);
    Qvar = sdpvar(nx^2,N);
    q = sdpvar(1,N);
    Rvar = sdpvar(nu^2,N-1);
    r = sdpvar(1,N-1);
    xmaxvar = sdpvar();
    umaxvar = sdpvar();
    x0var = sdpvar(nx,1);
    cnstr = xvar(:,1) == x0var;
    objval = 0;
    for k = 1:N-1
        cnstr = [cnstr;
                 xvar(:,k+1) == reshape(Avar(:,k),[nx,nx])*xvar(:,k) + reshape(Bvar(:,k),[nx,nu])*uvar(:,k);
                 norm(xvar(:,k),2) <= xmaxvar;
                 norm(uvar(:,k),inf) <= umaxvar;
                 norm(reshape(Qvar(:,k),[nx,nx])*xvar(:,k)) <= q(k);
                 q(k) >= 0;
                 norm(reshape(Rvar(:,k),[nu,nu])*uvar(:,k)) <= r(k);
                 r(k) >= 0];
        objval = objval + q(k)^2 + r(k)^2;  % Note that while using optimizer we need to reformulate the quadratic objective as cones
    end
    cnstr = [cnstr;norm(reshape(Qvar(:,N),[nx,nx])*xvar(:,N)) <= q(N)];
    objval = objval + q(N)^2;
    P = optimizer(cnstr,objval,sdpsettings('solver','ecos'),{Avar,Bvar,Qvar,Rvar,xmaxvar,umaxvar,x0var},{objval,xvar,uvar});
end

% Uses 3D arrays to store matrices
function P = parse_via_optimizer_3D(nx,nu,N)
    xvar = sdpvar(nx,N,'full');
    uvar = sdpvar(nu,N-1,'full');
    Avar = sdpvar(nx,nx,N-1,'full');
    Bvar = sdpvar(nx,nu,N-1,'full');
    Qvar = sdpvar(nx,nx,N,'diagonal');
    q = sdpvar(1,N);
    Rvar = sdpvar(nu,nu,N-1,'diagonal');
    r = sdpvar(1,N-1);
    xmaxvar = sdpvar();
    umaxvar = sdpvar();
    x0var = sdpvar(nx,1);
    cnstr = xvar(:,1) == x0var;
    objval = 0;
    for k = 1:N-1
        cnstr = [cnstr;
                 xvar(:,k+1) == Avar(:,:,k)*xvar(:,k) + Bvar(:,:,k)*uvar(:,k);
                 norm(xvar(:,k),2) <= xmaxvar;
                 norm(uvar(:,k),inf) <= umaxvar;
                 norm(Qvar(:,:,k)*xvar(:,k)) <= q(k);
                 q(k) >= 0;
                 norm(Rvar(:,:,k)*uvar(:,k)) <= r(k);
                 r(k) >= 0];
        objval = objval + q(k)^2 + r(k)^2;  % Note that while using optimizer we need to reformulate the quadratic objective as cones
    end
    cnstr = [cnstr; norm(Qvar(:,:,N)*xvar(:,N)) <= q(N)];
    objval = objval + q(N)^2;
    P = optimizer(cnstr,objval,sdpsettings('solver','ecos'),{Avar,Bvar,Qvar,Rvar,xmaxvar,umaxvar,x0var},{objval,xvar,uvar});
end