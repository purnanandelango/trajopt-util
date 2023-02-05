clear all
close all
clc

problem_data

cvx_solver mosek

NN = Nbar-N+1;

% contrainer for state and input at each MPC instant
X = zeros(nx,NN);
U = zeros(nu,NN);

fprintf("MPC simulation\n--------------\n")
for k=1:NN
    
    % record state 
    X(:,k) = x0; 
    
    % set reference
    y = ybar(:,k:k+N-1);
    
    % parse + solve
    cvx_begin quiet
        variable x(6,N)
        variable u(3,N-1)
        expression obj_expr(N-1)
        for i = 1:N-2
            obj_expr(i) = 0.5*quad_form(x(:,i+1)-y(:,i+1),Q) + 0.5*quad_form(u(:,i),R);
        end
        obj_expr(N-1) = 0.5*quad_form(x(:,N)-y(:,N),Qf) + 0.5*quad_form(u(:,N-1),R);
        minimize sum(obj_expr)
        subject to
            x(:,1) == x0
            for i = 1:N-1
                x(:,i+1) == Ad*x(:,i) + Bd*u(:,i) + gd;
                norm(u(:,i)) <= umax;
                norm(x(4:6,i+1)) <= Vmax;
                norm(x(1:2,i+1)-y(1:2,N)) - cvec*(x(1:3,i+1)-y(1:3,N)) <= 0;
                norm(u(1:2,i)) - dvec*u(:,i) <= 0;
                u(3,i) >= umin;
            end
    cvx_end
    
    fprintf("\nInstant %2d | CVX status: %s | Solve time: %.3f s",k,cvx_status,cvx_cputime);
    if ~(strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved'))
        fprintf("\n")
        error("Problem is not solved properly.")
    end
    
    % propgate system with first control input; disturbance is added
    x0 = Ad*x0 + Bd*u(1:nu,1) + gd + 0.1*randn(nx,1);
    
    % record input
    U(:,k) = u(1:nu,1);
    
    % warm-start is not possible
    
end
fprintf("\n")

tmpc = 0:Del:((NN-1)*Del);
plot_soln(X,U,NN,tmpc,gamgs,thetmax,Vmax,umax,umin)