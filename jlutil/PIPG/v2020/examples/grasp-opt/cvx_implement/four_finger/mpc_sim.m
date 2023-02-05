clear all
close all
clc

problem_data

cvx_solver ecos

NN = Nbar-N+1;

% contrainer for state and input at each MPC instant
X = zeros(nx,NN);
U = zeros(nu,NN);

fprintf("MPC simulation\n--------------\n")
for k=1:NN
    
    % record state 
    X(:,k) = x0; 
    
    % set reference
    yy = yybar(:,k:k+N-1);
    y = reshape(yy,[nx*N,1]);
    
    % parse + solve
    cvx_begin quiet
        variables x(nx*N) u(nu*(N-1))
        minimize  q_wt*sum_square(x(1:(N-1)*nx)-y(1:(N-1)*nx)) + qf_wt*sum_square(x((N-1)*nx+1:N*nx)-y((N-1)*nx+1:N*nx)) + r_wt*sum_square(u-uref)
        subject to
            for j=1:N-1
                x(j*nx+1:(j+1)*nx) == Ad*x((j-1)*nx+1:j*nx) + Bd*u((j-1)*nu+1:j*nu) + gd; % dynamics constraint
                Gam*u((j-1)*nu+1:j*nu) == 0;                                              % rotational equilibrium
                norm(u((j-1)*nu+2:(j-1)*nu+3)) <= mu1*u((j-1)*nu+1);                      % friction cone
                norm(u((j-1)*nu+5:(j-1)*nu+6)) <= -mu2*u((j-1)*nu+4);                     % friction cone            
                norm(u([(j-1)*nu+7,(j-1)*nu+9])) <= mu3*u((j-1)*nu+8);                    % friction cone            
                norm(u([(j-1)*nu+10,(j-1)*nu+12])) <= -mu4*u((j-1)*nu+11);                % friction cone   
                norm(u((j-1)*nu+1:(j-1)*nu+3)) <= F1max;                                  % grasp force magnitude
                norm(u((j-1)*nu+4:(j-1)*nu+6)) <= F2max;                                  % grasp force magnitude
                norm(u((j-1)*nu+7:(j-1)*nu+9)) <= F3max;                                  % grasp force magnitude
                norm(u((j-1)*nu+10:(j-1)*nu+12)) <= F4max;                                % grasp force magnitude
                norm(x(j*nx+4:j*nx+6)) <= Vmax;                                           % velocity magnitude                  
            end
            x(1:nx) == x0; 
    cvx_end
    
    fprintf("\nInstant %2d | CVX status: %s | Solve time: %.3f s",k,cvx_status,cvx_cputime);
    if ~strcmp(cvx_status,'Solved')
        fprintf("\n")
        error("Problem is not solved properly.")
    end
    
    % propgate system with first control input; disturbance is added
    x0 = Ad*x0 + Bd*u(1:nu) + gd + 0.001*randn(nx,1);
    
    % record input
    U(:,k) = u(1:nu);
    
    % warm-start is not possible
    
end
fprintf("\n")

plot_soln(X,U,yybar,Nbar-N,s1,s2,s3,s4,a)

