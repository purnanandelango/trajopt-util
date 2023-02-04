%% Oscillating masses - solve recursively

clc

ctl = mpc_ctl;  % create an instance of the class
ctl.conf.in_iter = 155;  % configure the algorithm
ctl.conf.ex_iter = 160;
ctl.conf.warmstart = true;

x0 = 1.5*[ones(nx/2,1);zeros(nx/2,1)];  % current state

noise_val = noise_std*randn(nx,Nmpc);
Umpc = zeros(nu,Nmpc);
solve_time(Nmpc) = 0;

cvx_solver mosek

for j=1:Nmpc
    ctl.x_ref = x_ref((j-1)*nx+1:(j-1)*nx+(N+1)*nx);
    ctl.u_ref = u_ref((j-1)*nu+1:(j-1)*nu+N*nu);
    tic
    ctl.solve_problem(x0);  % solve the QP via muaompc
    solve_time(j) = toc*1000;
    
    % solve QP via cvx - mosek 
    cvx_begin quiet
        variables x(nx,N+1) u(nu,N)
        expression obj_val(N+1)
        for k=1:N
            x(:,k+1) == Ad*x(:,k) + Bd*u(:,k);

            u_lb <= u(:,k) <= u_ub;
            e_lb <= x(:,k) <= e_ub;

            obj_val(k) = quad_form(x(:,k)-x_ref((j-1)*nx+(k-1)*nx+1:(j-1)*nx+k*nx),Q) + quad_form(u(:,k)-u_ref((j-1)*nu+(k-1)*nu+1:(j-1)*nu+k*nu),R); 
        end
        x(:,1) == x0;
        f_lb <= x(:,N+1) <= f_ub;
        obj_val(N+1) = quad_form(x(:,N+1)-x_ref((j-1)*nx+N*nx+1:(j-1)*nx+(N+1)*nx),Q);
        minimize sum(obj_val)
    cvx_end
    u = reshape(u,[nu*N,1]);
    
    x0 = Ad*x0 + Bd*ctl.u_opt(1:nu) + noise_val(:,j);

    assert(abs(norm(ctl.u_opt-u,Inf)/norm(u,Inf))<1e-3);
    fprintf('\nInstant %2d | Relative d2o : %s | Solve Time : %.3f ms',j,abs(norm(ctl.u_opt-u,Inf)/norm(u,Inf)),solve_time(j));
    
    Umpc(:,j) = ctl.u_opt(1:nu);
end

fprintf('\n')