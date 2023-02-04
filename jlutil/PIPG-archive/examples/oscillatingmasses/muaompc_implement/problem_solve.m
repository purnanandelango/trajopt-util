%% Oscillating masses - solve once
% nx = 20, nu = 10, in_iter ~ 155, ex_iter ~ 160 gives rd2o ~ 1e-3

ctl = mpc_ctl;  % create an instance of the class
ctl.x_ref = x_ref;
ctl.u_ref = u_ref;
ctl.conf.in_iter = 155;  % configure the algorithm
ctl.conf.ex_iter = 160;  % configure the algorithm
ctl.conf.warmstart = false;

x0 = 1.5*[ones(nx/2,1);zeros(nx/2,1)];  % current state

% forming a QP is only necessary if a different QP solver is used
ctl.form_qp(x0)                 % form the QP for the current state
qpx = ctl.qpx;                  % qpx contains the QP in standard form
u1 = quadprog(qpx.HoL, qpx.gxoL, [], [], [], [], qpx.u_lb, qpx.u_ub);

% the QP is automatically formed when using its own algorithm (ALM+FGM) 
ctl.solve_problem(x0);          % solve the MPC problem for current state

% solve problem via cvx
cvx_solver mosek
cvx_begin
    variables x(nx,N+1) u2(nu,N)
    expression obj_val(N+1)
    for j=1:N
        x(:,j+1) == Ad*x(:,j) + Bd*u2(:,j);

        u_lb <= u2(:,j) <= u_ub;
        e_lb <= x(:,j) <= e_ub;

        obj_val(j) = quad_form(x(:,j)-x_ref((j-1)*nx+1:j*nx),Q) + quad_form(u2(:,j)-u_ref((j-1)*nu+1:j*nu),R); 
    end
    x(:,1) == x0;
    f_lb <= x(:,N+1) <= f_ub;
    obj_val(N+1) = quad_form(x(:,N+1)-x_ref(N*nx+1:(N+1)*nx),Q);
    minimize sum(obj_val)
cvx_end
u2 = reshape(u2,[nu*N,1]);

% disp('Input computed by muaompc');
% disp(ctl.u_opt')  % display the computed input sequence
% ctl.u_opt is the MPC approximation of u
disp('Difference between input computed by muaompc and optimal input from quadprog');
disp(norm(u1 - ctl.u_opt,Inf)/norm(u1,Inf))  
disp('Difference between input computed by muaompc and optimal input from cvx');
disp(norm(u2 - ctl.u_opt,Inf)/norm(u2,Inf))  
