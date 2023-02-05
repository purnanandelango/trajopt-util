close all
clear all
clc

problem_data
cvx_solver gurobi

I6 = eye(6);

cvx_begin
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

% activity of constraints
tlt_cone = ones(1,N-1);    % tilt-angle cone
gs_cone = ones(1,N);       % glide-slope cone
thrst_ubnd = ones(1,N-1);  % thrust upper bound
thrst_lbnd = ones(1,N-1);  % thrust lower bound
V_ubnd = ones(1,N);        % speed upper bound
for j = 1:N
    gs_cone(j) = norm(x(1:2,j)-y(1:2,N)) - cvec*(x(1:3,j)-y(1:3,N));
    V_ubnd(j) = norm(x(4:6,j)) - Vmax;
    if j<N
        tlt_cone(j) = norm(u(1:2,j)) - dvec*u(:,j);
        thrst_ubnd(j) = norm(u(:,j)) - umax;
        thrst_lbnd(j) = umin - u(3,j);
    end
end

fprintf("\nSolution performance\n--------------------\nFinal position deviation: %.3f m\nFinal velocity deviation: %.3f m s^-1\n",norm(x(1:3,end)-y(1:3,end)),norm(x(4:6,end)-y(4:6,end)));
fprintf("\nConstraint activity\n-------------------\nGlide-slope cone: %.3f\nThrust tilt cone: %.3f\nThrust magnitude upper bound: %.3f\nThrust lower bound approx.: %.3f\nVelocity magnitude upper bound: %.3f\n",max(gs_cone),max(tlt_cone),max(thrst_ubnd),max(thrst_lbnd),max(V_ubnd));

plot_soln(x,u,N,tvec,gamgs,thetmax,Vmax,umax,umin)

