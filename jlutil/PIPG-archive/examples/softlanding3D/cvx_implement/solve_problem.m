close all
clear all
clc

problem_data;

cvx_solver mosek;

cvx_begin
    variables x(nx,N) f(nu,N-1)

    for j=1:N-1
        x(:,j+1) == Ad*x(:,j) + Bd*f(:,j) + gd;
        ( -mu1(j+1)*x(7,j+1) - f(4,j) )/norm([-mu1(j+1),-1]) <= -mu1(j+1)*(1+z0(j+1))/norm([-mu1(j+1),-1]);
        (  mu2(j+1)*x(7,j+1) + f(4,j) )/norm([ mu2(j+1), 1]) <=  mu2(j+1)*(1+z0(j+1))/norm([ mu2(j+1), 1]);
        norm(f(1:3,j)) <= f(4,j);
        norm(x(4:6,j+1)) <= Vmax;
        norm(E*x(1:3,j+1)) <= nhat'*x(1:3,j+1)*tangamgs;
        z0(j+1) <= x(7,j+1) <= z1(j+1);
    end

%     for j=1:N-1
%         
%         x(:,j+1) == Ad*x(:,j) + Bd*f(:,j) + gd;
%         norm(f(1:3,j)) <= f(4,j);
%         norm(x(4:6,j+1)) <= Vmax;
%         norm(E*x(1:3,j+1)) <= nhat'*x(1:3,j+1)*tangamgs;
%         z0(j+1) <= x(7,j+1) <= z1(j+1);
%         
%         if j<=N-2
%            -mu1(j+1)*x(7,j+1)-f(4,j+1) <= -mu1(j+1)*(1+z0(j+1));
%             mu2(j+1)*x(7,j+1)+f(4,j+1) <= mu2(j+1)*(1+z0(j+1));        
%         end
%         
%     end
%     mu1(1)*(1+z0(1)) - mu1(1)*x0(7) <= f(4,1) <= mu2(1)*(1+z0(1)) - mu2(1)*x0(7);

    x(:,1) == x0;
%     x(1:6,N) == [rf;vf];
    
    
%     minimize sum(f(4,:))
    expression obj_val(N-1)
    for j=1:N-1
        obj_val(j) = quad_form(x(:,j+1)-yref(:,j+1),Q) + quad_form(f(:,j),R);
    end
    minimize sum(obj_val)
cvx_end

m = exp(x(7,:));
T = zeros(3,N-1);
for j=1:N-1
    T(:,j) = f(1:3,j) * m(j);
end

lcvx_test(N-1) = 0;
for j=1:N-1
    lcvx_test(j) = m(j)*f(4,j);
end

plot_soln
