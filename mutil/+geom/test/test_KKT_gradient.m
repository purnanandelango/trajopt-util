% Test the estimate of gradient of the Lagrangian formed for the (nonconvex) ellipsoid projection problem

clearvars
close all
clc

% rng default
        
n = 6;        

fprintf("\n");
for k = 1:100
        
        A = randn(n,n);
        P = A*A';
        invP = inv(P);
        x = randn(n,1);

        In = eye(n);
        y = @(lam) (In - lam*invP)\x;
        func_g = @(lam) y(lam)'*invP*y(lam) - 1;
        grad_func_g = @(lam) 2*y(lam)'*invP'*(inv(In - lam*invP)')*invP*y(lam);

        N = 1000;
        lam_grid = linspace(-10,10,N);
        dlam = lam_grid(2) - lam_grid(1);

        true_grad_g(N-1) = 0;
        num_grad_g(N-1) = 0;
        g_val(N) = func_g(lam_grid(N));
        for j = 1:N-1
                num_grad_g(j) = ( func_g(lam_grid(j+1)) - func_g(lam_grid(j)) )/dlam;
                true_grad_g(j) = grad_func_g(lam_grid(j));
                g_val(j) = func_g(lam_grid(j));
        end

        opts = optimoptions('fsolve','Algorithm','levenberg-marquardt','FunctionTolerance',1e-9,'StepTolerance',1e-8,...
                            'Display','iter','SpecifyObjectiveGradient',true);
        [lam_str,~,exit_flag] = fsolve(@(lam) fun_solve(lam,func_g,grad_func_g),4,opts);

        fprintf("Dual variable = %.3f, Constraint violation = %.2e\n",lam_str,func_g(lam_str));

end

figure
plot(lam_grid(1:end-1),num_grad_g,'-b','DisplayName','Forward diff. derivative');
hold on
plot(lam_grid(1:N-1),true_grad_g,'--r','DisplayName','Exact derivative');
plot(lam_grid,g_val,'-k','DisplayName','Function value');
ylim([-10,10])
legend('AutoUpdate','on');

function [g_val,dg_val] = fun_solve(lam,func_g,grad_func_g)
    g_val = func_g(lam);
    dg_val = grad_func_g(lam);
end
