% Test estimation of signed distance to ellipsoid by solving an NLP via fmincon

clearvars
clc
close all

n = 6;
Q = qr(randn(n));
D = diag(randn(1,n) .^ 2);
A = Q*D;
invP = inv( A*A' );

% Store samples
xs = [];
ys = [];
y2s = [];
y3s = [];

for k = 1:25

    fprintf("\nSample %d\n",k);
    
    x = randn(n,1);
    xs = [xs,x];
    
    opts = optimoptions('fmincon','Algorithm','interior-point',...
                        'SpecifyConstraintGradient',true,'SpecifyObjectiveGradient',true,...
                        'CheckGradients',false,...
                        'MaxIterations',1000,...
                        'HessianFcn',@(y,lambda) hessianfunc(y,lambda,invP),'Display','none');
    [y,obj_val,exit_flag] = fmincon(@(y) obj_fun(y,x),x+0.1*rand(n,1),[],[],[],[],[],[],@(y) nonlin_con(y,invP),opts);

    [y2,obj_val2] = geom.sign_dist_ellip_solveNLP(x,A);
    [y3,obj_val3] = geom.sign_dist_ellip_solveKKT(x,A);

    assert(ismember(exit_flag,[1,2]))
    
    [~,cnstr_viol] = nonlin_con(y,invP);
    percent_error2 = 100*norm(y-y2)/norm(y2);
    percent_error3 = 100*norm(y-y3)/norm(y3);    

    fprintf("Constraint violation: %9.2e, Percent error (NLP): %9.2e, Percent error (KKT): %9.2e\n",cnstr_viol,percent_error2,percent_error3);

    ys = [ys, y];
    y2s = [y2s, y2];    
    y3s = [y3s, y3];        
end


% Functions required for solving the NLP by calling fmincon

function [obj_val,obj_grad] = obj_fun(y,x)
    obj_val = norm(y-x)^2;
    obj_grad = 2*(y-x);
end

function [c,ceq,dc,dceq] = nonlin_con(y,invP)
    ceq = y'*invP*y - 1;
    c = [];
    dc = [];
    dceq = 2*invP*y;
end

function H = hessianfunc(y,lambda,invP)
    n = length(y);
    H = 2*(eye(n) + lambda.eqnonlin*invP);
end
