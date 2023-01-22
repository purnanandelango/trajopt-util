Test estimation of signed distance to ellipsoid by solving an NLP via fmincon

clearvars
clc
close all

n = 6;

for k = 1:100

    fprintf("\nSample %d\n",k);

    Q = randn(n);
    Q = qr(Q);
    D = diag(randn(1,n) .^ 2);
    A = Q*D;
    
    I = eye(n);
    invP = (A\I)'/A;
    
    x = randn(n,1);
    
    opts = optimoptions('fmincon','Algorithm','interior-point',...
                        'SpecifyConstraintGradient',true,'SpecifyObjectiveGradient',true,...
                        'CheckGradients',false,...
                        'MaxIterations',1000,...
                        'HessianFcn',@(y,lambda) hessianfunc(y,lambda,invP),'Display','iter-detailed');
    [y,obj_val,exit_flag] = fmincon(@(y) obj_fun(y,x),x+0.1*rand(n,1),[],[],[],[],[],[],@(y) nonlin_con(y,invP),opts);

    assert(ismember(exit_flag,[1,2]))

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
