function [y,dist_val] = sign_dist_ellip_solveNLP(x,A,varargin)
% Compute projection of x onto ellipsoid defined by { z | norm(A\z) <= 1 }
% The ellipsoid can also be defined by { z | z'*inv(P)*z <= 1 }, where
%         P = inv( inv(A)'*inv(A) ) = A*A'
% The boundary of the ellipsoid is represented by y = A*u, where 
% u lies of the boundary of unit ball
%
% Solve optimization problem:
%   minimize   ||x-y||
%      y                
%  subject to  ||inv(A)*y||^2 = 1
%
% Use fmincon to solve the nonlinear equality-constrained optimization problem (NLP)
% Analytic gradients of constraint and objective, and Hessian of Lagrangian are provided 

    % Should diagnostic information be printed?
    if nargin == 2
        verbose = false;
    else
        verbose = varargin{1};
    end

    % If x is exterior to the ellipsoid
    if norm(A\x) >= 1

        [y,dist_val] = geom.sign_dist_ellip_ext(x,A,verbose);

    else % x is in the interior of the ellipsoid

        n = length(x);
        I = eye(n);
        invP = (A\I)'/A;

        if verbose
            display_flag = 'iter';
        else
            display_flag = 'none';
        end
        
        opts = optimoptions('fmincon','Algorithm','interior-point',...
                            'SpecifyConstraintGradient',true,'SpecifyObjectiveGradient',true,...
                            'CheckGradients',false,...
                            'MaxIterations',1000,...
                            'HessianFcn',@(y,lambda) hessianfunc(y,lambda,P),'Display',display_flag);
        [y,obj_val,exit_flag] = fmincon(@(y) obj_fun(y,x),x+0.1*rand(n,1),[],[],[],[],[],[],@(y) nonlin_con(y,P),opts); 
        
        switch exit_flag
            case {0,-1,-2}
                error('Problem is unsolved.');
        end
    
        dist_val = -sqrt(obj_val); % It is known that x is in the interior of the ellipsoid

    end

end

% Functions required for solving the NLP by calling fmincon

function [obj_val,obj_grad] = obj_fun(y,x)
    obj_val = norm(y-x)^2;
    obj_grad = 2*(y-x);
end

function [c,ceq,dc,dceq] = nonlin_con(y,P)
    ceq = y'*P*y - 1;
    c = [];
    dc = [];
    dceq = 2*P*y;
end

function H = hessianfunc(y,lambda,P)
    n = length(y);
    H = 2*(eye(n) + lambda.eqnonlin*P);
end
