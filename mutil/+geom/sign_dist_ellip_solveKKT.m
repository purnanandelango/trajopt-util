function [y,dist_val] = sign_dist_ellip_solveKKT(x,A,varargin)
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
% Use fsolve to compute the root of the function g(lam) = L(y,lam) s.t. y'*inv(P)*y = 1
% This method is not robust because the computed root might not correspond
% to a local minimum. KKT conditions are only necessary conditions for optimality

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

        if verbose
            display_flag = 'iter';
        else
            display_flag = 'none';
        end            

        In = eye(length(x));
        invP = (A\In)'/A;
        
        fun_y           = @(lam) (In - lam*invP)\x;
        fun_solve       = @(lam) fun_y(lam)'*(invP*fun_y(lam)) -  1;
        fun_solve_grad  = @(lam) 2*fun_y(lam)'*invP'*(inv(In - lam*invP)')*invP*fun_y(lam);    
        
        opts = optimoptions('fsolve','Algorithm','levenberg-marquardt',...
                            'FunctionTolerance',1e-9,'StepTolerance',1e-8,'Display',display_flag,...
                            'SpecifyObjectiveGradient',true);
    
        % Solve the KKT conditions multiple times and heuristically detect local minimum
        
        % Maximum number of times fsolve is called
        max_outer_iter = 20;
        cntr = 1; % Counter
    
        lam_curr = 1e3;
        dist_val_curr = 1e6;
    
        % This while loop will successfully terminate only if the local minimum
        % is found immediately before or after the local maximum
    
        while cntr <= max_outer_iter
            
            lam_prev = lam_curr;
    
            % Solve KKT conditions with uniform random initial guess in (0,1)
            [lam_curr,~,exit_flag] = fsolve(@(lam) fun_solve_main(lam,fun_solve,fun_solve_grad),rand(),opts);
    
            % Exit loop if fsolve is successful and if at least 5 iterations have occured
            if exit_flag == 1 && cntr >= 5
                break
            end
    
            % Current primal value (projection point)
            y = fun_y(lam_curr);
            dist_val_prev = dist_val_curr;
            % Current distance between x and projection point        
            dist_val_curr = norm(x-y);
    
            if verbose
                if cntr == 1
                    fprintf("\n")
                end
                fprintf("Itr. %02d projection distance: %.2f\n",cntr,dist_val_curr);
            end
            
            % If the current iterate distance exceeds the previous iterate
            % distance then exit loop
            if dist_val_curr > dist_val_prev
                lam_curr = lam_prev;
                break
            end
            % After at least 3 iterations if the previous iterate distance exceeds the current iterate
            % distance after then exit loop        
            if cntr>=2 && dist_val_prev > dist_val_curr
                break
            end
            cntr = cntr+1;
        end
        if exit_flag <= 0
            error("Signed-distance is not computed.")
        end
    
        y = fun_y(lam_curr);
        dist_val = -norm(x-y); % Since it is known that x is inside the ellipsoid

    end
end

function [f,df] = fun_solve_main(lam,fun_solve,fun_solve_grad)
    f = fun_solve(lam);
    df = fun_solve_grad(lam);
end
