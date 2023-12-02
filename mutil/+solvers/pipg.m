function [z_jp1,w_jp1,v_jp1,status] = pipg(model, options, ...
                                           varargin)
% PIPG with extrapolation
%
% Model:
%   - nx
%   - nu
%   - K
%   - Phat
%   - phat
%   - Gtil
%   - gtil
%   - Htil
%   - htil
%   - scl_bnd
%   - i_idx
%   - f_idx
%   - zhat_i
%   - zhat_f
% 
% Options:
%   - alpha
%   - beta
%   - rho
%   - eps_abs
%   - max_iter
%   - verbose
%   - test_termination

    if nargin == 5 % Warmstart both primal and dual
        zet_j = varargin{1};
        eta_j = varargin{2};
        chi_j = varargin{3};
    elseif nargin == 3 % Warmstart only primal
        zet_j = varargin{1};
        eta_j = zeros(size(model.Gtil,1),1);
        chi_j = zeros(size(model.Htil,1),1);
    elseif nargin == 2 % Start from zeros
        zet_j = zeros(size(model.Phat,1),1);
        eta_j = zeros(size(model.Gtil,1),1);
        chi_j = zeros(size(model.Htil,1),1);
    else
        error("Incorrect number of arguments passed.")
    end

    % Convenience definitions
    nxK     = model.nx*model.K;
    nuK     = model.nu*model.K;
    nxKm1   = model.nx*(model.K-1);
    nxnuK   = (model.nx+model.nu)*model.K;
    uhatmin = model.scl_bnd(1)*ones(nuK,1);
    uhatmax = model.scl_bnd(2)*ones(nuK,1);

    GtilT = model.Gtil';
    HtilT = model.Htil';

    % Allocation
    z_jp1 = zet_j;
    w_jp1 = eta_j;
    v_jp1 = chi_j;

    if options.verbose
        fprintf("\n");
        fprintf("+-----------------------------------------------------------------------+\n");
        fprintf("|                             ..:: PIPG ::..                            |\n");
        fprintf("+-----------+-----------+------------------+--------------+-------------+\n");
        fprintf("| Iteration | Objective | Constraint Viol. | Primal Diff. |  Dual Diff. |\n");
        fprintf("+-----------+-----------+------------------+--------------+-------------+\n");
    end

    for j = 1:options.max_iter

        z_j = z_jp1;
        w_j = w_jp1;
        v_j = v_jp1;
       
        % Gradient descent
        z_jp1 = zet_j - options.alpha * (model.Phat*zet_j + model.phat + GtilT*eta_j + HtilT*chi_j);
       
        % Projection
        z_jp1(model.i_idx)       = model.zhat_i;                    % Initial condition
        z_jp1(nxKm1+model.f_idx) = model.zhat_f;                    % Final condition
        z_jp1(nxK+1:nxnuK)       = max(uhatmin,min(uhatmax,...
                                       z_jp1(nxK+1:nxnuK)));        % Control input bounds
        z_jp1(nxnuK+1:end)       = max(0,z_jp1(nxnuK+1:end));       % Slack nonnegativity

        % Gradient ascent and projection on polar cone
        w_jp1 = eta_j + options.beta*(model.Gtil*(2*z_jp1 - zet_j) - model.gtil);
        v_jp1 = max(0,chi_j + options.beta*(model.Htil*(2*z_jp1 - zet_j) - model.htil));

        % Extrapolation
        zet_j = (1-options.rho)*zet_j + options.rho*z_jp1;
        eta_j = (1-options.rho)*eta_j + options.rho*w_jp1;
        chi_j = (1-options.rho)*chi_j + options.rho*v_jp1;
        
        % Test termination criteria
        if mod(j,options.test_termination) == 0
    
            obj_val        = 0.5*z_jp1'*model.Phat*z_jp1 + model.phat'*z_jp1;
            cnstr_viol     = norm([model.Gtil*z_jp1 - model.gtil;
                                   max(0,model.Htil*z_jp1 - model.htil)],'inf');
            primal_diff    = norm(z_jp1-z_j,'inf');
            dual_diff_eq   = norm(w_jp1-w_j,'inf');
            dual_diff_ineq = norm(v_jp1-v_j,'inf');
            
            if options.verbose
                fprintf("| %9.2e | %9.2e |     %9.2e    |   %9.2e  |  %9.2e  |\n",j,obj_val,cnstr_viol,primal_diff,max(dual_diff_eq,dual_diff_ineq));
            end

            if primal_diff < options.eps_abs ...
               && dual_diff_eq < options.eps_abs ...
               && dual_diff_ineq < options.eps_abs

                if options.verbose
                    fprintf("+-----------+-----------+------------------+--------------+-------------+\n");
                    fprintf("Converged!\n")
                    status = "solved";
                end
                break
            
            end
        end

    end
    if j == options.max_iter
        if options.verbose
            fprintf("+-----------+-----------+------------------+--------------+-------------+\n");  
            fprintf("Max iterations.\n");
        end
        status = "unsolved";
    end
end