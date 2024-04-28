function [xbar,ubar,cost_val,converged] = ctscvx_dvar_noparam(xbar,ubar,prb,sys_constr_cost_fun,varargin)
% ct-SCvx without parameters as decision variables and ZOH/FOH discretization
% Provision for updating problem parameters after each SCP iteration
% Scaling terms cx and cu are not used
% Exact penalty weight can be matrix-valued

    converged = false;
    K = prb.K;

    % assert(max(prb.cx == zeros(prb.nx,1)) && max(prb.cu == zeros(prb.nu,1)),"Scaling parameters cx and cu should be 0 for deviation variables.");
    if ~(max(prb.cx == zeros(prb.nx,1)) && max(prb.cu == zeros(prb.nu,1)))
        warning("Affine scaling offset terms cx and cu are nonzero.")
    end

    % Check if type of discretization computation is specified
    if isfield(prb,'foh_type')
        foh_type = string(prb.foh_type);
        assert(ismember(foh_type,["v1","v2","v3","v3_parallel"]),"Incorrect type of FOH discretization.");        
    else
        foh_type = "v3"; % Default
    end

    if isfield(prb,'zoh_type')
        zoh_type = string(prb.zoh_type);
        assert(ismember(zoh_type,["v1","v3","v3_parallel"]),"Incorrect type of ZOH discretization.");  
    else
        zoh_type = "v3"; % Default
    end

    % Exact penalty weight
    if isfield(prb,'w_ep')
        expnwt =  prb.w_ep;
    elseif isfield(prb,'W_ep')
        expnwt = prb.W_ep;
    end    
    
    fprintf("+--------------------------------------------------------------------------------------------------------+\n");
    fprintf("|                           ..::   ct-SCvx - Successive Convexification   ::..                           |\n");
    fprintf("+-------+------------+-----------+-----------+---------+---------+------------+---------+----------------+\n");
    fprintf("| Iter. | Prop. [ms] | Prs. [ms] | Slv. [ms] | log(px) | log(ep) |    Cost    |   ToF   | log(ep cnstr.) |\n");
    fprintf("+-------+------------+-----------+-----------+---------+---------+------------+---------+----------------+\n");
    
    for j = 1:prb.scp_iters
        
        yalmip clear

        % Variables
        dx = sdpvar(prb.nx,K);
        du = sdpvar(prb.nu,K);
        ep_minus = sdpvar(prb.nx,K-1);
        ep_plus = sdpvar(prb.nx,K-1);

        % Unscaled state and control input
        x_unscl = sdpvar(prb.nx,K);
        u_unscl = sdpvar(prb.nu,K);
        for k = 1:K
            x_unscl(:,k) = prb.Sx*dx(:,k) + xbar(:,k);
            u_unscl(:,k) = prb.Su*du(:,k) + ubar(:,k);
        end
        
        Jep = 0;
        cnstr = [];
        for k = 1:K-1
            % Virtual control penalty 
            Jep = Jep + sum(expnwt*(ep_plus(:,k) + ep_minus(:,k)));
            cnstr = [cnstr;
                     ep_plus(:,k) >= 0;
                     ep_minus(:,k) >= 0];
        end

        % Trust region penalty
        switch prb.px_norm
            case {2,inf}
                Jpx = sdpvar(1,K);        
                for k = 1:K
                    cnstr = [cnstr; norm([dx(:,k);du(:,k)],prb.px_norm) <= Jpx(k)]; 
                end                
            case 'quad'
                Jpx = 0;        
                for k = 1:K
                    Jpx = Jpx + 0.5*([dx(:,k);du(:,k)])'*([dx(:,k);du(:,k)]);
                end                
        end

        % Linearized dynamics constraint
        if prb.disc == "FOH"
            % Propagation
            tic
            if isfield(prb,'ode_solver')
                [Ak,Bmk,Bpk,~,~,xbarprop] = feval("disc.compute_foh_noparam_"+foh_type,prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize,prb.ode_solver);
            else
                [Ak,Bmk,Bpk,~,~,xbarprop] = feval("disc.compute_foh_noparam_"+foh_type,prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize);
            end
            propagate_time = toc*1000;

            for k = 1:K-1
                Ahatk   = prb.invSx*Ak(:,:,k)*prb.Sx;
                Bmhatk  = prb.invSx*Bmk(:,:,k)*prb.Su;
                Bphatk  = prb.invSx*Bpk(:,:,k)*prb.Su;
                whatk   = prb.invSx*(xbarprop(:,k+1) - xbar(:,k+1)); 

                % Row normalization
                % scl_mat = eye(prb.nx);                
                scl_mat = diag(vecnorm(...
                                       [-eye(prb.nx), eye(prb.nx), -eye(prb.nx), Ahatk, Bmhatk, Bphatk, whatk], ...
                                       Inf,2));                

                cnstr = [cnstr;
                         zeros(prb.nx,1) == scl_mat\(ep_minus(:,k) - ep_plus(:,k) - dx(:,k+1) +...
                                                     Ahatk*dx(:,k) +...
                                                     Bmhatk*du(:,k) +...
                                                     Bphatk*du(:,k+1) +...
                                                     whatk)];
            end                        
        elseif prb.disc == "ZOH"
            % Propagation
            tic
            if isfield(prb,'ode_solver')
                [Ak,Bk,~,~,xbarprop] = feval("disc.compute_zoh_noparam_"+zoh_type,prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize,prb.ode_solver);
            else
                [Ak,Bk,~,~,xbarprop] = feval("disc.compute_zoh_noparam_"+zoh_type,prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize);
            end
            propagate_time = toc*1000;

            for k = 1:K-1
                Ahatk   = prb.invSx*Ak(:,:,k)*prb.Sx;
                Bhatk  = prb.invSx*Bk(:,:,k)*prb.Su;
                whatk   = prb.invSx*(xbarprop(:,k+1) - xbar(:,k+1)); 

                % Row normalization
                % scl_mat = eye(prb.nx);                
                scl_mat = diag(vecnorm(...
                                       [-eye(prb.nx), eye(prb.nx), -eye(prb.nx), Ahatk, Bhatk, whatk], ...
                                       Inf,2));                

                cnstr = [cnstr;
                         zeros(prb.nx,1) == scl_mat\(ep_minus(:,k) - ep_plus(:,k) - dx(:,k+1) +...
                                                     Ahatk*dx(:,k) +...
                                                     Bhatk*du(:,k) +...
                                                     whatk)];
            end
            cnstr = [cnstr; du(:,K) == du(:,K-1)];             
        end
        
        % Constraints
        [cnstr_sys,cost_fun,ep_constr_term] = sys_constr_cost_fun(x_unscl,u_unscl,prb,...
                                                                  xbar,ubar);


        cnstr = [cnstr;cnstr_sys];
        
        % Objective
        obj_fun = Jep + prb.w_px*sum(Jpx) + cost_fun;            
        
        % Solve
        % Model = export(cnstr,obj_fun,prb.solver_settings); % Export input to solver from YALMIP
        yalmip_out = optimize(cnstr,obj_fun,prb.solver_settings);
        % assert(ismember(yalmip_out.problem,[0,3]),"Subproblem is unsolved.\nSolver message: %s",yalmiperror(yalmip_out.problem));
        if ~ismember(yalmip_out.problem,[0,4])
            fprintf("+--------------------------------------------------------------------------------------------------------+\n");
            fprintf('Subproblem is unsolved. Returning the previous iterate.\n'); 
            break
        end
        if yalmip_out.problem == 4
            warning("Solver numerical issues.");
        end
        
        % Post process
        solve_time = yalmip_out.solvertime*1000;
        parse_time = yalmip_out.yalmiptime*1000;            
        dx = value(dx);
        du = value(du);
        x_unscl = value(x_unscl);
        u_unscl = value(u_unscl);        
        cost_val = value(cost_fun)/prb.cost_factor;
        ep_term = value(Jep);
        ep_constr_term = value(ep_constr_term)/max(expnwt(:));

        % Ensure that the px value is always displayed consistently with infinity norm
        % Note that this is for display and for evaluation of termination criteria 
        Jpx_post_solve = zeros(1,K);        
        for k = 1:K
            Jpx_post_solve(k) = norm([dx(:,k);du(:,k)],'inf');
        end
        px_term = max(Jpx_post_solve);

        % Update reference trajectory
        xbar = x_unscl;
        ubar = u_unscl;     

        ToF = prb.time_of_maneuver(xbar,ubar);        
        
        % Console output
        fprintf('|  %02d   |   %7.1e  |  %7.1e  |  %7.1e  | %5.1f   | %5.1f   | %10.3e | %7.1e |    %5.1f       |\n',j,propagate_time,parse_time,solve_time,log10(px_term),log10(ep_term),cost_val,ToF,log10(ep_constr_term));
        
        if ep_term < prb.eps_ep && px_term < prb.eps_px
            converged = true;
            fprintf("+--------------------------------------------------------------------------------------------------------+\n")
            fprintf('Converged!\n')
            break
        end

        if nargin == 5 && j < prb.scp_iters % Update problem parameters
            prb = varargin{1}(prb,xbar,ubar);
        end
        
    end

end