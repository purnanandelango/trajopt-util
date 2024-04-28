function [xbar,ubar,pbar,cost_val,converged] = ctscvx(xbar,ubar,pbar,prb,sys_constr_cost_fun)
% ct-SCvx with parameters as decision variables (including time dilation) and ZOH/FOH discretization
% Exact penalty weight can be matrix-valued

    converged = false;
    K = prb.K;

    % Check if type of FOH computation is specified
    if isfield(prb,'foh_type')
        foh_type = string(prb.foh_type);
        assert(ismember(foh_type,["v1","v2","v3"]),"Incorrect type of FOH discretization.");
    else
        foh_type = "v3";
    end

    % Exact penalty weight
    if isfield(prb,'w_ep')
        expnwt =  prb.w_ep;
    elseif isfield(prb,'W_ep')
        expnwt = prb.W_ep;
    end    
    
    fprintf("+--------------------------------------------------------------------------------------------------------+\n");
    fprintf("|                                          ..::   ct-SCvx   ::..                                         |\n");
    fprintf("|                 Successive Convexification with Continuous-Time Constraint Satisfaction                |\n");
    fprintf("+-------+------------+-----------+-----------+---------+---------+------------+---------+----------------+\n");
    fprintf("| Iter. | Prop. [ms] | Prs. [ms] | Slv. [ms] | log(px) | log(ep) |    Cost    |   ToF   | log(ep cnstr.) |\n");
    fprintf("+-------+------------+-----------+-----------+---------+---------+------------+---------+----------------+\n");

    for j = 1:prb.scp_iters
        
        yalmip clear

        % Variables
        x = sdpvar(prb.nx,K);
        u = sdpvar(prb.nu,K);
        p = sdpvar(prb.np,1);
        ep_minus = sdpvar(prb.nx,K-1);
        ep_plus = sdpvar(prb.nx,K-1);

        % Unscaled state and control input
        x_unscl = sdpvar(prb.nx,K);
        u_unscl = sdpvar(prb.nu,K);
        p_unscl = prb.Sp*p + prb.cp;
        for k = 1:K
            x_unscl(:,k) = prb.Sx*x(:,k) + prb.cx;
            u_unscl(:,k) = prb.Su*u(:,k) + prb.cu;
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
        xubar_scl = sdpvar(prb.nx+prb.nu,K);
        switch prb.px_norm
            case {2,inf}
                Jpx = sdpvar(1,K+1);        
                for k = 1:K
                    xubar_scl(:,k) = [prb.invSx*(xbar(:,k)-prb.cx);
                                      prb.invSu*(ubar(:,k)-prb.cu)];                    
                    cnstr = [cnstr; norm([x(:,k);u(:,k)]-xubar_scl(:,k),prb.px_norm) <= Jpx(k)];        
                end
                pbar_scl = prb.invSp*(pbar-prb.cp);
                cnstr = [cnstr;norm(p-pbar_scl) <= Jpx(K+1)]; 
            case 'quad'
                Jpx = 0;        
                for k = 1:K
                    xubar_scl(:,k) = [prb.invSx*(xbar(:,k)-prb.cx);
                                      prb.invSu*(ubar(:,k)-prb.cu)];
                    Jpx = Jpx + 0.5*([x(:,k);u(:,k)]-xubar_scl(:,k))'*([x(:,k);u(:,k)]-xubar_scl(:,k));        
                end
                pbar_scl = prb.invSp*(pbar-prb.cp);
                Jpx = Jpx + (p-pbar_scl)'*(p-pbar_scl);
        end

        % Linearized dynamics constraint
        if prb.disc == "FOH"
            % Propagation
            tic
            if isfield(prb,'ode_solver')
                [Ak,Bmk,Bpk,Sk,wk] = feval("disc.compute_foh_"+foh_type,prb.tau,xbar,ubar,pbar,prb.h,prb.dyn_func,prb.dyn_func_linearize,prb.ode_solver);
            else 
                [Ak,Bmk,Bpk,Sk,wk] = feval("disc.compute_foh_"+foh_type,prb.tau,xbar,ubar,pbar,prb.h,prb.dyn_func,prb.dyn_func_linearize);    
            end
            propagate_time = toc*1000;

            for k = 1:K-1
                cnstr = [cnstr;
                         ep_plus(:,k) - ep_minus(:,k) == - x(:,k+1) - prb.invSx*prb.cx +...
                                                           prb.invSx*Ak(:,:,k)*(prb.Sx*x(:,k)+prb.cx) +...
                                                           prb.invSx*Bmk(:,:,k)*(prb.Su*u(:,k)+prb.cu) +...
                                                           prb.invSx*Bpk(:,:,k)*(prb.Su*u(:,k+1)+prb.cu) +...
                                                           prb.invSx*Sk(:,k)*(prb.Sp*p+prb.cp) +...
                                                           prb.invSx*wk(:,k)];
            end                        
        elseif prb.disc == "ZOH"
            % Propagation
            tic
            if isfield(prb,'ode_solver')
                [Ak,Bk,Sk,wk] = disc.compute_zoh(prb.tau,xbar,ubar,pbar,prb.h,prb.dyn_func,prb.dyn_func_linearize,prb.ode_solver);
            else
                [Ak,Bk,Sk,wk] = disc.compute_zoh(prb.tau,xbar,ubar,pbar,prb.h,prb.dyn_func,prb.dyn_func_linearize);
            end
            propagate_time = toc*1000;

            for k = 1:K-1
                cnstr = [cnstr;
                         ep_plus(:,k) - ep_minus(:,k) == - x(:,k+1) - prb.invSx*prb.cx +...
                                                           prb.invSx*Ak(:,:,k)*(prb.Sx*x(:,k)+prb.cx) +...
                                                           prb.invSx*Bk(:,:,k)*(prb.Su*u(:,k)+prb.cu) +...
                                                           prb.invSx*Sk(:,k)*(prb.Sp*p+prb.cp) +...
                                                           prb.invSx*wk(:,k)];
            end
            cnstr = [cnstr; u(:,K) == u(:,K-1)];             
        end        

        % Constraints
        [cnstr_sys,cost_fun,ep_constr_term] = sys_constr_cost_fun(x_unscl,u_unscl,p_unscl,prb,...
                                                                  xbar,ubar,pbar);
        
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
        x = value(x);
        u = value(u);
        p = value(p);
        x_unscl = value(x_unscl);
        u_unscl = value(u_unscl);
        p_unscl = value(p_unscl);        
        cost_val = value(cost_fun)/prb.cost_factor;
        ep_term = sum(value(Jep));
        ep_constr_term = value(ep_constr_term)/max(expnwt(:));
        
        % Ensure that the px value is always displayed consistently with infinity norm
        % Note that this is for display and for evaluation of termination criteria 
        Jpx_post_solve = zeros(1,K+1);        
        for k = 1:K
            xubar_scl = [prb.invSx*(xbar(:,k)-prb.cx);
                        prb.invSu*(ubar(:,k)-prb.cu)];                    
            Jpx_post_solve(k) = norm([x(:,k);u(:,k)]-xubar_scl,'inf');        
        end
        pbar_scl = prb.invSp*(pbar-prb.cp);
        Jpx_post_solve(K+1) = norm(p-pbar_scl,'inf');
        px_term = max(Jpx_post_solve);

        % Update reference trajectory
        xbar = x_unscl;
        ubar = u_unscl;
        pbar = p_unscl;        

        ToF = prb.time_of_maneuver(xbar,ubar,pbar);
        
        % Console output
        fprintf('|  %02d   |   %7.1e  |  %7.1e  |  %7.1e  | %5.1f   | %5.1f   | %10.3e | %7.1e |    %5.1f       |\n',j,propagate_time,parse_time,solve_time,log10(px_term),log10(ep_term),cost_val,ToF,log10(ep_constr_term))
        
        if ep_term < prb.eps_ep && px_term < prb.eps_px
            converged = true;
            fprintf("+--------------------------------------------------------------------------------------------------------+\n")
            fprintf('Converged!\n') 
            break
        end
        
    end

end