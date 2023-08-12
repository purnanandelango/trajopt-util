function [xbar,ubar,pbar,converged] = run_ptr(xbar,ubar,pbar,prb,sys_constr_cost_fun)
% PTR SCP with parameters as decision variables (including time dilation) and ZOH/FOH discretization

    converged = false;
    K = prb.K;

    % Check if type of FOH computation is specified
    if isfield(prb,'foh_type')
        foh_type = string(prb.foh_type);
        assert(ismember(foh_type,["v1","v2","v3"]),"Incorrect type of FOH discretization.");
    else
        foh_type = "v3";
    end
    
    fprintf("+------------------------------------------------------------------------------------------------------+\n");
    fprintf("|                                   ..:: Penalized Trust Region ::..                                   |\n");
    fprintf("+-------+------------+-----------+-----------+---------+---------+----------+---------+----------------+\n");
    fprintf("| Iter. | Prop. [ms] | Prs. [ms] | Slv. [ms] | log(TR) | log(VC) |   Cost   |   ToF   | log(VC cnstr.) |\n");
    fprintf("+-------+------------+-----------+-----------+---------+---------+----------+---------+----------------+\n");

    for j = 1:prb.scp_iters
        
        yalmip clear

        % Variables
        x = sdpvar(prb.nx,K);
        u = sdpvar(prb.nu,K);
        p = sdpvar(prb.np,1);
        vc_minus = sdpvar(prb.nx,K-1);
        vc_plus = sdpvar(prb.nx,K-1);
        
        Jvc = 0;
        cnstr = [];
        for k = 1:K-1
            % Virtual control penalty            
            Jvc = Jvc + sum(vc_plus(:,k) + vc_minus(:,k));
            cnstr = [cnstr;
                     vc_plus(:,k) >= 0;
                     vc_minus(:,k) >= 0];
        end

        % Trust region penalty
        switch prb.tr_norm
            case {2,inf}
                Jtr = sdpvar(1,K+1);        
                for k = 1:K
                    zbar_scl = [prb.invSx*(xbar(:,k)-prb.cx);
                                prb.invSu*(ubar(:,k)-prb.cu)];                    
                    cnstr = [cnstr; norm([x(:,k);u(:,k)]-zbar_scl,prb.tr_norm) <= Jtr(k)];        
                end
                pbar_scl = prb.invSp*(pbar-prb.cp);
                cnstr = [cnstr;norm(p-pbar_scl) <= Jtr(K+1)]; 
            case 'quad'
                Jtr = 0;        
                for k = 1:K
                    zbar_scl = [prb.invSx*(xbar(:,k)-prb.cx);
                                prb.invSu*(ubar(:,k)-prb.cu)];
                    Jtr = Jtr + ([x(:,k);u(:,k)]-zbar_scl)'*([x(:,k);u(:,k)]-zbar_scl);        
                end
                pbar_scl = prb.invSp*(pbar-prb.cp);
                Jtr = Jtr + (p-pbar_scl)'*(p-pbar_scl);
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
                         vc_plus(:,k) - vc_minus(:,k) == - x(:,k+1) - prb.invSx*prb.cx +...
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
                         vc_plus(:,k) - vc_minus(:,k) == - x(:,k+1) - prb.invSx*prb.cx +...
                                                           prb.invSx*Ak(:,:,k)*(prb.Sx*x(:,k)+prb.cx) +...
                                                           prb.invSx*Bk(:,:,k)*(prb.Su*u(:,k)+prb.cu) +...
                                                           prb.invSx*Sk(:,k)*(prb.Sp*p+prb.cp) +...
                                                           prb.invSx*wk(:,k)];
            end
            cnstr = [cnstr; u(:,K) == u(:,K-1)];             
        end        

        % Constraints
        [cnstr_sys,cost_fun,vc_constr_term] = sys_constr_cost_fun(x,u,p,prb,...
                                                                  xbar,ubar,pbar);
        
        cnstr = [cnstr;cnstr_sys];
        
        % Objective
        obj_fun = prb.wvc*sum(Jvc) + prb.wtr*sum(Jtr) + cost_fun;            
        
        % Solve
        % Model = export(cnstr,obj_fun,prb.solver_settings); % Export input to solver from YALMIP        
        yalmip_out = optimize(cnstr,obj_fun,prb.solver_settings); 
        % assert(ismember(yalmip_out.problem,[0,3]),"Subproblem is unsolved.\nSolver message: %s",yalmiperror(yalmip_out.problem));
        if ~ismember(yalmip_out.problem,[0,4])
            fprintf("+------------------------------------------------------------------------------------------------------+\n");
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
        cost_val = value(cost_fun);
        vc_term = sum(value(Jvc));
        if isfield(prb,'wvb')
            vc_constr_term = value(vc_constr_term)/prb.wvb;
        else
            vc_constr_term = value(vc_constr_term)/prb.wvc;
        end
        
        % Ensure that the TR value is always displayed consistently with 2-norm
        % Note that this is for display and for evaluation of termination criteria 
        switch prb.tr_norm
            case 2
                tr_term = sum(value(Jtr));
            case {'quad',inf}
                Jtr_post_solve = zeros(1,K+1);        
                for k = 1:K
                    zbar_scl = [prb.invSx*(xbar(:,k)-prb.cx);
                                prb.invSu*(ubar(:,k)-prb.cu)];                    
                    Jtr_post_solve(k) = norm([x(:,k);u(:,k)]-zbar_scl,2);        
                end
                pbar_scl = prb.invSp*(pbar-prb.cp);
                Jtr_post_solve(K+1) = norm(p-pbar_scl);
                tr_term = sum(Jtr_post_solve);
        end

        % Update reference trajectory
        for k = 1:prb.K
            xbar(:,k) = prb.Sx*x(:,k) + prb.cx;
            ubar(:,k) = prb.Su*u(:,k) + prb.cu;
        end
        pbar = prb.Sp*p + prb.cp;        

        ToF = prb.time_of_maneuver(xbar,ubar,pbar);
        
        % Console output
        fprintf('|  %02d   |   %7.1e  |  %7.1e  |  %7.1e  | %5.1f   | %5.1f   | %8.1e | %7.1e |    %5.1f       |\n',j,propagate_time,parse_time,solve_time,log10(tr_term),log10(vc_term),cost_val,ToF,log10(vc_constr_term))
        
        if vc_term < prb.epsvc && tr_term < prb.epstr
            converged = true;
            fprintf("+------------------------------------------------------------------------------------------------------+\n")
            fprintf('Converged!\n') 
            break
        end
        
    end

end