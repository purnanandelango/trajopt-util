function [xbar,ubar,cost_val,converged] = run_ptr_noparam(xbar,ubar,prb,sys_constr_cost_fun,varargin)
% PTR SCP without parameters as decision variables and ZOH/FOH discretization
% Provision for updating problem parameters after each SCP iteration
% Exact penalty weight can be matrix-valued

    converged = false;
    K = prb.K;

    % Check if type of FOH computation is specified
    if isfield(prb,'foh_type')
        foh_type = string(prb.foh_type);
        assert(ismember(foh_type,["v1","v2","v3","v3_parallel"]),"Incorrect type of FOH discretization.");        
    else
        foh_type = "v3"; % Default
    end

    % Exact penalty weight
    if isfield(prb,'wvc')
        expnwt =  prb.wvc;
    elseif isfield(prb,'Wvc')
        expnwt = prb.Wvc;
    end
    
    fprintf("+--------------------------------------------------------------------------------------------------------+\n");
    fprintf("|                                    ..:: Penalized Trust Region ::..                                    |\n");
    fprintf("+-------+------------+-----------+-----------+---------+---------+------------+---------+----------------+\n");
    fprintf("| Iter. | Prop. [ms] | Prs. [ms] | Slv. [ms] | log(TR) | log(VC) |    Cost    |   ToF   | log(VC cnstr.) |\n");
    fprintf("+-------+------------+-----------+-----------+---------+---------+------------+---------+----------------+\n");
    
    for j = 1:prb.scp_iters
        
        yalmip clear

        % Variables
        x = sdpvar(prb.nx,K);
        u = sdpvar(prb.nu,K);
        vc_minus = sdpvar(prb.nx,K-1);
        vc_plus = sdpvar(prb.nx,K-1);
    
        % Unscaled state and control input
        x_unscl = sdpvar(prb.nx,K);
        u_unscl = sdpvar(prb.nu,K);
        for k = 1:K
            x_unscl(:,k) = prb.Sx*x(:,k) + prb.cx;
            u_unscl(:,k) = prb.Su*u(:,k) + prb.cu;
        end
        
        Jvc = 0;
        cnstr = [];
        for k = 1:K-1
            % Virtual control penalty 
            Jvc = Jvc + sum(expnwt*(vc_plus(:,k) + vc_minus(:,k)));
            cnstr = [cnstr;
                     vc_plus(:,k) >= 0;
                     vc_minus(:,k) >= 0];
        end

        % Trust region penalty
        xubar_scl = sdpvar(prb.nx+prb.nu,K);
        switch prb.tr_norm
            case {2,inf}
                Jtr = sdpvar(1,K);        
                for k = 1:K
                    xubar_scl(:,k) = [prb.invSx*(xbar(:,k)-prb.cx);
                                      prb.invSu*(ubar(:,k)-prb.cu)];                    
                    cnstr = [cnstr; norm([x(:,k);u(:,k)]-xubar_scl(:,k),prb.tr_norm) <= Jtr(k)]; 
                end                
            case 'quad'
                Jtr = 0;        
                for k = 1:K
                    xubar_scl(:,k) = [prb.invSx*(xbar(:,k)-prb.cx);
                                      prb.invSu*(ubar(:,k)-prb.cu)];
                    Jtr = Jtr + 0.5*([x(:,k);u(:,k)]-xubar_scl(:,k))'*([x(:,k);u(:,k)]-xubar_scl(:,k));
                end                
        end

        % Linearized dynamics constraint
        if prb.disc == "FOH"
            % Propagation
            tic
            if isfield(prb,'ode_solver')
                [Ak,Bmk,Bpk,wk] = feval("disc.compute_foh_noparam_"+foh_type,prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize,prb.ode_solver);
            else
                [Ak,Bmk,Bpk,wk] = feval("disc.compute_foh_noparam_"+foh_type,prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize);
            end
            propagate_time = toc*1000;

            for k = 1:K-1
                cnstr = [cnstr;
                         vc_plus(:,k) - vc_minus(:,k) == - x(:,k+1) - prb.invSx*prb.cx +...
                                                           prb.invSx*Ak(:,:,k)*(prb.Sx*x(:,k)+prb.cx) +...
                                                           prb.invSx*Bmk(:,:,k)*(prb.Su*u(:,k)+prb.cu) +...
                                                           prb.invSx*Bpk(:,:,k)*(prb.Su*u(:,k+1)+prb.cu) +...
                                                           prb.invSx*wk(:,k)];
            end                        
        elseif prb.disc == "ZOH"
            % Propagation
            tic
            if isfield(prb,'ode_solver')
                [Ak,Bk,wk] = disc.compute_zoh_noparam(prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize,prb.ode_solver);
            else
                [Ak,Bk,wk] = disc.compute_zoh_noparam(prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize);
            end
            propagate_time = toc*1000;

            for k = 1:K-1
                cnstr = [cnstr;
                         vc_plus(:,k) - vc_minus(:,k) == - x(:,k+1) - prb.invSx*prb.cx +...
                                                           prb.invSx*Ak(:,:,k)*(prb.Sx*x(:,k)+prb.cx) +...
                                                           prb.invSx*Bk(:,:,k)*(prb.Su*u(:,k)+prb.cu) +...
                                                           prb.invSx*wk(:,k)];
            end
            cnstr = [cnstr; u(:,K) == u(:,K-1)];             
        elseif prb.disc == "Impulse"
            % Propagation
            tic
            if isfield(prb,'ode_solver')
                [Ak,wk] = disc.compute_impulse_noparam(prb.tau,xbar,ubar,prb.Eu2x,prb.h,prb.dyn_func,prb.dyn_func_linearize,prb.ode_solver);
            else
                [Ak,wk] = disc.compute_impulse_noparam(prb.tau,xbar,ubar,prb.Eu2x,prb.h,prb.dyn_func,prb.dyn_func_linearize);            
            end
            propagate_time = toc*1000;

            for k = 1:K-1
                cnstr = [cnstr;
                         vc_plus(:,k) - vc_minus(:,k) == - x(:,k+1) - prb.invSx*prb.cx +...
                                                           prb.invSx*Ak(:,:,k)*(prb.Sx*x(:,k)+prb.cx) +...
                                                           prb.invSx*Ak(:,:,k)*prb.Eu2x*(prb.Su*u(:,k)+prb.cu) +...
                                                           prb.invSx*wk(:,k)];                
            end
            cnstr = [cnstr; u(:,K) == u(:,K-1)];  
        end
        
        % Constraints
        [cnstr_sys,cost_fun,vc_constr_term] = sys_constr_cost_fun(x_unscl,u_unscl,prb,...
                                                                  xbar,ubar);


        cnstr = [cnstr;cnstr_sys];
        
        % Objective
        obj_fun = Jvc + prb.wtr*sum(Jtr) + cost_fun;            
        
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
        x_unscl = value(x_unscl);
        u_unscl = value(u_unscl);
        cost_val = value(cost_fun)/prb.cost_factor;
        vc_term = value(Jvc);
        vc_constr_term = value(vc_constr_term)/max(expnwt(:));

        % Ensure that the TR value is always displayed consistently with infinity norm
        % Note that this is for display and for evaluation of termination criteria 
        Jtr_post_solve = zeros(1,K);        
        for k = 1:K
            xubar_scl = [prb.invSx*(xbar(:,k)-prb.cx);
                        prb.invSu*(ubar(:,k)-prb.cu)];                    
            Jtr_post_solve(k) = norm([x(:,k);u(:,k)]-xubar_scl,'inf');
        end
        tr_term = max(Jtr_post_solve);

        % Update reference trajectory
        xbar = x_unscl;
        ubar = u_unscl;

        ToF = prb.time_of_maneuver(xbar,ubar);        
        
        % Console output
        fprintf('|  %02d   |   %7.1e  |  %7.1e  |  %7.1e  | %5.1f   | %5.1f   | %10.3e | %7.1e |    %5.1f       |\n',j,propagate_time,parse_time,solve_time,log10(tr_term),log10(vc_term),cost_val,ToF,log10(vc_constr_term));
        
        if vc_term < prb.epsvc && tr_term < prb.epstr
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