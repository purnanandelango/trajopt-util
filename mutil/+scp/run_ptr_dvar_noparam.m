function [xbar,ubar,converged] = run_ptr_dvar_noparam(xbar,ubar,prb,sys_constr_cost_fun,varargin)
% PTR SCP without parameters as decision variables and ZOH/FOH discretization
% Provision for updating problem parameters after each SCP iteration
% Scaling terms cx and cu are not used

    converged = false;
    K = prb.K;

    % Check if type of FOH computation is specified
    if isfield(prb,'foh_type')
        foh_type = string(prb.foh_type);
        assert(ismember(foh_type,["v1","v2","v3"]),"Incorrect type of FOH discretization.");        
    else
        foh_type = "v3"; % Default
    end
    
    fprintf("+------------------------------------------------------------------------------------------------------+\n");
    fprintf("|                                   ..:: Penalized Trust Region ::..                                   |\n");
    fprintf("+-------+------------+-----------+-----------+---------+---------+----------+---------+----------------+\n");
    fprintf("| Iter. | Prop. [ms] | Prs. [ms] | Slv. [ms] | log(TR) | log(VC) |   Cost   |   ToF   | log(VC cnstr.) |\n");
    fprintf("+-------+------------+-----------+-----------+---------+---------+----------+---------+----------------+\n");
    
    for j = 1:prb.scp_iters
        
        yalmip clear

        % Variables
        dx = sdpvar(prb.nx,K);
        du = sdpvar(prb.nu,K);
        vc_minus = sdpvar(prb.nx,K-1);
        vc_plus = sdpvar(prb.nx,K-1);

        x = dx + prb.invSx*xbar;
        u = du + prb.invSu*ubar;
        
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
                Jtr = sdpvar(1,K);        
                for k = 1:K
                    cnstr = [cnstr; norm([dx(:,k);du(:,k)],prb.tr_norm) <= Jtr(k)]; 
                end                
            case 'quad'
                Jtr = 0;        
                for k = 1:K
                    Jtr = Jtr + ([dx(:,k);du(:,k)])'*([dx(:,k);du(:,k)]);
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
                cnstr = [cnstr;
                         vc_minus(:,k) - vc_plus(:,k) == - dx(:,k+1) +...
                                                           prb.invSx*Ak(:,:,k)*prb.Sx*dx(:,k) +...
                                                           prb.invSx*Bmk(:,:,k)*prb.Su*du(:,k) +...
                                                           prb.invSx*Bpk(:,:,k)*prb.Su*du(:,k+1) +...
                                                           prb.invSx*(xbarprop(:,k+1) - xbar(:,k+1))];
            end                        
        elseif prb.disc == "ZOH"
            % Propagation
            tic
            if isfield(prb,'ode_solver')
                [Ak,Bk,~,~,xbarprop] = disc.compute_zoh_noparam(prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize,prb.ode_solver);
            else
                [Ak,Bk,~,~,xbarprop] = disc.compute_zoh_noparam(prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize);
            end
            propagate_time = toc*1000;

            for k = 1:K-1
                cnstr = [cnstr;
                         vc_minus(:,k) - vc_plus(:,k) == - dx(:,k+1) +...
                                                           prb.invSx*Ak(:,:,k)*prb.Sx*dx(:,k) +...
                                                           prb.invSx*Bk(:,:,k)*prb.Su*du(:,k) +...
                                                           prb.invSx*(xbarprop(:,k+1) - xbar(:,k+1))];
            end
            cnstr = [cnstr; u(:,K) == u(:,K-1)];             
        end
        
        % Constraints
        [cnstr_sys,cost_fun,vc_constr_term] = sys_constr_cost_fun(x,u,prb,...
                                                                  xbar,ubar);


        cnstr = [cnstr;cnstr_sys];
        
        % Objective
        obj_fun = prb.wvc*Jvc + prb.wtr*sum(Jtr) + cost_fun;            
        
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
        % x = value(x);
        % u = value(u);
        dx = value(dx);
        du = value(du);
        cost_val = value(cost_fun)/prb.cost_factor;
        vc_term = value(Jvc);
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
                Jtr_post_solve = zeros(1,K);        
                for k = 1:K
                    Jtr_post_solve(k) = norm([dx(:,k);du(:,k)],2);        
                end
                tr_term = sum(Jtr_post_solve);
        end

        % Update reference trajectory
        for k = 1:prb.K
            xbar(:,k) = prb.Sx*dx(:,k) + xbar(:,k);
            ubar(:,k) = prb.Su*du(:,k) + ubar(:,k);
        end        

        ToF = prb.time_of_maneuver(xbar,ubar);        
        
        % Console output
        fprintf('|  %02d   |   %7.1e  |  %7.1e  |  %7.1e  | %5.1f   | %5.1f   | %8.1e | %7.1e |    %5.1f       |\n',j,propagate_time,parse_time,solve_time,log10(tr_term),log10(vc_term),cost_val,ToF,log10(vc_constr_term));
        
        if vc_term < prb.epsvc && tr_term < prb.epstr
            converged = true;
            fprintf("+------------------------------------------------------------------------------------------------------+\n")
            fprintf('Converged!\n')
            break
        end

        if nargin == 5 && j < prb.scp_iters % Update problem parameters
            prb = varargin{1}(prb,xbar,ubar);
        end
        
    end

end