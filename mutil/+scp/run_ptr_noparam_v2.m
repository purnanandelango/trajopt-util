function [xbar,ubar,converged] = run_ptr_noparam_v2(xbar,ubar,prb,sys_constr_cost_fun,varargin)
% PTR SCP without parameters as decision variables and ZOH/FOH discretization
% Same as run_ptr_noparam except that VC weight is matrix valued
% Provision for updating problem parameters after each SCP iteration

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
        vc_minus = sdpvar(prb.nx,K-1);
        vc_plus = sdpvar(prb.nx,K-1);
        
        Jvc = 0;
        cnstr = [];
        for k = 1:K-1
            % Virtual control penalty            
            Jvc = Jvc + sum(prb.Wvc*(vc_plus(:,k) + vc_minus(:,k)));
            cnstr = [cnstr;
                     vc_plus(:,k) >= 0;
                     vc_minus(:,k) >= 0];
        end

        % Trust region penalty
        switch prb.tr_norm
            case {2,inf}
                Jtr = sdpvar(1,K);        
                for k = 1:K
                    zbar_scl = [prb.invSx*(xbar(:,k)-prb.cx);
                                prb.invSu*(ubar(:,k)-prb.cu)];                    
                    cnstr = [cnstr; norm([x(:,k);u(:,k)]-zbar_scl,prb.tr_norm) <= Jtr(k)]; 
                end                
            case 'quad'
                Jtr = 0;        
                for k = 1:K
                    zbar_scl = [prb.invSx*(xbar(:,k)-prb.cx);
                                prb.invSu*(ubar(:,k)-prb.cu)];
                    Jtr = Jtr + ([x(:,k);u(:,k)]-zbar_scl)'*([x(:,k);u(:,k)]-zbar_scl);
                end                
        end

        % Linearized dynamics constraint
        if prb.disc == "FOH"
            % Propagation
            tic
            [Ak,Bmk,Bpk,wk] = feval("disc.compute_foh_noparam_"+foh_type,prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize);
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
            [Ak,Bk,wk] = disc.compute_zoh_noparam(prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize);
            propagate_time = toc*1000;

            for k = 1:K-1
                cnstr = [cnstr;
                         vc_plus(:,k) - vc_minus(:,k) == - x(:,k+1) - prb.invSx*prb.cx +...
                                                           prb.invSx*Ak(:,:,k)*(prb.Sx*x(:,k)+prb.cx) +...
                                                           prb.invSx*Bk(:,:,k)*(prb.Su*u(:,k)+prb.cu) +...
                                                           prb.invSx*wk(:,k)];
            end
            cnstr = [cnstr; u(:,K) == u(:,K-1)];             
        end
        
        % Constraints
        [cnstr_sys,cost_fun,vc_constr_term] = sys_constr_cost_fun(x,u,prb,...
                                                                  xbar,ubar);


        cnstr = [cnstr;cnstr_sys];
        
        % Objective
        obj_fun = Jvc + prb.wtr*sum(Jtr) + cost_fun;            
        
        % Solve
        yalmip_out = optimize(cnstr,obj_fun,prb.solver_settings);
        % assert(ismember(yalmip_out.problem,[0,3]),"Subproblem is unsolved.\nSolver message: %s",yalmiperror(yalmip_out.problem));
        if ~ismember(yalmip_out.problem,[0,3])
            fprintf("+------------------------------------------------------------------------------------------------------+\n");
            fprintf('Subproblem is unsolved. Returning the previous iterate.\n'); 
            break
        end
        
        % Post process
        solve_time = yalmip_out.solvertime*1000;
        parse_time = yalmip_out.yalmiptime*1000;            
        x = value(x);
        u = value(u);
        cost_val = value(cost_fun);
        vc_term = 0;
        vc_constr_term = value(vc_constr_term);

        for k = 1:K-1
            % Virtual control penalty            
            vc_term = vc_term + sum(value(vc_plus(:,k)) + value(vc_minus(:,k)));
        end        


        % Ensure that the TR value is always displayed consistently with 2-norm
        % Note that this is for display and for evaluation of termination criteria 
        switch prb.tr_norm
            case 2
                tr_term = sum(value(Jtr));
            case {'quad',inf}
                Jtr_post_solve = zeros(1,K);        
                for k = 1:K
                    zbar_scl = [prb.invSx*(xbar(:,k)-prb.cx);
                                prb.invSu*(ubar(:,k)-prb.cu)];                    
                    Jtr_post_solve(k) = norm([x(:,k);u(:,k)]-zbar_scl,2);        
                end
                tr_term = sum(Jtr_post_solve);
        end

        % Update reference trajectory
        for k = 1:prb.K
            xbar(:,k) = prb.Sx*x(:,k) + prb.cx;
            ubar(:,k) = prb.Su*u(:,k) + prb.cu;
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