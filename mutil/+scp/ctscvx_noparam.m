function [xbar,ubar,cost_val,converged] = ctscvx_noparam(xbar,ubar,prb,sys_constr_cost_fun,varargin)
% ct-SCvx without parameters as decision variables and FOH/ZOH/FBP/Impulse discretization
% Provision for updating problem parameters after each SCP iteration
% Exact penalty weight can be matrix-valued

    converged = false;
    K = prb.K;

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

    if isfield(prb,'impulse_type')
        impulse_type = string(prb.impulse_type);
        assert(ismember(impulse_type,["v3","v3_parallel"]),"Incorrect type of impulse discretization.");  
    else
        impulse_type = "v3"; % Default
    end

    if isfield(prb,'fbp_type')
        fbp_type = string(prb.fbp_type);
        assert(ismember(fbp_type,["v3","v3_parallel"]),"Incorrect type of FBP discretization.");  
    else
        fbp_type = "v3"; % Default
    end    

    % Exact penalty weight
    if isfield(prb,'wvc')
        expnwt =  prb.wvc;
    elseif isfield(prb,'Wvc')
        expnwt = prb.Wvc;
    end
    
    fprintf("+--------------------------------------------------------------------------------------------------------+\n");
    fprintf("|                                          ..::   ct-SCvx   ::..                                         |\n");
    fprintf("|                 Successive Convexification with Continuous-Time Constraint Satisfaction                |");
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
                Ahatk   = prb.invSx*Ak(:,:,k)*prb.Sx;
                Bmhatk  = prb.invSx*Bmk(:,:,k)*prb.Su;
                Bphatk  = prb.invSx*Bpk(:,:,k)*prb.Su;
                whatk   = prb.invSx*(wk(:,k) + Ak(:,:,k)*prb.cx + Bmk(:,:,k)*prb.cu + Bpk(:,:,k)*prb.cu - prb.cx);

                % Row normalization
                % scl_mat = eye(prb.nx);                
                scl_mat = diag(vecnorm(...
                                       [-eye(prb.nx), eye(prb.nx), -eye(prb.nx), Ahatk, Bmhatk, Bphatk, whatk], ...
                                       Inf,2));


                cnstr = [cnstr;
                         zeros(prb.nx,1) == scl_mat\(- vc_plus(:,k) + vc_minus(:,k) - x(:,k+1) + Ahatk*x(:,k) + Bmhatk*u(:,k) + Bphatk*u(:,k+1) + whatk)];
            end                        
        elseif prb.disc == "ZOH"
            % Propagation
            tic
            if isfield(prb,'ode_solver')
                [Ak,Bk,wk] = feval("disc.compute_zoh_noparam_"+zoh_type,prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize,prb.ode_solver);
            else
                [Ak,Bk,wk] = feval("disc.compute_zoh_noparam_"+zoh_type,prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize);
            end
            propagate_time = toc*1000;

            for k = 1:K-1
                Ahatk   = prb.invSx*Ak(:,:,k)*prb.Sx;
                Bhatk   = prb.invSx*Bk(:,:,k)*prb.Su;
                whatk   = prb.invSx*(wk(:,k) + Ak(:,:,k)*prb.cx + Bk(:,:,k)*prb.cu - prb.cx);

                % Row normalization
                % scl_mat = eye(prb.nx);                
                scl_mat = diag(vecnorm(...
                                       [-eye(prb.nx), eye(prb.nx), -eye(prb.nx), Ahatk, Bhatk, whatk], ...
                                       Inf,2));


                cnstr = [cnstr;
                         zeros(prb.nx,1) == scl_mat\(- vc_plus(:,k) + vc_minus(:,k) - x(:,k+1) + Ahatk*x(:,k) + Bhatk*u(:,k) + whatk)];
            end
            cnstr = [cnstr; u(:,K) == u(:,K-1)];  
        elseif prb.disc == "FBP"
            % Propagation
            tic
            % In-build ODE solver is required
            [Ak,Bk,wk] = feval("disc.compute_fbp_noparam_"+fbp_type,prb.tau,xbar,ubar,prb.t_burn,prb.dyn_func,prb.dyn_func_linearize,prb.ode_solver);
            propagate_time = toc*1000;

            for k = 1:K-1
                Ahatk   = prb.invSx*Ak(:,:,k)*prb.Sx;
                Bhatk   = prb.invSx*Bk(:,:,k)*prb.Su;
                whatk   = prb.invSx*(wk(:,k) + Ak(:,:,k)*prb.cx + Bk(:,:,k)*prb.cu - prb.cx);

                % Row normalization
                % scl_mat = eye(prb.nx);                
                scl_mat = diag(vecnorm(...
                                       [-eye(prb.nx), eye(prb.nx), -eye(prb.nx), Ahatk, Bhatk, whatk], ...
                                       Inf,2));


                cnstr = [cnstr;
                         zeros(prb.nx,1) == scl_mat\(- vc_plus(:,k) + vc_minus(:,k) - x(:,k+1) + Ahatk*x(:,k) + Bhatk*u(:,k) + whatk)];
            end
            cnstr = [cnstr; u(:,K) == u(:,K-1)];            
        elseif prb.disc == "Impulse"
            % Propagation
            tic
            if isfield(prb,'ode_solver')
                [Ak,wk] = feval("disc.compute_impulse_noparam_"+impulse_type,prb.tau,xbar,ubar,prb.Eu2x,prb.h,prb.dyn_func,prb.dyn_func_linearize,prb.ode_solver);
            else
                [Ak,wk] = feval("disc.compute_impulse_noparam_"+impulse_type,prb.tau,xbar,ubar,prb.Eu2x,prb.h,prb.dyn_func,prb.dyn_func_linearize);            
            end
            propagate_time = toc*1000;

            for k = 1:K-1
                Ahatk   = prb.invSx*Ak(:,:,k)*prb.Sx;
                Bhatk   = prb.invSx*Ak(:,:,k)*prb.Eu2x*prb.Su;
                whatk   = prb.invSx*(wk(:,k) + Ak(:,:,k)*prb.cx + Ak(:,:,k)*prb.Eu2x*prb.cu - prb.cx);

                % Row normalization
                % scl_mat = eye(prb.nx);                
                scl_mat = diag(vecnorm(...
                                       [-eye(prb.nx), eye(prb.nx), -eye(prb.nx), Ahatk, Bhatk, whatk], ...
                                       Inf,2));


                cnstr = [cnstr;
                         zeros(prb.nx,1) == scl_mat\(- vc_plus(:,k) + vc_minus(:,k) - x(:,k+1) + Ahatk*x(:,k) + Bhatk*u(:,k) + whatk)];
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