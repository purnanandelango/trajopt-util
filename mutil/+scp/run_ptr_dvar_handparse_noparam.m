function [xbar,ubar,cost_val,converged] = run_ptr_dvar_handparse_noparam(xbar,ubar,prb)
% PTR SCP without parameters as decision variables and ZOH/FOH discretization
% Subproblem solver input is hand-parsed
% No provision for updating the problem parameters after each SCP iteration
% Exact penalty weight cannot be matrix-valued; only scalar wvc is allowed

    converged = false;
    K = prb.K;
    nx = prb.nx;
    nu = prb.nu;

    % Check if type of FOH computation is specified
    if isfield(prb,'foh_type')
        foh_type = string(prb.foh_type);
        assert(ismember(foh_type,["v1","v2","v3","v3_parallel"]),"Incorrect type of FOH discretization.");        
    else
        foh_type = "v3"; % Default
    end

    ni = size(prb.Ei,1);
    nf = size(prb.Ef,1);
    ny = size(prb.Ey,1);
    nz = (nx+nu)*K+2*nx*(K-1);    

    zhat_i          = (prb.Ei*prb.invSx*prb.Ei')*(prb.zi-prb.Ei*xbar(:,1));
    zhat_f          = (prb.Ef*prb.invSx*prb.Ef')*(prb.zf-prb.Ef*xbar(:,K)); 
    uhatmax         = zeros(nu*K,1);
    uhatmin         = zeros(nu*K,1);
    eps_hat         = zeros(ny*(K-1),1);
    for k = 1:K
        uhatmin((k-1)*nu+1:k*nu) = prb.invSu*(prb.umin-ubar(:,k));
        uhatmax((k-1)*nu+1:k*nu) = prb.invSu*(prb.umax-ubar(:,k));
        if k < K
            eps_hat((k-1)*ny+1:k*ny) = prb.eps_cnstr*prb.Ey*prb.invSx*prb.Ey'*ones(ny,1) - prb.Ey*prb.invSx*(xbar(:,k+1)-xbar(:,k));
        end
    end

    % Parse to QP canonical form
    Phat            = blkdiag(prb.wtr*speye((nx+nu)*K),sparse(2*nx*(K-1),2*nx*(K-1)));
    phat            = prb.cost_factor*[sparse(nx*(K-1),1); 
                                       prb.term_cost_vec;
                                       sparse(nu*K+2*nx*(K-1),1)] ...
                      + prb.wvc*[sparse((nx+nu)*K,1);
                                 ones(2*nx*(K-1),1)];
    Ghat_if         = [prb.Ei              sparse(ni,nx*(K-1));
                       sparse(nf,nx*(K-1)) prb.Ef];
    Ghat_if         = [Ghat_if sparse((ni+nf),nu*K+2*nx*(K-1))];
    Ghat_mu         = [speye(nx*(K-1)) -speye(nx*(K-1))];
    if prb.disc == "ZOH"
        ghat            = [zhat_i;
                           zhat_f;
                           zeros(nx*(K-1)+nu,1)];
        n_eq_cnstr      = ni + nf + nx*(K-1) + nu;
    else
        ghat            = [zhat_i;
                           zhat_f;
                           zeros(nx*(K-1),1)];        
        n_eq_cnstr      = ni + nf + nx*(K-1);
    end
    Hhat_u_mu       = [sparse(2*nu*K,nx*K)     [speye(nu*K); -speye(nu*K)] sparse(2*nu*K,2*nx*(K-1));
                       sparse(2*nx*(K-1),(nx+nu)*K) -speye(2*nx*(K-1))];
    Hhat_y          = - [kron(speye(K-1),prb.Ey) sparse(ny*(K-1),nx)] + [sparse(ny*(K-1),nx) kron(speye(K-1),prb.Ey)]; 
    Hhat            = [Hhat_u_mu;
                       Hhat_y sparse(ny*(K-1),nu*K+2*nx*(K-1))];
    % hhat            = [ uhatmax;
    %                    -uhatmin;
    %                     sparse(2*nx*(K-1),1);
    %                     eps_hat];
    n_ineq_cnstr    = 2*nu*K + 2*nx*(K-1) + ny*(K-1);

    % PIPG
    Htil            = [Hhat_y sparse(ny*(K-1),nu*K+2*nx*(K-1))];  
    % htil            = eps_hat;
    % Hhtil           = [Htil htil];
    % Hhtil_normal    = linalg.mat_normalize(Hhtil,'row');

    % Scale reference state and control
    zbar = zeros(nz,1);
    
    fprintf("+---------------------------------------------------------------------------------------+\n");
    fprintf("|                            ..:: Penalized Trust Region ::..                           |\n");
    fprintf("+-------+------------+-----------+-----------+---------+---------+------------+---------+\n");
    fprintf("| Iter. | Prop. [ms] | Prs. [ms] | Slv. [ms] | log(TR) | log(VC) |    Cost    |   ToF   |\n");
    fprintf("+-------+------------+-----------+-----------+---------+---------+------------+---------+\n");
    
    for j = 1:prb.scp_iters

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
            
            tic
            Ahat = [];
            Bmhat = [];
            Bphat = [];
            what = [];
            for k = 1:prb.K-1
                Ahat = blkdiag(Ahat,prb.invSx*Ak(:,:,k)*prb.Sx);                
                Bmhat = blkdiag(Bmhat,prb.invSx*Bmk(:,:,k)*prb.Su);                
                Bphat = blkdiag(Bphat,prb.invSx*Bpk(:,:,k)*prb.Su);
                what = [what;
                        prb.invSx*(xbarprop(:,k+1)-xbar(:,k+1))];
            end        
            
            Ghat_x = [Ahat  sparse(nx*(K-1),nx)] - [sparse(nx*(K-1),nx) speye(nx*(K-1))];
            Ghat_u = [Bmhat sparse(nx*(K-1),nu)] + [sparse(nx*(K-1),nu) Bphat];
            Ghat = [Ghat_if;
                    Ghat_x Ghat_u Ghat_mu];
            ghat((ni+nf)+1:end) = -what;  
            % phat = phat_cost_vc;
            parse_time = toc*1000;

        elseif prb.disc == "ZOH"
            % Propagation
            tic
            if isfield(prb,'ode_solver')
                [Ak,Bk,~,~,xbarprop] = disc.compute_zoh_noparam(prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize,prb.ode_solver);
            else
                [Ak,Bk,~,~,xbarprop] = disc.compute_zoh_noparam(prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize);
            end
            propagate_time = toc*1000;
            
            tic
            Ahat = [];
            Bhat = [];
            what = [];            
            for k = 1:prb.K-1
                Ahat = blkdiag(Ahat,prb.invSx*Ak(:,:,k)*prb.Sx);                
                Bhat = blkdiag(Bhat,prb.invSx*Bk(:,:,k)*prb.Su);
                what = [what;
                        prb.invSx*(xbarprop(:,k+1)-xbar(:,k+1))];                
            end

            Ghat_x = [Ahat sparse(nx*(K-1),nx)] - [sparse(nx*(K-1),nx) speye(nx*(K-1))];
            Ghat_u = [Bhat sparse(nx*(K-1),nu)]; 
            Ghat = [Ghat_if;
                    Ghat_x Ghat_u Ghat_mu;
                    sparse(nu,nx*K+nu*(K-2)) speye(nu) -speye(nu) sparse(nu,2*nx*(K-1))];
            ghat((ni+nf)+1:end) = [-what;sparse(nu,1)];  
            % phat = phat_cost_vc;
            parse_time = toc*1000;

        end

        zhat_i          = (prb.Ei*prb.invSx*prb.Ei')*(prb.zi-prb.Ei*xbar(:,1));
        zhat_f          = (prb.Ef*prb.invSx*prb.Ef')*(prb.zf-prb.Ef*xbar(:,K)); 
        for k = 1:K
            uhatmin((k-1)*nu+1:k*nu) = prb.invSu*(prb.umin-ubar(:,k));
            uhatmax((k-1)*nu+1:k*nu) = prb.invSu*(prb.umax-ubar(:,k));
            if k < K
                eps_hat((k-1)*ny+1:k*ny) = prb.eps_cnstr*prb.Ey*prb.invSx*prb.Ey'*ones(ny,1) - prb.Ey*prb.invSx*(xbar(:,k+1)-xbar(:,k));
            end
        end   

        ghat(1:ni+nf)   = [zhat_i;zhat_f];
        hhat            = [ uhatmax;
                           -uhatmin; 
                            sparse(2*nx*(K-1),1);
                            eps_hat]; 

        tic
        if prb.solver.name == "quadprog"
            z = quadprog(Phat,phat,Hhat,hhat,Ghat,ghat,[],[],[],optimoptions('quadprog','Algorithm','interior-point-convex',... 'active-set',...
                         'Display',prb.solver.Display,'ConstraintTolerance',prb.solver.ConstraintTolerance,'OptimalityTolerance',prb.solver.OptimalityTolerance));
        elseif prb.solver.name == "ecos"
            z = ecosqp(Phat,phat,Hhat,full(hhat),Ghat,full(ghat),ecosoptimset('verbose',prb.solver.verbose,'abstol',prb.solver.abstol,'reltol',prb.solver.reltol));
        elseif prb.solver.name == "piqp"
            solver = piqp('sparse');
            solver.update_settings('compute_timings',false,'verbose',prb.solver.verbose,...
                                   'eps_abs',prb.solver.eps_abs,'eps_rel',prb.solver.eps_rel,...
                                   'eps_duality_gap_abs',prb.solver.eps_duality_gap_abs,'eps_duality_gap_rel',prb.solver.eps_duality_gap_rel);
            solver.setup(Phat, phat, Ghat, ghat, Hhat, hhat, [], []);
            result = solver.solve();
            z = result.x;
        elseif prb.solver.name == "scs"
            cones = struct;
            cones.z = n_eq_cnstr;
            cones.l = n_ineq_cnstr;
            data = struct;
            data.P = Phat;
            data.c = full(phat);
            data.A = [Ghat;Hhat];
            data.b = full([ghat;hhat]);
            if j > 1
                data.x = z;
                data.y = y_scs;
                data.s = s_scs;
            end
            settings = struct('eps_abs',prb.solver.eps_abs,'eps_rel',prb.solver.eps_rel,'verbose',prb.solver.verbose);
            [z,y_scs,s_scs] = scs(data,cones,settings);
        elseif prb.solver.name == "mosek"
            blc = [ghat;-Inf(n_ineq_cnstr,1)];
            buc = [ghat;hhat];
            settings = struct('MSK_DPAR_INTPNT_QO_TOL_PFEAS',prb.solver.MSK_DPAR_INTPNT_QO_TOL_PFEAS,'MSK_DPAR_INTPNT_QO_TOL_DFEAS',prb.solver.MSK_DPAR_INTPNT_QO_TOL_DFEAS,'MSK_DPAR_INTPNT_QO_TOL_REL_GAP',prb.solver.MSK_DPAR_INTPNT_QO_TOL_REL_GAP);
            res = mskqpopt(Phat,phat,[Ghat;Hhat],blc,buc,[],[],settings,'minimize info echo(0)');
            z = res.sol.itr.xx;
        elseif prb.solver.name == "osqp"
            solver = osqp;
            settings = struct('eps_abs',prb.solver.eps_abs,'eps_rel',prb.solver.eps_rel,'max_iter',prb.solver.max_iter,'verbose',prb.solver.verbose);
            solver.setup(Phat,phat,[Ghat;Hhat],[ghat;-Inf(2*nu*K+2*nx*(K-1)+ny*(K-1),1)],[ghat;hhat],settings);
            if j > 1
                solver.warm_start('x',z,'y',y_osqp);
            end            
            res = solver.solve();
            z = res.x; 
            y_osqp = res.y;
        elseif prb.solver.name == "gurobi"
            model = struct;
            model.Q = Phat/2;
            model.obj = full(phat);
            model.A = [Ghat;Hhat];
            model.rhs = full([ghat;hhat]);
            model.lb = -Inf((nx+nu)*K+2*nx*(K-1),1);
            % model.ub = [];
            model.sense = [repmat('=',[n_eq_cnstr,1]);
                           repmat('<',[n_ineq_cnstr,1])];
            params = struct;
            params.OptimalityTol = prb.solver.OptimalityTol;
            params.FeasibilityTol = prb.solver.FeasibilityTol;
            params.OutputFlag = prb.solver.verbose;
            result = gurobi(model,params);
            z = result.x;
        elseif prb.solver.name == "pipg"

            htil         = eps_hat;
            Hhtil        = [Htil htil];
            Hhtil_normal = linalg.mat_normalize(Hhtil,'row');  
            % Hhtil_normal = Hhtil;
            
            Ggtil        = [Ghat(ni+nf+1:end,:) ghat(ni+nf+1:end)];
            Ggtil_normal = linalg.mat_normalize(Ggtil,'row');
            % Ggtil_normal = Ggtil;

            model = struct;
            model.nx = nx;
            model.nu = nu;
            model.K = K;
            model.Phat = prb.solver.lambda*Phat;
            model.phat = prb.solver.lambda*phat;
            model.Gtil = Ggtil_normal(:,1:end-1);
            model.gtil = Ggtil_normal(:,end);
            model.Htil = Hhtil_normal(:,1:end-1);
            model.htil = Hhtil_normal(:,end);
            model.uhatmin = uhatmin;
            model.uhatmax = uhatmax;
            model.i_idx = prb.i_idx;
            model.f_idx = prb.f_idx;
            model.zhat_i = ghat(1:ni);
            model.zhat_f = ghat(ni+1:ni+nf);

            sGtilHtil = svd(full([model.Gtil;model.Htil]));

            options = struct;
            options.alpha = 2/(sqrt(prb.wtr^2 + 4*prb.solver.omega*(sGtilHtil(1)^2))+prb.wtr);
            options.beta = prb.solver.omega*options.alpha;
            options.rho = prb.solver.rho;
            options.eps_abs = prb.solver.eps_abs;
            options.max_iter = prb.solver.max_iter;
            options.verbose = prb.solver.verbose;
            options.test_termination = prb.solver.test_termination;

            if j > 1
                [z,eta_pipg,chi_pipg] = solvers.pipg(model,options,zbar,eta_pipg,chi_pipg);
            else
                [z,eta_pipg,chi_pipg] = solvers.pipg(model,options,zbar);
            end
        end
        solve_time = toc*1000;

        x_scl    = reshape(z(1:nx*K),[nx,K]);
        u_scl    = reshape(z(nx*K+1:(nx+nu)*K),[nu,K]);
        x        = prb.Sx*x_scl + xbar;
        u        = prb.Su*u_scl + ubar;
        cost_val = prb.term_cost_vec'*prb.invSx*x(:,K);
        vc_term  = [sparse(1,(nx+nu)*K), ones(1,2*nx*(K-1))]*z;

        % Ensure that the TR value is always displayed consistently with infinity norm
        % Note that this is for display and for evaluation of termination criteria 
        Jtr_post_solve = zeros(1,K);        
        for k = 1:K
            Jtr_post_solve(k) = norm([x_scl(:,k);u_scl(:,k)],'inf');
        end
        tr_term = max(Jtr_post_solve);

        % Update reference trajectory
        xbar     = x;
        ubar     = u;
        zbar     = z;

        ToF = prb.time_of_maneuver(xbar,ubar);        
        
        % Console output
        fprintf('|  %02d   |   %7.1e  |  %7.1e  |  %7.1e  | %5.1f   | %5.1f   | %10.3e | %7.1e |\n',j,propagate_time,parse_time,solve_time,log10(tr_term),log10(vc_term),cost_val,ToF);
        
        if vc_term < prb.epsvc && tr_term < prb.epstr
            converged = true;
            fprintf("+---------------------------------------------------------------------------------------+\n");
            fprintf('Converged!\n')
            break
        end
        
    end

end