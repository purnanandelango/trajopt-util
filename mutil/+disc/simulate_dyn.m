function [sol_t,sol_x,sol_u] = simulate_dyn(x0,ucell,dyn_func,tspan,N,disc_flg,varargin)
% Simulate a trajectory of nonlinear system from given initial condition and control inputs
%   [sol_t,sol_x,sol_u] = simulate_dyn(x0,ucell,dyn_func,tspan,N,disc_flg)   
%   x0 is initial condition
%   ucell = {tvec,uvec}
%   The control inputs uvec defined at discrete nodes tvec are interpolated according to flag 
%       disc_flg \in {"FOH", "ZOH", "FBP", "Impulse"}
%
%   Optional argument to specify in-built MATLAB ode solver and options
%       varargin{1} if disc_flg \in {"FOH", "ZOH"}
%       varargin{2} if disc_flg \in {"Impulse", "FBP"}
%   Required argument Eu2x: input-to-state matrix
%       varargin{1} if disc_flg = "Impulse"
%   Required argument t_burn: finite-burn pulse duration
%       varargin{1} if disc_flg = "FBP"

    tvec = ucell{1};
    uvec = ucell{2};
    [nu,K] = size(uvec);

    if disc_flg == "Impulse"
        h = (tvec(2)-tspan(1))/(N-1);        
        Eu2x = varargin{1};
        
        if K == N % Make integration step finer
            h = h/10;
        end
        zeros_nu = zeros(size(uvec,1),1);
        x_init = x0 + Eu2x*uvec(:,1);
        x = x0;
        t = tspan(1);
        for k = 1:K-1
            if nargin == 8
                [sol_t,sol_x] = feval(varargin{2}{1},@(t,x) dyn_func(t,x,zeros_nu),tvec(k:k+1),x_init,varargin{2}{2});
                sol_t = sol_t';
                sol_x = sol_x';
            else                        
                [sol_t,sol_x] = disc.rk4_march(@(t,x,u) dyn_func(t,x,u),tvec(k:k+1),x_init,h,@(t) zeros_nu);        
            end
            if k < K-1
                x_init = sol_x(:,end) + Eu2x*uvec(:,k+1);
            end
            if K ~= N
                x = [x sol_x(:,2:end)];
                t = [t sol_t(2:end)];
            else % Return open-loop propagated trajectory on the same grid that defines the impulses
                x = [x sol_x(:,end)];
                t = [t sol_t(end)];
            end
        end
        sol_t = linspace(tspan(1),tspan(2),N);
        sol_x = interp1(t',x',sol_t')';
        sol_u = nan(size(uvec,1),N);
        for k=1:K
            [~,idx] = min(abs(tvec(k)-sol_t));
            sol_u(:,idx) = uvec(:,k);
        end
    elseif disc_flg == "FOH" || disc_flg == "ZOH"
        if disc_flg == "FOH"
            u_func = @(t) disc.u_foh(t,uvec,tvec); 
        elseif disc_flg == "ZOH"
            u_func = @(t) disc.u_zoh(t,uvec,tvec);
        end
    
        if nargin == 7
            [sol_t,sol_x] = feval(varargin{1}{1},@(t,x) dyn_func(t,x,u_func(t)),linspace(tspan(1),tspan(end),N),x0,varargin{1}{2});
            sol_t = sol_t';
            sol_x = sol_x';
        else
            h = (tspan(2)-tspan(1))/(N-1);        
            [sol_t,sol_x] = disc.rk4_march(@(t,x,u) dyn_func(t,x,u),tspan,x0,h,u_func);        
        end
        
        sol_u = u_func(sol_t);
    
    elseif disc_flg == "FBP"
        t_burn = varargin{1};
        assert(nargin==8,"In-built ODE solver should be specified for FBP discretization.");
        ode_solver = varargin{2};

        x_init = x0;
        t = [];
        x = [];
        u = [];
        for k = 1:K-1
            % High-precision integration during the burn
            tspan_burn = tvec(k) + [0, t_burn];
            [sol_t,sol_x] = ode45(@(t,x) dyn_func(t,x,uvec(:,k)),tspan_burn,x_init,...
                                  odeset('AbsTol',1e-5,'RelTol',1e-4));
            M = length(sol_t);
            t = [t; sol_t(1:M-1)];
            u = [u; repmat(uvec(:,k)',[M-1,1])];
            x = [x; sol_x(1:M-1,:)];

            x_init = sol_x(M,:)';

            % Specified-precision integration for the inter-burn free-drift
            tspan_drift = [tvec(k) + t_burn, tvec(k+1)];
            [sol_t,sol_x] = feval(ode_solver{1},@(t,x) dyn_func(t,x,zeros(nu,1)),tspan_drift,x_init,...
                                  ode_solver{2});
            M = length(sol_t); 
            t = [t; sol_t(1:M-1)];
            u = [u; zeros(M-1,nu)];
            x = [x; sol_x(1:M-1,:)];   

             x_init = sol_x(M,:)';
        end
        t = [t; sol_t(end)];
        u = [u; zeros(1,nu)];
        x = [x; sol_x(end,:)];

        sol_t = linspace(tspan(1),tspan(2),N);
        sol_x = interp1(t,x,sol_t')';
        sol_u = interp1(t,u,sol_t','nearest')';        
    else
        error("Unsupported discretization flag");
    end

end