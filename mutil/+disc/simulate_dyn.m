function [sol_t,sol_x,sol_u] = simulate_dyn(x0,ucell,dyn_func,tspan,N,disc_flg)
% Simulate a trajectory of nonlinear system from given initial condition and control inputs
%   [sol_t,sol_x,sol_u] = simulate_dyn(x0,ucell,dyn_func,tspan,N,disc_flg)   
%   x0 is initial condition
%   ucell = {tvec,uvec}
%   The control inputs uvec defined at discrete nodes tvec are interpolated
%   according to flag disc_flg \in {"FOH","ZOH"}
    tvec = ucell{1};
    uvec = ucell{2};

    if disc_flg == "FOH"
        u_func = @(t) disc.u_foh(t,uvec,tvec); 
    elseif disc_flg == "ZOH"
        u_func = @(t) disc.u_zoh(t,uvec,tvec);
    end

    h = (tspan(2)-tspan(1))/(N-1);

    [sol_t,sol_x] = disc.rk4_march(@(t,x,u) dyn_func(t,x,u),tspan,x0,h,u_func);
    sol_u = u_func(sol_t);
end