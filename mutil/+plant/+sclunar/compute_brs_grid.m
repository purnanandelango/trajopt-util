% Compute the backward reachable sets for a safety horizon ts, over the planning horizon tp
function [BRS,A] = compute_brs_grid(x0,tp,ts,H,astro)

    Ns = length(ts);
    Np = length(tp);

    xbar = plant.sclunar.propagate_dyn_func_inert(x0,tp,astro,1,0);
    
    BRS = cell(Ns,Np);
    A = zeros(6,6,Np-1);
    
    for j = 1:Np
        if j<Np
            A(:,:,j) = disc.stm_ltv(tp(j),tp(j+1),xbar(:,j),@(t,x) plant.sclunar.dyn_func_inert_SRP(t,x,astro), ...
                                                            @(t,x) plant.sclunar.dyn_func_inert_SRP_jac(t,x,astro,true));
        end
        BRS{1,j} = H;
        for k = 2:Ns
            STM = disc.stm_ltv(tp(j),tp(j) + ts(k),xbar(:,j),@(t,x) plant.sclunar.dyn_func_inert_SRP(t,x,astro), ...
                                                             @(t,x) plant.sclunar.dyn_func_inert_SRP_jac(t,x,astro,true));
            BRS{k,j} = H*STM;
        end
    end

    save('brs_data','BRS','A');
    
end