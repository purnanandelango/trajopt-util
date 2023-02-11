function [H_grid,g_grid] = assemble_ltv_back_reach_polytopes(H_final,g_final,A_grid,Np,Ns)
% Assemble the backward reachable polytopes that should be avoid at each
% time of a planning horizon 1:Np+1 to ensure passive safety over a safety
% horizon 1:Ns+1
% Np : planning horizon length
% Ns : safety horizon length (the number of obstacles to avoid at each time instant is Ns+1)
    M = length(A_grid);
    assert(M>=3,'The grid of LTV matrices should have at least 3 elements.');
    assert(Np>1 && Ns>1,'Both safety and planning horizon lengths should be greater than 1.');
    assert(M>Np+Ns,'The size of the grid of LTV matrices should be at least the sum of the planning and safety horizon lengths.');
    H_grid = cell(Ns+1,Np+1);
    g_grid = cell(Ns+1,Np+1);
    for t = 1:Np+1
        [H_grid(:,t),g_grid(:,t)] = pplant.construct_back_reach_polytopes(H_final,g_final,A_grid(t:t+Ns-1));
    end
end
