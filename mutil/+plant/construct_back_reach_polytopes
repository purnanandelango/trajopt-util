function [H_grid,g_grid] = construct_back_reach_polytopes(H_final,g_final,A_grid,varargin)
% Compute the sequence of backward reachable polytopes for a discrete-time
% LTI system defined by a single matrix A
% LTV system defined by a grid of matices A_grid
% Representation of final polytope
% { x | H_final*x <= g_final }

    if iscell(A_grid)
        N = length(A_grid);             % Horizon length
        H_grid = cell(1,N+1);           % Containers for holding BRS polytopes
        g_grid = cell(1,N+1);           
        H_grid{1} = H_final;
        g_grid{1} = g_final;
        if N>0
            STM = eye(size(H_final,2));
            for t = 1:N
               STM = A_grid{t}*STM;
               H_grid{t+1} = H_final*STM;
               g_grid{t+1} = g_final;
            end
        end
    else
        if nargin ~= 4
            error('Need three arguments when second argument is not a cell.');
        else
            A = A_grid;
            N = varargin{1};
            H_grid = cell(1,N+1);           % Containers for holding BRS polytopes
            g_grid = cell(1,N+1);           
            H_grid{N+1} = H_final;
            g_grid{N+1} = g_final;
            STM = eye(size(H_final,2));            
            for t = N:-1:1
                STM = STM*A;
                H_grid{t} = H_final*STM;
                g_grid{t} = g_final;
            end
        end
    end
end
