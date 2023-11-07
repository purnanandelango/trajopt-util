function [xbar,ubar,pbar] = create_initialization(prb,flg,varargin)
    switch flg
        case 1 % Straight initialization
            
            K = prb.K;
            xbar = zeros(prb.nx,K);
            ubar = zeros(prb.nu,K);
            if isfield(prb,'p')
                pbar = prb.p;
            else
                pbar = [];
            end

            for k = 1:K
                xbar(:,k) = ( (prb.tau(end)-prb.tau(k))*prb.x1 + (prb.tau(k)-prb.tau(1))*prb.xK )/(prb.tau(end)-prb.tau(1)); 
                ubar(:,k) = ( (prb.tau(end)-prb.tau(k))*prb.u1 + (prb.tau(k)-prb.tau(1))*prb.uK )/(prb.tau(end)-prb.tau(1));
            end

            % If the indices of the quaternions are known, then perform SLERP using q_lib/q_slerp
            if isfield(prb,'quaternion_idx')
                for k = 1:K
                    xbar(prb.quaternion_idx,k) = qlib.q_slerp(prb.x1(prb.quaternion_idx),prb.xK(prb.quaternion_idx),prb.tau(k));
                end
            end

        case 2 % Warmstart by providing previous solution as guess

            x = varargin{1};
            u = varargin{2};
            if isempty(varargin{3})
                tau = linspace(prb.tau(1),prb.tau(end),size(x,2));
            else
                tau = varargin{3};
            end

            xbar = interp1(tau,x',prb.tau')';
            ubar = interp1(tau,u',prb.tau')';
            
            if nargin == 6
                pbar = varargin{4};
            else
                pbar = [];
            end

        otherwise

            error('Invalid flag for initialization.')
    end
end