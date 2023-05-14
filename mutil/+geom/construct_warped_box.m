function [H,h] = construct_warped_box(box_size,box_center,varargin)
% Construct a linear inequality representation of a warped and offset cube
    if nargin >= 3
        rng(varargin{1}); % Set random number generator seed if provided
        if nargin == 4
            scl = varargin{2}; % Set scaling factor if provided
        else
            scl = 3;
        end
    end
    n = length(box_center);
    In = eye(n);
    [Q,~] = qr(randn(n)); % Construct an orthonormal basis
    T = scl*Q'*diag(min(1/scl,rand(1,n)))*Q;
    H = [ In;
         -In ]*T;
    h = H*box_center + box_size*ones(2*n,1);
end