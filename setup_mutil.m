function [] = setup_mutil(varargin)
% Add or remove utilities for trajectory optimization
    if nargin == 1
        switch varargin{1}
            case 'add'
                addpath('mutil');
                savepath
            case 'remove'
                rmpath('mutil');
                savepath                
        end
    else
        addpath('mutil');
        savepath;
    end
end