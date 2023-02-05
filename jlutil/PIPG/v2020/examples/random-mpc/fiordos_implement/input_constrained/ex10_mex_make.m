function ex10_mex_make(varargin)
%     floatprec   - precision of float data type
%                   allowed values: 'double', 'single'
%                   default: 'double'

% 
% This file is generated by FiOrdOs, a program licensed under GPL
% by copyright holder Automatic Control Laboratory, ETH Zurich.
% 
% If you are interested in using this file commercially,
% please contact the copyright holder.
% 

ip=inputParser();
ip.addParamValue('floatprec','double', @(val)(ismember(val,{'double','single'})));
ip.parse(varargin{:})
use_realtype_single = strcmp(ip.Results.floatprec,'single');

if use_realtype_single
    mex -v -DUSE_REALTYPE_SINGLE ex10_mex.c ex10_solver.c
else
    mex -v ex10_mex.c ex10_solver.c
end
end
