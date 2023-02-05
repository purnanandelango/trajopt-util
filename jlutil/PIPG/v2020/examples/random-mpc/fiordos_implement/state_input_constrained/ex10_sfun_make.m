function ex10_sfun_make(varargin)
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
    mex -v -DUSE_REALTYPE_SINGLE ex10_sfun.c
else
    mex -v ex10_sfun.c
end


if bdIsLoaded('ex10_sfun_lib'),
    bdclose('ex10_sfun_lib');
end
if exist('ex10_sfun_lib.mdl','file')==4
    delete('ex10_sfun_lib.mdl');
end
if exist('ex10_sfun_lib.slx','file')==4
    delete('ex10_sfun_lib.slx');
end
new_system('ex10_sfun_lib','Library');

add_block('built-in/SubSystem','ex10_sfun_lib/ex10_sfun');
set_param('ex10_sfun_lib/ex10_sfun','MaskPromptString','algoInner.maxit source:| algoInner.maxit value [1x1]:|algoOuter.maxit source:| algoOuter.maxit value [1x1]:|approach.warmstartInner|algoInner.init source:| algoInner.init value [120x1]:|algoInner.stopgEps source:| algoInner.stopgEps value [1x1]:|algoInner.stopgStride [1x1]:|algoOuter.init source:| algoOuter.init value [80x1]:|algoOuter.stopgEps source:| algoOuter.stopgEps value [1x1]:|algoOuter.stopgStride [1x1]:|Sample time (>0):');
set_param('ex10_sfun_lib/ex10_sfun','MaskVariables','d_algoInner_maxit=@1;s_algoInner_maxit=@2;d_algoOuter_maxit=@3;s_algoOuter_maxit=@4;s_approach_warmstartInner=@5;d_algoInner_init=@6;s_algoInner_init=@7;d_algoInner_stopgEps=@8;s_algoInner_stopgEps=@9;s_algoInner_stopgStride=@10;d_algoOuter_init=@11;s_algoOuter_init=@12;d_algoOuter_stopgEps=@13;s_algoOuter_stopgEps=@14;s_algoOuter_stopgStride=@15;ts=@16;');
set_param('ex10_sfun_lib/ex10_sfun','MaskStyleString','popup(internal|external),edit,popup(internal|external),edit,checkbox,popup(internal|external),edit,popup(internal|external),edit,edit,popup(internal|external),edit,popup(internal|external),edit,edit,edit');
set_param('ex10_sfun_lib/ex10_sfun','MaskValueString','internal|5000|internal|5000|on|internal|[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]|internal|0.0050000000000000001|1|internal|[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]|internal|0.002|1|1');
set_param('ex10_sfun_lib/ex10_sfun','MaskTunableValueString','off,off,off,off,off,off,off,off,off,off,off,off,off,off,off,off');
set_param('ex10_sfun_lib/ex10_sfun','MaskVisibilityString','on,on,on,on,on,on,on,on,on,on,on,on,on,on,on,on');
set_param('ex10_sfun_lib/ex10_sfun','MaskTabNameString','general,general,general,general,init,init,init,algoInner.stop,algoInner.stop,algoInner.stop,init,init,algoOuter.stop,algoOuter.stop,algoOuter.stop,general');
set_param('ex10_sfun_lib/ex10_sfun','MaskCallBackString','ex10_sfun_lib_aux(''cb'',gcb,1)||ex10_sfun_lib_aux(''cb'',gcb,3)|||ex10_sfun_lib_aux(''cb'',gcb,6)||ex10_sfun_lib_aux(''cb'',gcb,8)|||ex10_sfun_lib_aux(''cb'',gcb,11)||ex10_sfun_lib_aux(''cb'',gcb,13)|||');

add_block('built-in/S-Function','ex10_sfun_lib/ex10_sfun/ex10_sfun0');
set_param('ex10_sfun_lib/ex10_sfun/ex10_sfun0','Position',[250 50 550 250]);
set_param('ex10_sfun_lib/ex10_sfun','Position',[50 50 350 250]);
set_param('ex10_sfun_lib/ex10_sfun/ex10_sfun0','Parameters','s_approach_warmstartInner,s_algoInner_stopgStride,s_algoOuter_stopgStride,ts');
set_param('ex10_sfun_lib/ex10_sfun/ex10_sfun0','FunctionName','ex10_sfun');
drawCmds={
    'port_label(''input'',1,''g [120x1]'')'
    'port_label(''input'',2,''c [1x1]'')'
    'port_label(''input'',3,''be [80x1]'')'
    'port_label(''input'',4,''\italgoInner.init [120x1]'',''texmode'',''on'')'
    'port_label(''input'',5,''\italgoInner.maxit [1x1]'',''texmode'',''on'')'
    'port_label(''input'',6,''\italgoInner.stopgEps [1x1]'',''texmode'',''on'')'
    'port_label(''input'',7,''\italgoOuter.init [80x1]'',''texmode'',''on'')'
    'port_label(''input'',8,''\italgoOuter.maxit [1x1]'',''texmode'',''on'')'
    'port_label(''input'',9,''\italgoOuter.stopgEps [1x1]'',''texmode'',''on'')'
    'port_label(''output'',1,''la [80x1]'')'
    'port_label(''output'',2,''x [120x1]'')'
    'port_label(''output'',3,''d [1x1]'')'
    'port_label(''output'',4,''iter [1x1]'')'
    'port_label(''output'',5,''exitflag [1x1]'')'
    'disp(''FiOrdOs'')'
};
set_param('ex10_sfun_lib/ex10_sfun/ex10_sfun0','MaskDisplay',char(drawCmds));

add_block('built-in/Inport','ex10_sfun_lib/ex10_sfun/g [120x1]','Position',[50 50 80 64]);
add_line('ex10_sfun_lib/ex10_sfun','g [120x1]/1','ex10_sfun0/1');
add_block('built-in/Inport','ex10_sfun_lib/ex10_sfun/c [1x1]','Position',[50 89 80 103]);
add_line('ex10_sfun_lib/ex10_sfun','c [1x1]/1','ex10_sfun0/2');
add_block('built-in/Inport','ex10_sfun_lib/ex10_sfun/be [80x1]','Position',[50 128 80 142]);
add_line('ex10_sfun_lib/ex10_sfun','be [80x1]/1','ex10_sfun0/3');

add_block('built-in/Constant','ex10_sfun_lib/ex10_sfun/algoInner.init [120x1]','Position',[50 167 80 181],'Value','s_algoInner_init');
add_line('ex10_sfun_lib/ex10_sfun','algoInner.init [120x1]/1','ex10_sfun0/4');
add_block('built-in/Constant','ex10_sfun_lib/ex10_sfun/algoInner.maxit [1x1]','Position',[50 206 80 220],'Value','s_algoInner_maxit');
add_line('ex10_sfun_lib/ex10_sfun','algoInner.maxit [1x1]/1','ex10_sfun0/5');
add_block('built-in/Constant','ex10_sfun_lib/ex10_sfun/algoInner.stopgEps [1x1]','Position',[50 245 80 259],'Value','s_algoInner_stopgEps');
add_line('ex10_sfun_lib/ex10_sfun','algoInner.stopgEps [1x1]/1','ex10_sfun0/6');
add_block('built-in/Constant','ex10_sfun_lib/ex10_sfun/algoOuter.init [80x1]','Position',[50 284 80 298],'Value','s_algoOuter_init');
add_line('ex10_sfun_lib/ex10_sfun','algoOuter.init [80x1]/1','ex10_sfun0/7');
add_block('built-in/Constant','ex10_sfun_lib/ex10_sfun/algoOuter.maxit [1x1]','Position',[50 323 80 337],'Value','s_algoOuter_maxit');
add_line('ex10_sfun_lib/ex10_sfun','algoOuter.maxit [1x1]/1','ex10_sfun0/8');
add_block('built-in/Constant','ex10_sfun_lib/ex10_sfun/algoOuter.stopgEps [1x1]','Position',[50 362 80 376],'Value','s_algoOuter_stopgEps');
add_line('ex10_sfun_lib/ex10_sfun','algoOuter.stopgEps [1x1]/1','ex10_sfun0/9');

add_block('built-in/Outport','ex10_sfun_lib/ex10_sfun/la [80x1]','Position',[700 50 730 64]);
add_line('ex10_sfun_lib/ex10_sfun','ex10_sfun0/1','la [80x1]/1');
add_block('built-in/Outport','ex10_sfun_lib/ex10_sfun/x [120x1]','Position',[700 89 730 103]);
add_line('ex10_sfun_lib/ex10_sfun','ex10_sfun0/2','x [120x1]/1');
add_block('built-in/Outport','ex10_sfun_lib/ex10_sfun/d [1x1]','Position',[700 128 730 142]);
add_line('ex10_sfun_lib/ex10_sfun','ex10_sfun0/3','d [1x1]/1');
add_block('built-in/Outport','ex10_sfun_lib/ex10_sfun/iter [1x1]','Position',[700 167 730 181]);
add_line('ex10_sfun_lib/ex10_sfun','ex10_sfun0/4','iter [1x1]/1');
add_block('built-in/Outport','ex10_sfun_lib/ex10_sfun/exitflag [1x1]','Position',[700 206 730 220]);
add_line('ex10_sfun_lib/ex10_sfun','ex10_sfun0/5','exitflag [1x1]/1');

set_param('ex10_sfun_lib/ex10_sfun','MaskSelfModifiable','on');
set_param('ex10_sfun_lib/ex10_sfun','MaskInitialization','ex10_sfun_lib_aux(''init'',gcb)');

save_system('ex10_sfun_lib','ex10_sfun_lib');
set_param('ex10_sfun_lib', 'Lock','on');
open_system('ex10_sfun_lib');

end
