function [] = setfig(varargin)
% 02/11/20
% Inspired from: https://www.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex

% fnm = 'Latin Modern Roman';
fnm = 'Palatino';
set(groot,'defaultAxesFontName',fnm,'defaultColorbarFontName',fnm,...
    'defaultGeoaxesFontName',fnm,'defaultGraphplotEdgeFontName',fnm,...
    'defaultGraphplotNodeFontName',fnm,'defaultLegendFontName',fnm,...
    'defaultPolaraxesFontName',fnm,'defaultTextFontName',fnm,...
    'defaultTextarrowshapeFontName',fnm,'defaultTextboxshapeFontName',fnm)

fsz = 42;
set(groot,'defaultAxesFontSize',fsz,'defaultColorbarFontSize',fsz,...
    'defaultGeoaxesFontSize',fsz,'defaultGraphplotEdgeFontSize',fsz,...
    'defaultGraphplotNodeFontSize',fsz,'defaultLegendFontSize',fsz,...
    'defaultPolaraxesFontSize',fsz,'defaultTextFontSize',fsz,...
    'defaultTextarrowshapeFontSize',fsz,'defaultTextboxshapeFontSize',fsz)

set(groot,'defaultAxesTitleFontWeight','normal','defaultAxesTitleFontSizeMultiplier',1,'defaultAxesLabelFontSizeMultiplier',1,'defaultAxesBox','off');

if nargin == 1
    inp = varargin{1};
else
    inp = 'latex';
end

set(groot,'defaultAxesTickLabelInterpreter',inp,'defaultColorbarTickLabelInterpreter',inp,...
    'defaultGraphplotInterpreter',inp,'defaultLegendInterpreter',inp,'defaultPolaraxesTickLabelInterpreter',inp,...
    'defaultTextInterpreter',inp,'defaultTextarrowshapeInterpreter',inp,'defaultTextboxshapeInterpreter',inp);

lw = 4;
set(groot,'defaultAxesLineWidth',lw,'defaultGraphplotLineWidth',lw,'defaultLineLineWidth',lw,...
    'defaultBarLineWidth',lw,'defaultGeoaxesLineWidth',lw,'defaultLegendLineWidth',lw,...
    'defaultLineshapeLineWidth',lw,'defaultStairLineWidth',lw);

% set(groot,'defaultAxesXGrid','on','defaultAxesXMinorGrid','on',...
%     'defaultAxesYGrid','on','defaultAxesYMinorGrid','on',...
%     'defaultAxesXMinorGridMode','manual','defaultAxesMinorGridAlpha',0.,...
%     'defaultAxesYMinorGridMode','manual');

set(groot,'defaultAxesXGrid','off','defaultAxesXMinorGrid','off',...
    'defaultAxesYGrid','off','defaultAxesYMinorGrid','off',...
    'defaultAxesXMinorGridMode','manual','defaultAxesMinorGridAlpha',0.,...
    'defaultAxesYMinorGridMode','manual');

mrksz = 45;
pos = [100,100,1000,600]; % [left bottom width height]
set(groot,'defaultFigurePosition',pos);
set(groot,'defaultLineMarkerSize',mrksz);

set(0, 'DefaultFigureRenderer', 'opengl');


end