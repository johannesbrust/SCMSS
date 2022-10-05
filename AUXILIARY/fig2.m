%--------------------------- fig2.m --------------------------------------%
%
% Plot outcomes of a comparison of L-MSSM
%
% Script to replot experiments with data from the "\DATA" folder
% for each solver. This is used for comparing the best L-MSSM solver
%
% For Figure 2
%
%-------------------------------------------------------------------------%
% 03/03/21, J.B., Initial version
% 03/05/21, J.B.,
% 03/08/21, J.B., Plot memory outcomes per solver
% 03/10/21, J.B., Plots with low memory values
% 04/06/21, J.B., Preparation of plots with 4 different "m" values
% 04/08/21, J.B., Comparison Init. 2
% 04/23/21, J.B., Re-naming of file and preparation for release
% 01/18/22, J.B., Modification of script to plot additional figures
%                   as part of the revision process 
% 01/28/22, J.B., Updating of figure names after revision
% 06/09/22, J.B., Plot outcomes
% 08/26/22, J.B., Best m outcomes
% 10/05/22, J.B., Preparation for release

clc;
clear;

addpath(genpath('../ALGS'));
figpath = fullfile(pwd,'..','/FIGS/');
datapath = fullfile(pwd,'..','/DATA/');

% Memory values from experiments (with non-constant initialization)
% 62 problems used
ms = [3 3];
qs = [5 5];

hasD = [1 0];

% Additional settings based on reviewer comments

setting = 2; % 0,1,2
fignm1 = 'fig2_a';
fignm3 = 'fig2_b';

indAlg = [1 2];
np = 60; % 62
nms = length(ms);
ex = zeros(np,nms);
t_aver = zeros(np,nms);
numit = zeros(np,nms);
numf = zeros(np,nms);
nsol = 2;

% Plotting labels, and additional parameters
mleg={   'L-MSSM:SC-L2-D',...
         'L-MSSM:SC-L2 '};

%mleg = cell(nms);
types.markers   = ['o' 'o' 'o' 'o' 'o' 'o' 'o']; %'s' 's'

% Initial line types
types.colors    = ['b' 'k' 'r' 'g' 'k' 'c' 'y']; %'k' 'y'
types.lines     = {'-', '-.', '-','-.','-','-.','-'}; %'-',   '-'

leglocation = 'SouthEast';
isLog = 1;
legFontSize = 8;
scaleMaxRatio = 4.0;
typesp.colors = [];
typesp.lines = {'-',':','-.','-','--'};
typesp.markers =['o' 'o' 'o' 'o' 'o'];

for j=1:nms
    %mleg{j} = [leg{j},' ($m=',num2str(ms(j)),', q=',num2str(qs(j)),'$)'];    
%     name_EXPERIMENT = ['EXPERIMENT_FINAL_VARINIT2_',...
%         'm',num2str(ms(j)),'_','q',num2str(qs(j)),'_LMSSM'];    
    %name_EXPERIMENT = 'EXPERIMENT_FINAL_COMP_LSR1_LMSSM_SC2';
    name_EXPERIMENT = ['EXPERIMENT_LMSS_',...
        'm',num2str(ms(j)),'_','q',num2str(qs(j))];
    if hasD(j)==1
        name_EXPERIMENT =[name_EXPERIMENT,'_DENSE']; %#ok<AGROW>
    end
    
    data = load([datapath,name_EXPERIMENT],'ex','numit','t_aver','numf');
    ex(:,j) = data.ex(1:np,2);
    numit(:,j) = data.numit(1:np,2);
    t_aver(:,j) = data.t_aver(1:np,2);
    numf(:,j) = data.numf(1:np,2);
    typesp.colors = [typesp.colors,types.colors(j)];
end

% Additional plots possible, by un-commenting the next lines

% perf(ex(:,indAlg),t_aver(:,indAlg),mleg(indAlg),1,typesp,...
%     leglocation,isLog,legFontSize,scaleMaxRatio);
% box on; grid on;
% title('$\textnormal{Time (Init. 2)}$','Interpreter','latex');
% fig                     = gcf;
% fig.PaperPositionMode   = 'auto';
% fig_pos                 = fig.PaperPosition;
% fig.PaperSize           = [fig_pos(3) fig_pos(4)];
% 
% figname = ['time_I_SOLV_ALL_BEST','.pdf']; %'time_EX_CUTEST_INIT_1';
% print(fig,'-dpdf',fullfile(figpath,figname));
% 
% perf(ex(:,indAlg),numit(:,indAlg),mleg(indAlg),1,typesp,leglocation,isLog,...
%     legFontSize,scaleMaxRatio);
% box on; grid on;
% title('$\textnormal{Iter (Init. 2)}$','Interpreter','latex');
% fig                     = gcf;
% fig.PaperPositionMode   = 'auto';
% fig_pos                 = fig.PaperPosition;
% fig.PaperSize           = [fig_pos(3) fig_pos(4)];
% 
% figname = ['iter_I_SOLV_ALL_BEST','.pdf'];
% print(fig,'-dpdf',fullfile(figpath,figname));

% Extended performance profiles 
ticks   = -2; % -1
ticke   = 5; % 2
XTick   = 2.^(ticks:ticke);
XLim    = [XTick(1),XTick(end)];
leglocation1 = 'SouthEast';
legFontSize1 = 7.5;
markerSize = 0;
lineWidths = [1.8, 1.6, 1.4, 1.2];

typesp.markers =['o' 'o' 'd' 'h' 'o'];

perf_ext_fnc(ex(:,indAlg),t_aver(:,indAlg),mleg(indAlg),1,typesp,...
    leglocation1,XTick,XLim,legFontSize1,markerSize,lineWidths);

% On first plot "grid on" does not produce desired result,
% thus toggle once more
grid on; grid off;
grid on;
title(['$\textnormal{Time}$'],'Interpreter','latex');

% Modified printing
fig                     = gcf;
fig.PaperPositionMode   = 'auto';
fig_pos                 = fig.PaperPosition;
fig.PaperSize           = [fig_pos(3) fig_pos(4)];

figname = [fignm1,'.pdf']; %'time_EX_CUTEST_INIT_1';

print(fig,'-dpdf',fullfile(figpath,figname));

close ALL;

% Plot of function calls
perf_ext_fnc(ex(:,indAlg),numf(:,indAlg),mleg(indAlg),1,typesp,...
    leglocation1,XTick,XLim,legFontSize1,markerSize,lineWidths);

% On first plot "grid on" does not produce desired result,
% thus toggle once more
grid on; grid off;
grid on;
%title('$\textnormal{Function Calls}$','Interpreter','latex');
title(['$\textnormal{Function Calls}$'],'Interpreter','latex');

% Modified printing
fig                     = gcf;
fig.PaperPositionMode   = 'auto';
fig_pos                 = fig.PaperPosition;
fig.PaperSize           = [fig_pos(3) fig_pos(4)];
figname = [fignm3,'.pdf'];

print(fig,'-dpdf',fullfile(figpath,figname));

close ALL;
