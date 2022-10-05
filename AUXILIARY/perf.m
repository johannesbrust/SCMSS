function perf(ex,T,legstr, varargin)
%PERF    Performace profiles
%
% PERF(T,logplot)-- produces a performace profile as described in
%   Benchmarking optimization software with performance profiles,
%   E.D. Dolan and J.J. More', 
%   Mathematical Programming, 91 (2002), 201--213.
% Each column of the matrix T defines the performance data for a solver.
% Failures on a given problem are represented by a NaN.
%
% This function is based on the perl script of Liz Dolan.
%
% Jorge J. More', June 2004
%
% Modified by Spartak Zikrin, November 2015

% 12/11/16, J.B. extension for plotting on beamer slides.
% 07/12/17, J.B. extension for passing line types.
% 05/02/19, J.B. extension for passing legend location

[np,ns] = size(T);


if nargin < 3
    leg=1:ns;
    leg=leg';
    legstr=num2str(leg);
end

isbeam = 0;
if nargin > 3
    
    isbeam = varargin{1};
    
end

colors  = ['b' 'r' 'k' 'm' 'c' 'g' 'y'];   
lines   = {'-' '--' '-.'};
markers = [ 's' 'o' '^' 'v' 'p' '<' 'x' 'h' '+' 'd' '*' '<' ];

leglocation = 'SouthEast';

hastypes = 0;
if nargin > 4
   
    hastypes = 1;
    types   = varargin{2};
    colors  = types.colors;
    lines   = types.lines;
    markers = types.markers;
    
end
if nargin > 5
    
    leglocation = varargin{3};
    
end

isLog = 0;
if nargin > 6
    
    isLog = varargin{4};
    
end

if ns > 3
    legFontSize = 10;
else
    legFontSize = 14;
end
if nargin > 7
    
    legFontSize = varargin{5};
    
end
scaleMaxRatio = 1.05;
if nargin > 8
    
    scaleMaxRatio = varargin{6};
    
end


T(ex~=1)=NaN;

% Minimal performance per solver
minperf = min(T,[],2);

% Compute ratios and divide by smallest element in each row.
r = zeros(np,ns);
for p = 1: np
    r(p,:) = T(p,:)/minperf(p);
end

max_ratio = max(max(r));

% Replace all NaN's with twice the max_ratio and sort.
r(isnan(r)) = 2*max_ratio;
r = sort(r);

% if isLog == 1
%     r = log(r);
% end

% Plot stair graphs with markers.
% Set figure properties
if isbeam == 1
    x0 = 0; y0 = 0;
    width =4; height =3;

    figure('Units','inches',...
        'Position',[x0 y0 width height],...
        'PaperPositionMode','auto');
else
    figure;
end

hold on;

% legend(legstr,'Location','SouthEast');
ymax=1;
ymin=1;
for s = 1: ns
    [xs,ys] = stairs(r(:,s),[1:np]/np);
    if ys(end)==1
         ys(end+1)=1;
         ymax=1.01;
         xs(end+1)=1.05*max_ratio;
    end
    
    if hastypes == 0
        sl = mod(s-1,3) + 1; 
        sc = mod(s-1,7) + 1;
    else
        sl = s;
        sc = s;
    end
    
    option = [char(lines(sl)) colors(sc)]; 

    i1=find(xs==1,1,'last');
    if isempty(i1)
     i1=1;
    end
    xs=xs(i1:end); 
    ys=ys(i1:end);
    ymin=max(0,min(ymin, ys(1)-0.05));
    if isLog == 1
        semilogx(xs,ys,option,'LineWidth',2);
    else
        plot(xs,ys,option,'LineWidth',2); %2
    end
    %
end

for s = 1: ns
    [xs,ys] = stairs(r(:,s),[1:np]/np);
    i1=find(xs==1,1,'last');
    if isempty(i1), i1=1; end
    xs=xs(i1:end); 
    ys=ys(i1:end);
    if hastypes == 0
        sl = mod(s-1,3) + 1; 
        sc = mod(s-1,7) + 1;
    else
        sl = s;
        sc = s;
    end
    options0=[char(lines(sl))  colors(sc) markers(2)]; % markers(sm)
    if isLog == 1
        semilogx(xs(1),ys(1),options0,'MarkerFaceColor',colors(sc),'MarkerSize',8);
    else
        plot(xs(1),ys(1),options0,'MarkerFaceColor',colors(sc),'MarkerSize',8);
    end
end

% Axis properties are set so that failures are not shown,
% but with the max_ratio data points shown. This highlights
% the "flatline" effect.

if isLog == 1
    smax_ratio = scaleMaxRatio*max_ratio; % scaled max ratio
    if log(smax_ratio) < 1
        axis([ 1 smax_ratio 0 1 ]);
    else
        axis([ 1 log(smax_ratio) 0 1 ]);
    end
else
    axis([ 1 scaleMaxRatio*max_ratio 0 1 ]);
end

% Legends and title should be added.
if ns > 3
    %legend(legstr,'Location','SouthEast','Fontsize',7);
    %legend(legstr,'Location','SouthEast','Fontsize',7,'interpreter','latex');
    
%elseif ns == 3
    legend(legstr,'Location',leglocation,'Fontsize',legFontSize,'interpreter','latex'); % SouthEast
else
    legend(legstr,'Location',leglocation,'Fontsize',legFontSize,'interpreter','latex');
end
xmin=1;
ymin=max(0,ymin);

if isLog == 1
    smax_ratio = scaleMaxRatio*max_ratio;
    if log(smax_ratio) < 1
        axis([ 1 smax_ratio 0 1 ]);
    else
        axis([ 1 log(smax_ratio) 0 1 ]);
    end
    %axis([ xmin log(1.05*max_ratio) ymin ymax]);
    xlabel('$\tau$','FontSize',14,'interpreter','latex');
else
	axis([xmin scaleMaxRatio*max_ratio ymin ymax]);
    xlabel('$\tau$','FontSize',14,'interpreter','latex');
end

ylabel('$\rho_s(\tau)$','FontSize',14,'interpreter','latex');

