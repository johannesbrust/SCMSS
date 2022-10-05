%---------------------------- EX_COMP_LMSSM_m_q5_DENSE -------------------%
%
% This script tests four L-MSS trust-region subproblem solvers
% implemented in a trust-region algorithm (LMSS_SC) on a collection of
% CUTEst problems.
% The subproblem solvers are:
%
% sc_sr1_infty.m        - Shape-changing infinity norm
% sc_sr1_2.m            - Shape-changing 2 norm
% obs.m                 - L2 norm (external solver)
% (lstrs.m               - LSTRS (external solver))[not included here]
% st_truncated_CG_qn    - Truncated CG
%
% The objective function, f(x): R^n -> R, are read from the file
% AUXILIARY/cutest_list.txt
%
% This script uses combinations of m,q and then runs the CUTEst problem.
% Here m=[3,5,7] and q = 5
%
%-------------------------------------------------------------------------%
%NOTE(S): This Experiment was carried out using CUTEst for macOS
% 04/07/21, J.B., initial version
% 08/18/22, J.B., Preparation of experiments of L-MSS
clc;
clear;

addpath(genpath('../../ALGS'));
addpath(genpath('../../EXTERNAL'));
addpath(genpath('../../AUXILIARY'));

wtest       = warning('off','all');
currentpath = pwd;

datapath    = fullfile(currentpath,'../..','/DATA/');
figpath     = fullfile(currentpath,'../..','/FIGS/');
probpath    = fullfile(currentpath,'../..','/AUXILIARY/');

rng(090317);

fprintf('---------------- EXPERIMENT CUTEst ------------------------\n');
tEX = tic; % Time experiment

% Trust-region algorithm parameters
% Detailed description of the method is in LMSS_SC.m
% Description of inputs
% x     := Initial point
% func  := Objective function; f = func(x)
% grad  := Gradient function; g = grad(x)
% pars  := Struct with parameters
%   pars.tol    := Tolerance; Stop if norm(gk,'inf') < tol
%   pars.maxiter:= Maximum iterations
%   pars.print  := Flag to pring iteration outputs
%   pars.LS     := Flag to chose line search
%   pars.LSftol := cvsrch line search parameter
%   pars.LSgtol := cvsrch line search parameter
%   pars.LSxtol := cvsrch line search parameter
%   pars.LSstpmin := cvsrch line search parameter
%   pars.LSstpmax := cvsrch line search parameter
%   pars.LSmaxfev := cvsrch line search parameter
%   pars.m          := Limited memory parameter

pars.c1     = 9.e-4; % 9.e-2
pars.c2     = 0.75; %
pars.tol    = 5e-4;
pars.print  = 0;
pars.maxiter= 5000;
pars.m      = 5; % Limited memory (low memory) m = 3
pars.gammaInit = 1.0; % Gamma initialization gam = 10
pars.storePsiPsi = 1;
pars.whichInit = 4;
pars.whichDenseInit = 4;
pars.q = 5;%floor(pars.m*1.5);
pars.SR1tolAdj = 1e-10;
% Default values of LSTRS and truncated CG

% Legend items for plotting
leg={'L-MSS:SC-INF (D)',...
    'L-MSS:SC-L2 (D)',...
    'L-MSS:L2 (D)',...
    'L-MSS:trCG (D)'};

types.colors    = ['b' 'r' 'm' 'g' 'k' 'c' 'y']; %'k' 'y'
types.lines     = {'-', '-.', ':','--','-','-.','-'}; %'-',   '-'
types.markers   = ['o' 'o' 'o' 'o' 'o' 'o' 'o']; %'s' 's'
indAlg = [1 2 3 4]; % 5 6

legLoc = 'SouthEast';
isLog = 1;


fprintf('----------- Running EX_COMP_LMSSM_m_q5 ----------- \n');
fprintf(['n \t TR:SC-INF        \t TR:SC-L2         \t TR:L2',...
    '\t tCG \n']);
fprintf(['- \t Time     norm(g) \t Time     norm(g) \t Time     norm(g) \t',...
    ' Time     norm(g) \n']);

CUTEst_init  % initialize CUTEr, see appropriate documentation

sfil    = 'TestResults';
probIdx = (1:61)';

% Initialize storage for output information
numRuns         = 1; % 3
numAlgorithms   = 4;
numProblems     = length(probIdx);
%numProblems     = 5;

pmax = 61;

% Two outer loops
ms =[3,5,7];

for i=1:length(ms)
    
    pars.m = ms(i);
    
        % Initialize results containers
        ex              = zeros(numProblems,numAlgorithms);
        numf            = zeros(numProblems,numAlgorithms);
        numg            = zeros(numProblems,numAlgorithms);
        numit           = zeros(numProblems,numAlgorithms);
        tcpu            = zeros(numProblems,numRuns,numAlgorithms);
        t_aver          = zeros(numProblems,numAlgorithms);
        tract           = zeros(numProblems,numAlgorithms);
        numrst          = zeros(numProblems,numAlgorithms);
        ngs             = zeros(numProblems,numAlgorithms);
        outs            = cell(numProblems,numAlgorithms);
        
        nms             = zeros(numProblems,2);
        
        % Reset problem counters
        p=1; % 1
        pout = 1;
        
        % read file that contains CUTEr test problems
        fid = fopen(fullfile(probpath,'cutest_list.txt'),'r');  
        tline = fgets(fid);
        while ischar(tline)
            tline = fgets(fid);
            
            if ((~strcmp(tline(1),'%'))  && (ischar(tline)) && (p <= pmax)...
                    && sum(pout == probIdx))
                
                if isunix == 1 && ismac == 0
                    eval(['!runcutest -p matlab -D ' tline]);
                else
                    tline_ = tline;
                    tline = strtrim(tline);
                    cmdcutest   = ['cutest2matlab_osx $MASTSIF/' tline];
                    unix(cmdcutest);
                end
                
                prob            = cutest_setup();
                x0              = prob.x;
                params.trradb   = max(norm(x0),1);
                n               = size(x0,1);
                clc;
                
                %mm = 10;
                %mm              = ceil(0.25*n);
                nms(p,1)        = 0;
                nms(p,2)        = n;
                
                func            = @cutest_fun;
                grad            = @cutest_gradient;
                
                %         if sum(pout == probIdxGam)
                %             pars.gammaInit = gams(pgam);
                %             pgam = pgam + 1;
                %         else
                %             pars.gammaInit = 1.0;
                %         end
                
                % Using sc_sr1_infty.m (Shape-changing infinity norm solver)  
                s=1;
                pars.whichSub       = 1;    
                pars.whichDenseInit = 4;
                [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),~,~,out1]=...
                    runAlgorithm(@LMSS_SC,func,grad,x0,pars,numRuns);
                
                ngs(p,s) = out1.ng;
                outs{p,s} = out1;
                
                % Using sc_sr1_2.m (Shape-changing 2 norm solver)
                s=s+1;
                pars.whichSub       = 2;
                pars.whichDenseInit = 4;
                [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),~,~,out2]=...
                    runAlgorithm(@LMSS_SC,func,grad,x0,pars,numRuns);
                
                ngs(p,s) = out2.ng;
                outs{p,s} = out2;
                
                % Using obs.m (L2 norm solver)
                s=s+1;
                pars.whichSub       = 3;
                pars.whichDenseInit = 1;
                [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),~,~,out3]=...
                    runAlgorithm(@LMSS_SC,func,grad,x0,pars,numRuns);
                
                ngs(p,s) = out3.ng;
                outs{p,s} = out3;
                
                % Using lstrs.m
                %         s=s+1;
                %         pars.whichSub       = 4;
                %         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),~,~,out4]=...
                %                  runAlgorithm(@LMSS_SC,func,grad,x0,pars,numRuns);
                %
                %         ngs(p,s) = out4.ng;
                %         outs{p,s} = out4;
                
                % Using st_truncated_CG_qn.m
                s=s+1;
                pars.whichSub       = 5;
                pars.whichDenseInit = 1;
                [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),~,~,out5]=...
                    runAlgorithm(@LMSS_SC,func,grad,x0,pars,numRuns);
                
                ngs(p,s) = out5.ng;
                outs{p,s} = out5;
                
                fprintf(['%i \t %3.1e  %3.1e \t %3.1e  %3.1e \t %3.1e  %3.1e \t',...
                    '%3.1e  %3.1e \n'],n,...
                    out1.ctime,out1.ng,out2.ctime,out2.ng,out3.ctime,out3.ng,...
                    out5.ctime,out5.ng);
                
                % Average CPU time
                if p==1 && numRuns > 2
                    for si=1:s
                        t_aver(p,si) = sum(tcpu(si,3:numRuns,p))/(numRuns-2);
                    end
                elseif numRuns == 2
                    for si=1:s
                        t_aver(p,si) = sum(tcpu(si,2:numRuns,p))/(numRuns-1);
                    end
                else
                    for si=1:s
                        t_aver(p,si) = tcpu(si,1:1,p);
                    end
                end
                
                cutest_terminate();
                
                p=p+1;
            end
            
            if (~strcmp(tline(1),'%'))  && (ischar(tline))
                pout = pout + 1;
            end
            
            delete( '*.d',...
                '*.o',...
                '*.dylib',...
                '*.f',...
                'mcutest.*');
            
        end
        
        % Close problem file
        fclose(fid);
        
        perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1,types,legLoc,isLog);
        box on; grid on;
        
        title(['$\textnormal{Time } (m=',num2str(pars.m),...
            ', q=',num2str(pars.q),')$'],'Interpreter','latex');
        
        % Modified printing
        fig                     = gcf;
        fig.PaperPositionMode   = 'auto';
        fig_pos                 = fig.PaperPosition;
        fig.PaperSize           = [fig_pos(3) fig_pos(4)];
                
        figname =['time_EX_LMSS_m',num2str(pars.m),...
            '_q',num2str(pars.q),'_DENSE.pdf'];       
        print(fig,'-dpdf',fullfile(figpath,figname));
        
        perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1,types,legLoc,isLog);
        box on; grid on;
        
        title(['$\textnormal{Iter } (m=',num2str(pars.m),...
            ', q=',num2str(pars.q),')$'],'Interpreter','latex');
        
        figname =['iter_EX_LMSS_m',num2str(pars.m),...
            '_q',num2str(pars.q),'_DENSE.pdf'];
        
        fig                     = gcf;
        fig.PaperPositionMode   = 'auto';
        fig_pos                 = fig.PaperPosition;
        fig.PaperSize           = [fig_pos(3) fig_pos(4)];
        
        print(fig,'-dpdf',fullfile(figpath,figname));
        
        % Function evaluations
        
        perf(ex(:,indAlg),numf(:,indAlg),leg(indAlg),1,types,legLoc,isLog);
        box on; grid on;
        
        title(['$\textnormal{Function Evaluations } (m=',num2str(pars.m),...
            ', q=',num2str(pars.q),')$'],'Interpreter','latex');
        
        % Modified printing
        fig                     = gcf;
        fig.PaperPositionMode   = 'auto';
        fig_pos                 = fig.PaperPosition;
        fig.PaperSize           = [fig_pos(3) fig_pos(4)];
                
        figname =['feval_EX_LMSS_m',num2str(pars.m),...
            '_q',num2str(pars.q),'_DENSE.pdf'];
        
        print(fig,'-dpdf',fullfile(figpath,figname));
                
        save(fullfile(datapath,['EXPERIMENT_LMSS_m',...
            num2str(pars.m),'_q',num2str(pars.q),'_DENSE']),'ex','numit',...
            't_aver','numf','numg','pars','tract', 'nms', 'outs','ngs');
        
end

delete('*.d',...
         '*.o',...
         '*.dylib',...
         '*.f',...
         'mcutest.*');

close ALL;

% Restore warning settings
warning(wtest);

tEX = toc(tEX); % Time to run experiment