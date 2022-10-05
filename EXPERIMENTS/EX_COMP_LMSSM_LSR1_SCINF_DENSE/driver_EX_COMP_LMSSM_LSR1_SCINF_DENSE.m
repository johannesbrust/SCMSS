%----------------- driver_EX_COMP_LMSSM_LSR1_m_q7_SCINF_DENSE ------------%
%
% This script tests the Shape-changing infinity norm trust-region
% solvers using the L-MSSM and L-SR1 matrix on a collection of 
% CUTEst problems.
% The subproblem solver is:
%
% sc_mssm_infty.m        - Shape-changing infinity norm
%
% The objective functions, f(x): R^n -> R, are read from the file
% AUXILIARY/cutest_list.txt
%
% The limited-memory parameter is not the same for all solvers
%-------------------------------------------------------------------------%
%NOTE(S): This Experiment was carried out using CUTEst for macOS
% 04/07/21, J.B., initial version
% 04/08/21, J.B., modifications for comparisons
% 01/20/22, J.B., including LBFGSB
% 06/09/22, J.B., including SC_LMSS
% 08/16/22, J.B., Comparison with L-SR1
% 08/18/22, J.B., Including the dense initialization
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
% Detailed description of the method is in LSR1_SC.m
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
pars.print  = 0; % 1
pars.maxiter= 50000;
pars.m      = 5; % Limited memory (low memory) m = 3
pars.gammaInit = 1.0; % Gamma initialization gam = 10
pars.storePsiPsi = 1;
pars.whichInit = 4;
pars.q = 7;%floor(pars.m*1.5);
pars.SR1tolAdj = 1e-10;
% Default values of LSTRS and truncated CG

% Parameters of L-BFGS-B
parsLB = struct('factr', 0,... % Suppress this stopping condition
                'pgtol', 5e-4,...
                'm', 5,...                
                'maxIts',25000,...
                'maxTotalIts',25000,...
                'printEvery', Inf);

% Experiment parameters            
qs =[7,5];
ms =[5,3];            
            
% Legend items for plotting
leg={'$\textnormal{L-SR1:SC-INF (D)}(m=',num2str(ms(1)),',q=',num2str(qs(1)),'$)',...
    '$\textnormal{L-MSS:SC-INF (D)}(m=',num2str(ms(1)),',q=',num2str(qs(1)),'$)',...
    '$\textnormal{L-MSS:SC-INF (D)}(m=',num2str(ms(2)),',q=',num2str(qs(2)),'$)'};

types.colors    = ['r' 'k' 'b' 'g' 'k' 'c' 'y']; %'k' 'y'
types.lines     = {'-', '-.', '-','-.','-','-.','-'}; %'-',   '-'
types.markers   = ['o' 'o' 'o' 'o' 'o' 'o' 'o']; %'s' 's'
indAlg = [1 2 3]; % 5 6
legLoc = 'SouthEast';
isLog = 1;


fprintf('----------- Running driver_EX_COMP_LMSSM_LSR1_m_q7 ----------- \n');
fprintf('n \t L-SR1            \t L-MSS(5)             \t L-MSS(3) \n');
fprintf('- \t Time     norm(g) \t Time     norm(g)  \t Time     norm(g) \n');

CUTEst_init  % initialize CUTEr, see appropriate documentation
sfil    = 'TestResults';

probIdx = (1:61)';

% Initialize storage for output information
numRuns         = 1; % 3
numAlgorithms   = 3;
numProblems     = length(probIdx);
%numProblems     = 5; % Selected set of problems

pmax = 61;%5 % 61;

for i=1:(length(ms)-1)
    
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
                
                nms(p,1)        = 0;
                nms(p,2)        = n;
                
                func            = @cutest_fun;
                grad            = @cutest_gradient;
                
                l = -inf(n,1);
                u = inf(n,1);
                parsLB.x0 = x0;
                funcLB = @(x)fminunc_wrapper( x, func, grad);
                                
                % Call to L-SR1 solver
                s=1;
                pars.whichSub       = 1;
                pars.m              = ms(1);
                pars.q              = qs(1);
                [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),~,~,out1]=...
                    runAlgorithm(@LSR1_SC,func,grad,x0,pars,numRuns);
                
                ngs(p,s) = out1.ng;
                outs{p,s} = out1;
                
                
                % Call to L-MSS, m=7
                s=s+1;
                pars.whichSub       = 1;
                pars.m              = ms(1);
                pars.q              = qs(1);
                pars.whichDenseInit = 4;
                [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),~,~,out2]=...
                    runAlgorithm(@LMSS_SC,func,grad,x0,pars,numRuns);
                
                ngs(p,s) = out2.ng;
                outs{p,s} = out2;
                
                
                % Call to L-MSS, m=3
                s=s+1;
                pars.whichSub       = 1;
                pars.m              = ms(2);
                pars.q              = qs(2);
                pars.whichDenseInit = 4;
                [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),~,~,out3]=...
                    runAlgorithm(@LMSS_SC,func,grad,x0,pars,numRuns);
                
                ngs(p,s) = out3.ng;
                outs{p,s} = out3;
                
%                 % Call to L-BFGS-B
%                 s=s+1;
%                 
%                 try
%                     tstart = tic;
%                     [xk1,fk1,out2] = lbfgsb(funcLB, l, u, parsLB);
%                     tend = toc(tstart);
% 
%                     ngLB = norm(grad(xk1),'inf');
%                     if ngLB < pars.tol
%                         ex(p,s) = 1;
%                     elseif out2.iterations == parsLB.maxIts || ...
%                             out2.totalIterations == parsLB.maxTotalIts
%                         ex(p,s) = 2;
%                     end
%                     ngs(p,s) = ngLB;
%                     outs{p,s} = out2;
%                     numit(p,s) = out2.iterations;
%                     numf(p,s) = out2.totalIterations;
%                     tcpu(s,1:1,p) = tend;
%                     %t_aver(p,s) = tend;
%                     outs{p,s} = out2;
%                 catch ME
%                     tend = -1;
%                     ngLB = -1;
%                 end
                
                fprintf('%i \t %3.1e  %3.1e \t %3.1e  %3.1e \t %3.1e  %3.1e \n',n,...
                    out1.ctime,out1.ng,out2.ctime,out2.ng,out3.ctime,out3.ng);
                
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
        
%         title(['$\textnormal{Time } (m=',num2str(pars.m),...
%             ', q=',num2str(pars.q),')$'],'Interpreter','latex');
        
        title('$\textnormal{Time} $','Interpreter','latex');
        
        % Modified printing
        fig                     = gcf;
        fig.PaperPositionMode   = 'auto';
        fig_pos                 = fig.PaperPosition;
        fig.PaperSize           = [fig_pos(3) fig_pos(4)];
                
        figname = 'time_EX_COMP_LSR1_LMSSM_SCINF_DENSE.pdf';
        
        print(fig,'-dpdf',fullfile(figpath,figname));
        
        perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1,types,legLoc,isLog);
        box on; grid on;
        
        title('$\textnormal{Iter} $','Interpreter','latex');
        
%         title(['$\textnormal{Iter } (m=',num2str(pars.m),...
%             ', q=',num2str(pars.q),')$'],'Interpreter','latex');

        figname = 'iter_EX_COMP_LSR1_LMSSM_SCINF_DENSE.pdf';
        
%         figname =['iter_EX_FINAL_VARINIT2_m',num2str(pars.m),...
%             '_q',num2str(pars.q),'_LMSSM.pdf'];
        
        fig                     = gcf;
        fig.PaperPositionMode   = 'auto';
        fig_pos                 = fig.PaperPosition;
        fig.PaperSize           = [fig_pos(3) fig_pos(4)];
        
        print(fig,'-dpdf',fullfile(figpath,figname));
        
        % Function evaluations
        
        perf(ex(:,indAlg),numf(:,indAlg),leg(indAlg),1,types,legLoc,isLog);
        box on; grid on;

        title('$\textnormal{Function Evaluations} $','Interpreter','latex');
        
%         title(['$\textnormal{Function Evaluations } (m=',num2str(pars.m),...
%             ', q=',num2str(pars.q),')$'],'Interpreter','latex');
        
        % Modified printing
        fig                     = gcf;
        fig.PaperPositionMode   = 'auto';
        fig_pos                 = fig.PaperPosition;
        fig.PaperSize           = [fig_pos(3) fig_pos(4)];
                
%         figname =['feval_EX_FINAL_VARINIT2_m',num2str(pars.m),...
%             '_q',num2str(pars.q),'_LMSSM.pdf'];

        figname = 'feval_EX_COMP_LSR1_LMSSM_SCINF_DENSE.pdf';
        
        print(fig,'-dpdf',fullfile(figpath,figname));
        
        
                
        save(fullfile(datapath,'EXPERIMENT_FINAL_COMP_LSR1_LMSSM_SCINF_DENSE'),'ex','numit',...
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