%----------------------------- example_3 ---------------------------------%
%
% This script runs the shape changing L-MSS trust-region method with 
% four different initialization strategies
%
% This example is intended for selecting different initializations.
% You can make a selection like in lines: 100, 104, 108 (in this file):
%
% pars.whichInit = 1;
% pars.whichInit = 2;
% pars.whichInit = 3;
% pars.whichInit = 4;
%
% Details on the initializations are in 
% 
% Algorithm xxx: SC-SR1: MATLAB Software for Limited-Memory 
%   SR1 Trust-Region Methods,
% Johannes J. Brust, Oleg P. Burdakov, Jennifer B. Erway, Roummel F. Marcia
% (2022)
%
% The subproblem solver is:
%
% sc_mssm_infty.m    - Shape-changing infinity norm (pars.whichSub = 1)
%
% The objective function, f(x): R^n -> R, is the rosenbrock function. The
% dimensions of the problems are: 
%
% n in [500, 1000, 5000, 10000, 50000, 100000, 300000]
%
% This script includes the option to print the results to the
% file example_3.txt in the "DATA/" folder, with format for LaTeX
%
%-------------------------------------------------------------------------%
% 04/27/21, J.B, Preparation for release
% 06/10/22, J.B., Update of example for MSSM matrix

clc;
clear;
warning('off','MATLAB:nearlySingularMatrix');

addpath(genpath('../ALGS'));
addpath(genpath('../EXTERNAL'));
addpath(genpath('../AUXILIARY'));

printFile   = 1;
fname       = '../DATA/example_3.txt';

% Rosenbrock objective function and gradient
func = @(x)( rosen_obj(x) );
grad = @(x)( rosen_grad(x) );

% Problem dimensions
ns = [500, 1000, 5000, 10000, 50000, 100000, 300000];

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

pars.c1     = 9.e-2; % 9.e-4
pars.c2     = 0.75; % 
pars.tol    = 1e-4;
pars.print  = 0;
pars.maxiter= 200;
pars.m      = 5; % Limited memory (low memory) m = 3
pars.gammaInit = 1; % Gamma initialization gam = 10
pars.whichSub  = 1;

% Open file
if printFile == 1
    fres = fopen(fname, 'w');
end

fprintf('Example 3 ##########################################\n');    
fprintf('Rosenbrock objective: f(x)                          \n');
fprintf('n = [500, 1000, 5000, 10000, 50000, 100000, 300000] \n');
fprintf('\n');
fprintf('L-MSS-SC (Shape-changing)                           \n');
fprintf('Sub. Algorithms: TR:SC-INF                          \n');
fprintf('Initializations: Init.1, Init.2, Init.3, Init.4     \n');
fprintf('####################################################\n');
fprintf('\n');
fprintf(['n \t Init.1        	 \t Init.2        	 \t Init.3        	 \t',...
    ' Init.4         \n']);
fprintf(['- \t Time     norm(g) \t Time     norm(g) \t Time     norm(g) \t',...
    ' Time     norm(g) \n']);

% Loop over problem dimensions
for i = 1:length(ns)
    
    % Initial point
    n       = ns(i);
    x0      = zeros(n,1);
    x0(1)   = 30; % 30, 2
    
    % Using sc_mssm_infty.m (Shape-changing infinity norm solver, Alg. 3)
    pars.whichInit      = 1;
    [xk1,gk1,fk1,out1]  = LMSS_SC(x0,func,grad,pars);
    
    % Using sc_mssm_2.m (Shape-changing 2 norm solver, Alg. 4)
    % Can select "4", too for LSTRS
    pars.whichInit      = 2;
    [xk2,gk2,fk2,out2]  = LMSS_SC(x0,func,grad,pars);

    % Using obs_mssm.m (L2 norm solver)
    % Can select "5", too for Steihaug-Toint
    pars.whichInit      = 3;
    [xk3,gk3,fk3,out3]  = LMSS_SC(x0,func,grad,pars);
    
    pars.whichInit      = 4;
    [xk4,gk4,fk4,out4]  = LMSS_SC(x0,func,grad,pars);
    
    fprintf('%i \t %3.1e  %3.1e \t %3.1e  %3.1e \t %3.1e  %3.1e \t %3.1e  %3.1e \n',n,...
        out1.ctime,out1.ng,out2.ctime,out2.ng,out3.ctime,out3.ng,out4.ctime,out4.ng);
    
    if printFile == 1
        fprintf(fres, [' \\texttt{%i} & \\texttt{%i} &  \\texttt{%3.2e} & \\texttt{%3.2e} &', ...
                       ' \\texttt{%i} & \\texttt{%i} &  \\texttt{%3.2e} & \\texttt{%3.2e} &', ...
                       ' \\texttt{%i} & \\texttt{%i} &  \\texttt{%3.2e} & \\texttt{%3.2e} &', ...
                       ' \\texttt{%i} & \\texttt{%i} &  \\texttt{%3.2e} & \\texttt{%3.2e} \\\\ \n'],...
                       out1.numiter,out1.numf,out1.ctime,out1.ng,...
                       out2.numiter,out2.numf,out2.ctime,out2.ng,...
                       out3.numiter,out3.numf,out3.ctime,out3.ng,...
                       out4.numiter,out4.numf,out4.ctime,out4.ng);                   
    end    
end

% Save file
if printFile == 1
    fclose(fres);
end

warning('on','MATLAB:nearlySingularMatrix');




