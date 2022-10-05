%----------------------------- example_2 ---------------------------------%
%
% This script runs three L-MSS trust-region methods. 
% The subproblem solvers are two shape-changing norms 
% and the external solver "obs_mssm.m" (L2 norm) 
%
% This example is intended for selecting different subproblems solvers.
% You can make a selection like in lines: 100, 104, 108 (in this file):
%
% pars.whichSub = 1;
% pars.whichSub = 2;
% pars.whichSub = 3;
%
% The subproblem solvers are implemented in:
%
% sc_mssm_infty.m    - Shape-changing infinity norm (pars.whichSub = 1)
% sc_mssm_2.m        - Shape-changing 2 norm (pars.whichSub = 2)
% obs_mssm.m         - L2 norm (external solver) (pars.whichSub = 3)
%
% In addition two more solvers are implemented with
%
% pars.whichSub = 3; - LSTR
% pars.whichSub = 4; - Steihaug-Toint
%
% The objective function, f(x): R^n -> R, is the rosenbrock function. The
% dimensions of the problems are: 
%
% n in [500, 1000, 5000, 10000, 50000, 100000, 300000]
%
% This script includes the option to print the results to the
% file example_2.txt in the "DATA/" folder, with format for LaTeX
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
fname       = '../DATA/example_2.txt';

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
pars.whichInit = 4;

% Open file
if printFile == 1
    fres = fopen(fname, 'w');
end

fprintf('Example 2 ##########################################\n');    
fprintf('Rosenbrock objective: f(x)                          \n');
fprintf('n = [500, 1000, 5000, 10000, 50000, 100000, 300000] \n');
fprintf('\n');
fprintf('L-MSS-SC (Shape-changing)                          \n');
fprintf('Sub. Algorithms: TR:SC-INF, TR:SC-L2, TR:L2        \n');
fprintf('####################################################\n');
fprintf('\n');
fprintf('n \t TR:SC-INF        \t TR:SC-L2         \t TR:L2 \n');
fprintf('- \t Time     norm(g) \t Time     norm(g) \t Time     norm(g) \n');

% Loop over problem dimensions
for i = 1:length(ns)
    
    % Initial point
    n       = ns(i);
    x0      = zeros(n,1);
    x0(1)   = 30; % 30, 2
    
    % Using sc_mssm_infty.m (Shape-changing infinity norm solver, Alg. 3)
    pars.whichSub       = 1;
    [xk1,gk1,fk1,out1]  = LMSS_SC(x0,func,grad,pars);
    
    % Using sc_mssm_2.m (Shape-changing 2 norm solver, Alg. 4)
    % Can select "4", too for LSTRS
    pars.whichSub       = 2; % 4
    [xk2,gk2,fk2,out2]  = LMSS_SC(x0,func,grad,pars);

    % Using obs_mssm.m (L2 norm solver)
    % Can select "5", too for Steihaug-Toint
    pars.whichSub       = 3; % 5
    [xk3,gk3,fk3,out3]  = LMSS_SC(x0,func,grad,pars);
    
    fprintf('%i \t %3.1e  %3.1e \t %3.1e  %3.1e \t %3.1e  %3.1e \n',n,...
        out1.ctime,out1.ng,out2.ctime,out2.ng,out3.ctime,out3.ng);
    
    if printFile == 1
        fprintf(fres, [' \\texttt{%i} & \\texttt{%i} &  \\texttt{%3.2e} & \\texttt{%3.2e} &', ...
                       ' \\texttt{%i} & \\texttt{%i} &  \\texttt{%3.2e} & \\texttt{%3.2e} &', ...
                       ' \\texttt{%i} & \\texttt{%i} &  \\texttt{%3.2e} & \\texttt{%3.2e} \\\\ \n'],...
                       out1.numiter,out1.numf,out1.ctime,out1.ng,...
                       out2.numiter,out2.numf,out2.ctime,out2.ng,...
                       out3.numiter,out3.numf,out3.ctime,out3.ng);
    end    
end

% Save file
if printFile == 1
    fclose(fres);
end

warning('on','MATLAB:nearlySingularMatrix');




