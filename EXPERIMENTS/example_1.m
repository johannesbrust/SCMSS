%--------------------------- example_1 -----------------------------------%
%
% This script runs LMSS_SC (a limited-memory Shape-Changing Multipoint-
% Symmetric Secant method 
%
% LMSS-SC: MATLAB Software for Limited-Memory MSS Trust-Region Methods,
% (Brust,Erway,Marcia,2022)                               
%
% The algorithm can be used for general unconstrained minimization problems
%
%                           minimize f(x),
%
% where f:R^n -> R. In this example the objective function is the
% Rosenbrock objective and n = 10000
%
% To run this script type:
%
% > example_1
%
%-------------------------------------------------------------------------%
% 04/27/21, J.B, Preparation for release
% 06/02/22, J.B., Initial tests

clc;
clear;

addpath(genpath('../ALGS'));
addpath(genpath('../EXTERNAL'));
addpath(genpath('../AUXILIARY'));

% Rosenbrock objective function and gradient
% Inputs to LSR_SC
func = @(x)( rosen_obj(x) );
grad = @(x)( rosen_grad(x) );

% Problem dimension
n = 100;

% This example uses default parameters, which are described in
% the comments of LSR1_SC.m. Only exception is to print the outputs
pars.print = 1;

% Initial guess (solution is xsol=ones(n,1), f(xsol)=0)
x0 = zeros(n,1);
x0(1) = 30;

% Example
fprintf('Example 1 ##########################################\n');    
fprintf('Rosenbrock objective: f(x)                         \n');
fprintf('                      n = %i                       \n',n);
fprintf('                     f0 = %0.5g                    \n',func(x0));
fprintf('               norm(g0) = %0.5g                    \n',norm(grad(x0)));
fprintf('\n');
fprintf('Algorithm: LMSS_SC                                  \n');
fprintf('####################################################\n');
fprintf('\n');
% Call to the optimization algorithm

% Memory parameter and scaling gamma
pars.m = 6;
pars.whichInit = 4;

[xk1,gk1,fk1,out1] = LMSS_SC(x0,func,grad,pars);
