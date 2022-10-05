%------------------------ test_lbfgsb_cutest -----------------------------%
%
% Initial test using the L-BFGS-B method on a problem loaded from CUTEST
%
%-------------------------------------------------------------------------%
% 01/20/22, J.B., initial implementation

clc;
clear;

%addpath(genpath('../ALGS'));
addpath(genpath('../EXTERNAL'));
addpath(genpath('../AUXILIARY'));

% Setup CUTEst problem
probname = 'BOX';
CUTEst_init
cmdcutest = ['cutest2matlab_osx $MASTSIF/' probname];
unix(cmdcutest);
prob = cutest_setup();
x0 = prob.x;
n = size(x0,1);
clc;

% Function and gradient
func = @cutest_fun;
grad = @cutest_gradient;

% Example
fprintf('Example LBB CT #####################################\n');    
fprintf('Function objective:   f(x)                         \n');
fprintf('                      n = %i                       \n',n);
fprintf('                     f0 = %0.5g                    \n',func(x0));
fprintf('               norm(g0) = %0.5g                    \n',norm(grad(x0)));
fprintf('\n');
fprintf('Algorithm: L-BFGS-B                                \n');
fprintf('####################################################\n');
fprintf('\n');

% Call to the optimization L-BFGS-B algorithm
% Setting options
l = -inf(n,1);
u = inf(n,1);
funcLB = @(x)fminunc_wrapper( x, func, grad);
pars = struct(  'factr', 0,... % Suppress this stopping condition
                'pgtol', 5e-4,...
                'm', 5,...
                'x0',x0,...
                'maxIts',5000,...
                'maxTotalIts',25000);

% Call to solver            
[xk1,fk1,out1] = lbfgsb(funcLB, l, u, pars);

% Delete library files from CUTEst
delete( '*.d',...
        '*.o',...
        '*.dylib',...
        '*.f',...
        'mcutest.*');
