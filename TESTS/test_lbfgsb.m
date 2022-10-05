%------------------------ test_lbfgsb ------------------------------------%
%
% Initial test using the L-BFGS-B method
% (on Rosenbrock objective)
%
%-------------------------------------------------------------------------%
% 01/20/22, J.B., initial implementation

clc;
clear;

%addpath(genpath('../ALGS'));
addpath(genpath('../EXTERNAL'));
addpath(genpath('../AUXILIARY'));

% Rosenbrock objective function and gradient
% Inputs to LSR_SC
func = @(x)( rosen_obj(x) );
grad = @(x)( rosen_grad(x) );

% Problem dimension
n = 10000;

% Initial guess (solution is xsol=ones(n,1), f(xsol)=0)
x0 = zeros(n,1);
x0(1) = 30;

% Example
fprintf('Example LBB ########################################\n');    
fprintf('Rosenbrock objective: f(x)                         \n');
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
                'maxIts',100000,...
                'maxTotalIts',1000000);

% Call to solver            
[xk1,fk1,out1] = lbfgsb(funcLB, l, u, pars);

