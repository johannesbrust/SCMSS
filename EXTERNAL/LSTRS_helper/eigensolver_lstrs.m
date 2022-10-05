function [ NCONV,LAMBDA1,Y1,LAMBDA2,Y2,V1 ] = eigensolver_lstrs( varargin )
% eigensolver_lstrs is used as the eigensolver to be used within LSTRS and
% the limited memory matrix B. This function uses matrix-vector products
% based on the compact representation of B:
% Bx = gammaIx + Psi M \ Psi'x
%
% This function's interface is taken from "lstrs.m" for compatibility.
% By default is calls MATLAB's 'eigs' function, in order to find the
% two smallest algebraic eigenvalues.
%
%   Input Parameters:
%   H:              Hessian matrix or matrix-vector multiplication routine
%   ALPHA:          key scalar parameter in LSTRS
%   G:              gradient vector in objective function 
%   HPAR:           parameters for matrix-vector multiplication routine H
%                   See 1. above
%   EIGENSOLVERPAR: structure containing parameters for eigensolver routine
%   MESSAGE_LEVEL:  as above
% 
%   Output Parameters:
%   NCONV:          number of converged eigenvalues
%   LAMBDA1, Y1:    smallest eigenvalue and corresponding eigenvector
%   LAMBDA2, Y2:    another eigenvalue and corresponding eigenvector
%                   (in general, the second smallest eigenvalue or close to it)
%   V1:             First column of Lanczos basis matrix (or first eigenvector)
%
% Part of article SC-SR1: MATLAB Software for Solving Shape-Changing L-SR1
% Trust-Region Subproblems, (Brust,Burdakov,Erway,Marcia,2018-2021)     [1]
%
%--------------------------------------------------------------------------    
% 02/10/21, J.B., Implementation of "custom" eigensolver

% H,alpha,...
%     g,Hpar,eigensolverpar,message_level

switch nargin
    case 2
        Balpha = varargin{1};
        eigensolverpar = varargin{2};
        H = Balpha.H;
        Hpar = Balpha.Hpar;
        message_level = 0;
        g = Balpha.g;
        alpha = Balpha.alpha;
    case 6
        H = varargin{1};
        alpha = varargin{2};
        g = varargin{3};
        Hpar = varargin{4};
        eigensolverpar = varargin{5};
        message_level = varargin{6};
end

% Initializations
if isfield(eigensolverpar,'K') 
    K = eigensolverpar.K;
else
    K = 2; 
end
if -1 < message_level && message_level < 3
    OPTS.disp = message_level;
else
    OPTS.disp = 0;
end

n = length(g);

% Matrix-vector product (bordered matrix)
AFUN = @(x)(    [alpha*x(1) + g'*x(2:end);...
                H(x(2:end),Hpar) + g*x(1) ]);

% EIGS options
OPTS.issym = true;
OPTS.isreal = true;
sigma = 'SA'; % smallest algebraic eigenvalues

% Call to EIGS
[V,D,FLAG] = eigs(AFUN,(n+1),K,sigma,OPTS);

% Allocating outputs
d = diag(D);
[ds,idx] = sort(d);

% Convergence output
if FLAG == 0; NCONV = 2; else NCONV = 0; end;

LAMBDA1 = ds(1);
Y1 = V(:,idx(1));

LAMBDA2 = ds(2);
Y2 = V(:,idx(2));

V1 = V(:,idx(1));


