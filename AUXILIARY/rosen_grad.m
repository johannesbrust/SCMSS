%-------------------------- rosen_grad.m ---------------------------------%
%
% Evaluation of the rosenbrock gradients. Here:
%
% f(x)  = sum (even(x) - odd(x)^2)^2 + (1 - odd(x))^2,
%
% x*    = 1 (solution)
%
% Input: 
% x 
%
% Output: 
% g (gradient)
%
% This function is modified from a function by Burdakov et al. (LMTR,'16)
%-------------------------------------------------------------------------%
% 08/09/19, J.B. 
% 04/22/20, J.B., Modification of output values (only gradient)

function [g] = rosen_grad(x,varargin)
% Rosenbrock function
[n,m] = size(x);

if n<m
    x=x';
    n=m;
end

if mod(n,2)
    error('problem dimension is not an even number!');
end

eidx = 2:2:n;
oidx = 1:2:(n-1);

f1   = x(eidx)-x(oidx).^2;
f2   = 1-x(oidx);
%f    = sum(f1.^2+f2.^2);


g = zeros(n,1);

g(eidx) = 2*f1;
g(oidx) = -4*x(oidx).*f1 - 2*f2;

% if nargout>2 % compute Hessian also
%     %H = zeros(n,n);
%     
%     Hd = zeros(n,1);
%     Hd(eidx) = 2;
%     Hd(oidx) = -4*f1 + 8*x(oidx).^2;
%     
%     Ho = zeros(n-1,1);
%     Ho(oidx) = -4*x(oidx);
%     
%     H = diag(Hd,0) + diag(Ho,1) + diag(Ho,-1);
%     
% %     H(eidx,eidx) = 2;
% %     H(eidx,oidx) = repmat(-4*x(oidx)',n/2,1);
% %     H(oidx,eidx) = repmat(-4*x(oidx)',n/2,1);
% %     H(oidx,oidx) = repmat(-4*f1' + 8*x(oidx)'.^2 + 2,n/2,1);
%     
% end
    
