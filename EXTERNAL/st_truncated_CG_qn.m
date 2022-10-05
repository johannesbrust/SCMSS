function [s1,As1,itn,exitflag] = st_truncated_CG_qn(g,Psi,invM,gamma_k,...
    delta,sr1_mem,convergtol,itn_max,show,varargin)    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Algorithm 7.5.1 in Trust-Region Methods by Conn, Gould, and Toint
%
% implemented by J.B. Erway, June 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% User must supply:
%
%           *       g  - the gradient!
%           *       Psi - Rectangular matrix for LSR1
%           *       invM - Square (middle) matrix for LSR1
%           *       gamma_k - Scaling value
%           *       delta - trust-region radius
%           *       sr1_mem - memory value
%           *       convergtol - tolerance for convergence
%           *       show - Flag to print output
%           *       itn_max - Maximum iterations
%           *       varargin (optional: instead of Psi, Sk = varargin{1},
%   `                   Yk = varargin{2} may be supplied
% Output: 
%
%           * s1 - approximate solution to the trsp
%           * As1 - Final matrix-vector product
%           * itn - number of iterations
%           * exitflag - Indicator of success
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 02/20/21, J.Br., Modifications for use in LSR1_SC.m
% 03/02/21, J.Br., Extension for varargin (to enable component-wise product
%                   with Psi)

%global S Y;


%%   TO DO: fix input and output descritpion; line things up

n=length(g);

itn = 1; %iteration number
itn_max = min(n,itn_max); %max number of iterations

%initialize
s0 = zeros(n,1);  %initial guess
As0= zeros(n,1);  %B*s0
g0 = g;
v0 = g0;
p0 = -v0;
exitflag = 0;
%show = 0; 

%set up Psi for compact formulation to perform mults with B
% SY          = S'*Y;
% SS          = S'*S;	
% invM        = tril(SY)+ tril(SY,-1)'-gamma_k.*SS;
% invM        = (invM+(invM)')./2;  %symmetrize it
% Psi         = Y-gamma_k.*S;



while ((itn<itn_max) & (exitflag==0))

  if nargin==9  
    Ap0 = (gamma_k)*p0  + Psi(:,1:sr1_mem) * (invM(1:sr1_mem,1:sr1_mem) \ (Psi(:,1:sr1_mem)'*p0));  %test this somehow
  else
      p1_ = varargin{2}'*p0 -gamma_k.*(varargin{1}'*p0);
      p2_ = invM(1:sr1_mem,1:sr1_mem) \ p1_;
      Ap0 = (gamma_k)*p0  + varargin{2}*p2_ -varargin{1}*(gamma_k.*p2_);
  end
  kappa = p0'*Ap0;
  if (kappa<=0)
    sigmaST = -(s0'*p0)+sqrt((s0'*p0)^2+(norm(p0)^2)*(delta^2-norm(s0)^2));
    sigmaST = sigmaST/(norm(p0)^2);
    s1 = s0+sigmaST*p0;
    As1 = As0+sigmaST*Ap0;
    return;
  end

  alpha= g0'*v0/kappa;
  ip = s0+alpha*p0;
  if ( ip'*ip >= (delta^2) )
    sigmaST = -(s0'*p0)+sqrt((s0'*p0)^2+(norm(p0)^2)*(delta^2-norm(s0)^2));
    sigmaST = sigmaST/(norm(p0)^2);
    s1 = s0+sigmaST*p0;
    As1 = As0+sigmaST*Ap0;
    return;
  end

  s1 = s0+alpha*p0;
  As1 = As0+alpha*Ap0;
  g1 = g0 + alpha*Ap0;
  v1 = g1;
  betaST = (g1'*v1)/(g0'*v0);
  p1 = -v1+betaST*p0;

  %test termination criteria
  if norm(g1)<=convergtol
    break;
  end
  
  %shift indices
  itn=itn+1;  %after gone through loop 1x, itn=2
  s0 = s1;
  As0 = As1;
  g0 = g1;
  v0 = v1;
  p0 = p1;


  
end


%test if exceeded max iteration
  if itn  >=itn_max
    exitflag=1;
  end
  


if show >=1
   fprintf('\nST:  Itn     gnorm       norm(s)');
   fprintf('     gTp');
   fprintf('\n     ---   ----------  ----------  -------- ');
   fprintf('\n %2d %11.4e %11.4e %7.1e ', itn, norm(g1), norm(s1),g'*s1)
end
