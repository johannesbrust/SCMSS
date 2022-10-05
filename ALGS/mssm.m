function [p_star]= mssm(prob,mem);
%
% Copyright (2021-2022): Jennifer Erway
%
%
% The technical report and software are available at 
% https://urldefense.com/v3/__http://www.wfu.edu/*erwayjb/publications.html__;fg!!Mih3wA!CsTSnLbzeha5hZsJUaTfRFdhYgtDNR5H3OYvBVXerQn_rmRMONA_fXEcshN4yyeiNPAB_q0eX-UUWrC1NQ$ 
% https://urldefense.com/v3/__http://www.wfu.edu/*erwayjb/software.html__;fg!!Mih3wA!CsTSnLbzeha5hZsJUaTfRFdhYgtDNR5H3OYvBVXerQn_rmRMONA_fXEcshN4yyeiNPAB_q0eX-WqKPz2SA$ 
%
%
% This code is distributed under the terms of the GNU General Public
% License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire 
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose.
%---------------------------------------------------------------------------    
% 
% Inputs: 
%         g: the gradient vector and right-hand side of the first
%                optimality condition
%
% Output: p_star  solution
%         iExit:  0=no termination condition met, 1=normG<grad_tol, 
%                 2=itn>=itn_max, 3= delta too small
%                 7 = error in inputs

myfunc = @cutest_obj;
mygrad = @cutest_grad;
x0 = prob.x

iExit       = 0; 

if (nargin~=4) 
    fprintf('\nError in the inputs\n\n');
    iExit = 7;
    p_star = 0;
    return;
end

%constants
n=size(x0,1);
mem_max = mem; 

%% memory allocation
Psi             = zeros(n,2*mem_max);
S               = zeros(n,mem_max);
Y               = zeros(n,mem_max);
SS              = zeros(mem_max,mem_max);
SY              = zeros(mem_max,mem_max);
BigG            = zeros(mem_max,mem_max);
YY              = zeros(mem_max,mem_max);


%tolerance used for rank and  y's
tol     = 1e-4; 


%trust-region constants
eta_1          =  0.01; 
eta_2          =  0.75; 
gamma_1        =  0.5;  
gamma_2        =  2.0;  
delta_min = 100*eps;  

%initializations
x = x0;
gamma_k = 1; 
gamma_perp = 1; 
delta = 1;
k = 0; 
itn = 0; 
f = myfunc(x);
g = mygrad(x);

%constants
itmax = 2*n;
func_max = 100*n;


%convergence
gtola   = 1e-05;
normG = norm(g);
grad_tol  = max(gtola*norm(g),gtola); 

%% initial test
if  (normG<gtola)
  fprintf('gradient is already too close to zero to proceed');
  p_star = x;
  iExit=1;
  return;
end

normx = norm(x);
update_pairs = 0;


while update_pairs==0 & iExit==0
  p  = -(1/gamma_k)*g; 
  Bp = -g;

  pnorm = norm(p); 
  xa = x + p; 
  fa = myfunc(xa); 
  ga = mygrad(xa); 
  
  ratio_test = (f-fa)/(-g'*p-.5*(p'*Bp));
  
  update_pairs = 0;
  if (ratio_test >= eta_1) 
    
    s = p;
    y = ga-g;
    update_pairs = 1;
    gtp=g'*p; 
    x= xa; f = fa; g = ga;
    normx = norm(x); normG = norm(g);
    
    if (ratio_test >= eta_2)
      if norm(p)<=0.8*delta
      else
	delta = gamma_2*delta;
      end
    else
    end 
  else     %reject iterate, reduce delta
    s = p;
    y = ga-g;
    update_pairs = 1;
    delta = gamma_1*delta;
  end
  

    if update_pairs==1

      S(:,1)=s;
      Y(:,1)=y;
      SS(1,1) = s'*s;
      SY(1,1) = s'*y;
      YY(1,1) = y'*y;
      
      
      %gamma and gamma perp cannot be negative
      gamma_k = (y'*y)/(y'*s);   
      gamma_perp = gamma_k;  
      if (gamma_k<tol) |  (gamma_k>(1/tol)) | (~isfinite(gamma_k))
	gamma_k = 1;
      end
      if (gamma_perp<tol) |  (gamma_perp>(1/tol)) | (~isfinite(gamma_perp))
	gamma_perp = 1;
      end
      
      Psi(:,1:2) = [s (y-gamma_k.*s)];
      iSS = 1/(SS(1,1));
      BigG(1,1) = s'*y;
      M = [-gamma_k*(iSS)-(iSS(1,1)*BigG(1,1)*iSS(1,1))  iSS; iSS 0];   
      rankS = 1;
      k = k+1;
      
    end 
 
    itn = itn+1;  
    
    
    if (normG<grad_tol)
      iExit = 1;
    elseif  (itn>=itmax)
      iExit = 2;
    elseif (delta<delta_min)
      iExit=3;
    else
    end
    
end %update_pairs==0 & iExit==0

first_flag = 1;

%main loop

while iExit==0 & (itn<itmax) 

if first_flag==0 & update_pairs==1


    if k<mem_max 
      k = k+1;
      S(:,1:k) = [s S(:,1:k-1)]; 
      Y(:,1:k) = [y Y(:,1:k-1)];
    else
      S = [s S(:,1:k-1)];
      Y = [y Y(:,1:k-1)];
    end

    %Ensure that S is full rank so that S'S in Psi is invertible
    [qS,rS,eS] = ldl(S(:,1:k)'*S(:,1:k),'vector');
    maskS = find(abs(diag(rS))>(tol*max(abs(diag(rS)))));    
    rankS = size(maskS,1);
    if rankS<k
      eS = eS(:,maskS);
      Stilde = S(:,1:k);
      Stilde = Stilde(:,eS);
      Ytilde = Y(:,1:k);
      Ytilde = Ytilde(:,eS);
    else
      Stilde = S(:,1:k);
      Ytilde = Y(:,1:k);
    end

    SS(1:rankS,1:rankS) = Stilde'*Stilde;
    SY(1:rankS,1:rankS) = Stilde'*Ytilde;
    YY(1:rankS,1:rankS) = Ytilde'*Ytilde;
    
    % form middle matrix
    BigG(1:rankS,1:rankS) = triu(SY(1:rankS,1:rankS))+triu(SY(1:rankS,1:rankS),1)';
    iSS= inv(SS(1:rankS,1:rankS));

    %initial B0
    gamma_k = (y'*y)/(y'*s);
    gamma_perp = gamma_k;
    if (gamma_k<tol) |  (gamma_k>(1/tol)) | (~isfinite(gamma_k))
      gamma_k = 1;
    end
    if (gamma_perp<tol) |  (gamma_perp>(1/tol)) |  (~isfinite(gamma_perp))
      gamma_perp = 1;
    end
    

    %compact formulation
     M = [-gamma_k*iSS-iSS*BigG(1:rankS,1:rankS))*iSS	iSS; iSS zeros(rankS,rankS)];
     Psi(:,1:2*rankS) = [Stilde(:,1:rankS)  Ytilde(:,1:rankS)];

end %end of update_pairs==1

if update_pairs ==1 | first_flag==1

   %fix psi'psi without forming psi
   PsiPsi = Psi(:,1:2*rankS)'*Psi(:,1:2*rankS);
   [Ldl,Dldl,Pldl] = ldl(PsiPsi,'vector');
   Dldl=diag(Dldl);  
   maskldl = find(Dldl>tol*max(Dldl)); 
   rankPsi = length(maskldl);     
   Pldl_fullRank = Pldl(maskldl);  
   R(1:rankPsi,1:2*rankS) = diag(sqrt(Dldl(maskldl)))*Ldl(:,maskldl)';
   
   % Eigenvalue and partial eigenbasis computation
   invP = [ ];
   invP(Pldl)=1:2*rankS; 
   RMR = R(1:rankPsi,invP)*M*R(1:rankPsi,invP)'; 
   
   %form RP^TMPR
   RMR         = (RMR+(RMR)')./2;
   [U,D]      = eig(RMR);
   
   %find Lambda
   [D,ind]     = sort(real(diag(D)));  
   U           = U(:,ind);            
   Lambda      = D + gamma_k.*ones(rankPsi,1);
   

    % for solver
    P_parallel = Psi(:,Pldl_fullRank)*inv(R(1:rankPsi,maskldl))*U;
    
end %first_flag==1 & update_pairs==1


%Now solve the TRSP
%Johannes: available variables....
%%%  1. Stilde (full rank S) and Ytilde (corresponding Y)
%%%  2. SY,YY,SS (products with S and Y)
%%%  3. P_parallel
%%%  4. iSS  = (SS)^{-1}
%%%  5. M = middle matrix in compact formulation
%%%  6. Psi
[p,Bp]= CALL_TRSP_SOLVER_here


%%% possibly update x,f,g,delta
pnorm = norm(p); 
xa = x + p; 
fa = myfunc(xa); 
ga = mygrad(xa); 

ratio_test = (f-fa)/(-g'*p-.5*(p'*Bp));

update_pairs = 0;
if (ratio_test >= eta_1) 
  s = p;
  y = ga-g;
  update_pairs = 1;
  gtp=g'*p;
  x= xa; f = fa; g = ga;
  normx = norm(x); normG = norm(g);
  
  if (ratio_test >= eta_2)
    if norm(p)<=0.8*delta
    else
      delta = gamma_2*delta;
    end
  else
  end
else    
  s = p;
  y = ga-g;
  update_pairs = 1;
  delta = gamma_1*delta;
end

itn = itn+1;         
first_flag=0;


if (normG<grad_tol)
 iExit = 1;
 time_counter=toc;
elseif  (itn>=itmax)
iExit = 2;
elseif (delta<delta_min)
  iExit=3;
else 
end



end % main while loop    


p_star = x;

