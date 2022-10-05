function [p_star,iExit,Bp_star]= sc_mssm_2(g,S_or_Psi,Y_or_invM, gamIn,...
    delta,flag,show,varargin)
%sc_mssm_2: Implementation of a the shape-changing l2 norm trust-region
% subproblem solver
%%
% INPUTS:
%         g: the gradient vector and right-hand side of the first
%                optimality condition
%         (S, Y) or (Psi, invM): used to define the quasi-Newton matrix
%         gamIn: scaling of (dense) initial Hessian matrix approximation B0,
%                i.e., B0 = P diag(gamma_k*I,gamma_p*I) P'
%         delta: the trust-region radius
%         flag:  flag==0 means given (S,Y); flag=1 means (Psi,invM)
%         show:   0=Runs silently; 1=Verbose setting
%         varargin: Optional argument of precomputed
%                   Psi'*Psi (nargin=8) and/or
%                   invM (nargin=9)
%
% OUTPUTS: p_star: the optimal solution to the trust-region subproblem
%         iExit:  0=no errors; 1=error
%         Bp_star: B*p_star when a dense initialization is used
%
%--------------------------------------------------------------------------
%
% Copyright (2022): Johannes J. Brust, Jennifer Erway, Roummel Marcia
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
%
%--------------------------------------------------------------------------
% 04/27/21, J.B., Preparation for release
% 06/09/22, J.B., Update for the MSS matrix
% 08/11/22, J.B., Improving the computation of B p_star with the
%                   dense initialization

tic
iExit       = 0;     % runtime error flag
optL.SYM = true;

%Check the input parameters for errors
if ((nargin<7) || (9<nargin)) || (delta<0) || ...
        ((flag~=0) && (flag~=1)) || ((show~=0) && (show~=1))
    
    fprintf('\nError in the inputs\n\n');
    iExit = 1;
    p_star = 0;
    return;
end

%if there are no stored updates, p=-g;
sr1_mem=size(S_or_Psi,2);
if sr1_mem == 0
    p_star  = - g/norm(g);
    return;
end

% Dense initialization
if length(gamIn) == 1
    gamma_k  = gamIn(1);
    gamma_kp = gamma_k;
    hasDense = 0;
elseif length(gamIn) == 2
    gamma_k  = gamIn(1);
    gamma_kp = gamIn(2);
    hasDense = 1;
end

if nargin == 7||nargin==8
    if flag == 1  % handed Psi and M
        Psi = S_or_Psi;
        invM = Y_or_invM;
    else  %handed S and Y
        %set up Psi and M^{-1}
        %S=S_or_Psi; Y=Y_or_invM
        SY          = S_or_Psi'*Y_or_invM;
        SS          = S_or_Psi'*S_or_Psi;
        invM        = tril(SY)+ tril(SY,-1)'-gamma_k.*SS;
        invM        = (invM+(invM)')./2;  %symmetrize it
        Psi         = Y_or_invM-gamma_k.*S_or_Psi;
    end
    if nargin==7
        PsiPsi = Psi(:,1:sr1_mem)'*Psi(:,1:sr1_mem);
    else
        PsiPsi = varargin{1};
    end
elseif nargin == 9
    PsiPsi = varargin{1};
    invM = varargin{2};
end


%% Initialization of constants
maxIter = 100;   % maximum number of iterations for Newton's method
N_iter = 0;     % number of Newton iterations
tau = sqrt(eps);     %tolerance for eigenvalues
Psi_m = size(g,1);  %problem dimension


%% Memory allocation
% Names of variables correspond to the names in the manuscript.
% Many variables are allocated as function outputs, as suggested by
% "Code Analyzer" (MATLAB)
R = zeros(sr1_mem,sr1_mem);
w_star = zeros(Psi_m,1);
e1 = zeros(sr1_mem,1);  e1(1)=1; % First basis vector

% Initializing optimal Lagrange multiplier.
sigma_parallel = 0;  %#ok<NASGU>

%% Algorithm 4 begins

%%Algorithm 1 in the paper: ComputeSpectral

[Ldl,Dldl,Pldl] = ldl(PsiPsi,'vector');  %ldl decomposition of PsiPsi'
Dldl=diag(Dldl);  %vectorize Dldl
maskldl = find(Dldl>(sqrt(eps)*max(Dldl)));  %generate mask to find lin indep columns of Psi

rankPsi = length(maskldl);     %compute rank of Psi
Pldl_fullRank = Pldl(maskldl);  %index of lin indep columns of Psi (used later)
R(1:rankPsi,:) = diag(sqrt(Dldl(maskldl)))*Ldl(:,maskldl)';

% Eigenvalue and partial eigenbasis computation
invP(Pldl)=1:sr1_mem;
%RMR = R(1:rankPsi,invP)*(linsolve(invM(1:sr1_mem,1:sr1_mem),R(1:rankPsi,invP)',optL)); %form RP^TMPR
RMR = R(1:rankPsi,invP)*(invM(1:sr1_mem,1:sr1_mem)*R(1:rankPsi,invP)'); %form RP^TMPR

% Uncommenting next line computes the matrix RMR using the "backslash"
% operation. Currently this implementation uses the "linsolve"

%RMR = R(1:rankPsi,invP)* (invM(1:sr1_mem,1:sr1_mem)\R(1:rankPsi,invP)');

% Form RP^TMPR
RMR = (RMR+(RMR)')./2;  %symmetrize RMR' to avoid complex
RMR(isnan(RMR)) = 0.0;
%roundoff error

% Check on well definedness (un-commenting the next lines will include
% a safeguard against "NaN" values in RMR)

% if sum(sum(isnan(RMR))) > 0
%   p_star  = - g/norm(g);
%   return;
% end

[U,D]      = eig(RMR);

%% Eliminate complex roundoff error then sort eigenvalues and eigenvectors
[D,ind]     = sort(real(diag(D)));  %D is now a vector; low to high
U           = U(:,ind);


% Compute Lambda as in Equation (7)
Lambda      = D + gamma_k.*ones(rankPsi,1);
Lambda      = Lambda.*(abs(Lambda)>tau); %thresholds

%% end of algorithm 1


%% Compute g_parallel using (14) or (13)
% Product computed using Yk,Sk if appropriate number of inputs present
% Otherwise Psi is directly used
if nargin == 9 % (14)
    g_parallel  = U'*(R(1:rankPsi,maskldl)'\...
        ((Y_or_invM(:,Pldl_fullRank)'*g)-gamma_k.*(S_or_Psi(:,Pldl_fullRank)'*g)));
else % (13)
    g_parallel  = U'*(R(1:rankPsi,maskldl)'\Psi(:,Pldl_fullRank)'*g);
end



%% Compute ||g_perp||
a_kp2       = sqrt(abs(g'*g-g_parallel'*g_parallel));
if a_kp2 < tau  %fix small values of a_kp2
    a_kp2 = 0;
end

if (Lambda(1)>tau)  %line 17 of Algorithm 4
    
    %case 1
    if (norm(g_parallel./Lambda)<=delta)  %try unconstrained minimizer
        sigma_parallel = 0;
    else
        [sigma_parallel,N_iter] = Newton(0,maxIter,tau,delta,Lambda,g_parallel);
    end
    v_parallel = -g_parallel./(Lambda+sigma_parallel);
    
elseif (abs(Lambda(1))<tau)    %line 24 of Algorithm 4
    
    %case 2
    index_distinct = find(abs(Lambda)>tau,1);
    
    case2_test1 = sum(abs(g_parallel(1:index_distinct-1))<tau) == (index_distinct-1);
    case2_test2 = (norm(g_parallel(index_distinct:end)./Lambda(index_distinct:end))<=delta);
    if isempty(case2_test1); case2_test1 = 0; end;
    if isempty(case2_test2); case2_test2 = 0; end;
    if case2_test1 && case2_test2
        sigma_parallel = 0;
        v_parallel = [zeros(index_distinct-1,1);-g_parallel(index_distinct:end)./Lambda(index_distinct:end)];
    else
        sigma_hat = max((abs(g_parallel)/delta) - Lambda);
        [sigma_parallel,N_iter]= Newton(max(sigma_hat,0),maxIter,tau,delta,Lambda,g_parallel);
        v_parallel = -g_parallel./(Lambda+sigma_parallel);
    end
    
else  %Lambda(1)<tau (line 25 of Algorithm 4)
    
    %case 3
    index_distinct = find(abs(Lambda-Lambda(1))>tau,1);
    if sum(abs(g_parallel(1:index_distinct-1))<tau) == (index_distinct-1)
        sigma_parallel=-Lambda(1);
        v = [zeros(index_distinct-1,1);-g_parallel(index_distinct:end)./(Lambda(index_distinct: ...
            end)-Lambda(1))];
        if (norm(v)<=delta)
            %fix for the "hard case"
            alpha = sqrt(delta^2-norm(v)^2);
            v_parallel = v + alpha*e1(1:length(v));
            %v_parallel = v + alpha*e1;
        else  %norm(v)>delta
            sigma_hat = max((abs(g_parallel)/delta) - Lambda);
            [sigma_parallel,N_iter]=  Newton(max(sigma_hat,0),maxIter,tau,delta,Lambda,g_parallel);
            v_parallel = -g_parallel./(Lambda+sigma_parallel);
        end
        
    else  %(g_parallel)_i ~= 0 for some 1<=i<=r
        sigma_hat = max((abs(g_parallel)/delta) - Lambda);
        [sigma_parallel,N_iter]= Newton(max(sigma_hat,0),maxIter,tau,delta,Lambda,g_parallel);
        v_parallel = -g_parallel./(Lambda+sigma_parallel);
    end
end %end of test on (Lambda(1)>tau)

if N_iter>=maxIter
    iExit = 1;
end



%% Algorithm 2: ComputeW and computing p_star from  (Section 3.5)
hasBeta = 0;
beta = 0.0;
if (gamma_k>0) && (a_kp2<=delta*gamma_k)
    hasBeta   = 1;
    beta      = -(1/gamma_k);
    w_star    = beta.*g;
elseif (gamma_k<=0) && (a_kp2<tau)
    for i=1:sr1_mem+1
        ei= (double(1:Psi_m==i))';  %generates ith standard basis vector
        % Branching according to number of inputs
        if nargin == 9
            P_par_ei = U'*(R(1:rankPsi,maskldl)'\...
                ((Y_or_invM(:,Pldl_fullRank)'*ei)-gamma_k.*(S_or_Psi(:,Pldl_fullRank)'*ei)));
        else
            P_par_ei = U'*(R(1:rankPsi,maskldl)'\(Psi(:,Pldl_fullRank)'*ei));
        end
        P_par_ei_norm = norm(P_par_ei);
        if (P_par_ei_norm<(1-tau) )
            Pei = sqrt(1-P_par_ei_norm^2);
            w_star = -(gamma_k/(Pei))*ei;
            break;
        end
    end
else
    hasBeta = 1;
    beta    = -(delta/(a_kp2));
    w_star  = beta.*g;
end  %end of computing w_star (equation (38))

% Updated p_star computation based on whether (14) or (13) is used
% The intermediate variable "p_star_w1" is used to simplify the
% computations. Mathematically, it is "p_star_w1 = P_parallel'*w_star".
if hasBeta == 1
    p_star_w1 = beta.*g_parallel;
else
    if nargin==9 % (14)
        p_star_w1 = R(1:rankPsi,maskldl)'\...
            ((Y_or_invM(:,Pldl_fullRank)'*w_star)-gamma_k.*(S_or_Psi(:,Pldl_fullRank)'*w_star));
    else % (13)
        p_star_w1 = R(1:rankPsi,maskldl)'\((Psi(:,Pldl_fullRank)')*w_star);
    end
end
if nargin==9 % (14)
    p1(1:rankPsi,1) = R(1:rankPsi,maskldl)\(U*(v_parallel-p_star_w1));
    p_star_parallel = Y_or_invM(:,Pldl_fullRank)*p1(1:rankPsi,1) - ...
        S_or_Psi(:,Pldl_fullRank)*(gamma_k.*p1(1:rankPsi,1));
else % (13)
%     p_star_parallel = Psi(:,Pldl_fullRank)*(R(1:rankPsi,maskldl)\...
%         (U*(v_parallel-p_star_w1)));

    if hasDense == 1
        
        %vmp = v_parallel-p_star_w1;
        V2  = [v_parallel-p_star_w1, Lambda.*v_parallel - gamma_kp.*p_star_w1];
        V2P = Psi(:,Pldl_fullRank)*(R(1:rankPsi,maskldl)\...
            (U*(V2)));
        
        p_star_parallel = V2P(:,1);
        
        Bp_star = gamma_kp.*w_star + V2P(:,2);
    
    else
        
        p_star_parallel = Psi(:,Pldl_fullRank)*(R(1:rankPsi,maskldl)\...
            (U*(v_parallel-p_star_w1)));
        
        Bp_star = 0;
        
        %er1 = norm(p_star_parallel - p_star_parallel_);        
    
    end
end

p_star = p_star_parallel + w_star;

%% End of Algorithm 4

time_counter=toc;


%%% optimality check for (P,2)-norm
% Note that this is typically only exectured for obtaining information
if show==1
    fprintf('\nVerbose mode on:');
    if flag==1
        fprintf('\n1. inputs were Psi and M^{-1}');
    else
        fprintf('\n1. inputs were S and Y; converted to compact form.');
    end
    
    Bp_star2 = (gamma_k)*p_star+Psi*(inv(invM)*(Psi'*p_star)); %#ok<MINV> %(equation 5)
    
    Bp_star = gamma_k*p_star+Psi(:,Pldl_fullRank)*(inv(R(1:rankPsi,maskldl))*(RMR* ...
        (inv(R(1:rankPsi,maskldl))'*(Psi(:,Pldl_fullRank)'*p_star)))); %#ok<MINV>
    fprintf('\nDifference between norm(Bp_star2-Bp_star)=%8.3e', norm(Bp_star2-Bp_star));
    sigma_perp  = 0;
    if (gamma_k<0) || (a_kp2> abs(gamma_k)*delta)
        sigma_perp = (a_kp2/delta) - gamma_k;
    end
    PpStar      = U'*(R(1:rankPsi,maskldl)'\(Psi(:,Pldl_fullRank)'*p_star));
    CpStar      = sigma_perp*p_star+(sigma_parallel-sigma_perp).*...
        (Psi(:,Pldl_fullRank)*(R(1:rankPsi,maskldl)\(U*(PpStar))));
    
    opt1        = norm(Bp_star+CpStar + g);
    norm_PpStar   = norm(PpStar);
    opt2       =  abs(sigma_parallel*(norm_PpStar-delta));
    norm_PperpStar = sqrt(norm(p_star)^2-norm(PpStar)^2);
    opt3       = abs(sigma_perp*(norm_PperpStar-delta));
    spd_check       = zeros(3,1);
    spd_check(1)    = Lambda(1) + sigma_parallel;
    spd_check(2)    = gamma_k + sigma_perp;
    spd_check(3)    = min(spd_check(1),spd_check(2));
    fprintf('\n2. Optimality condition #1: %8.3e', opt1);
    fprintf('\n3. Optimality condition #2: %8.3e', opt2);
    fprintf('\n4. Optimality condition #3: %8.3e', opt3);
    fprintf('\n5. spd check: min eig(B+C_par) = %8.2e', spd_check(3));
    fprintf('\n6. sigma_parallel= %5.3e', sigma_parallel);
    fprintf('\n7. sigma_perp =%5.3e', sigma_perp);
    fprintf('\n8. gamma =%5.3e', gamma_k);
    %print to a file
    fres = fopen('res1000out.text', 'a+');
    fprintf(fres, ' {%3.2e}', Psi_m);
    fprintf(fres, '& {%3.2e} & {%3.2e} &  {%3.2e} ', opt1,opt2,opt3  );
    fprintf(fres, '& {%3.2e} & {%3.2e} & {%3.2e}', spd_check(3), ...
        sigma_parallel, sigma_perp);
    fprintf(fres, '& {%3.2e} ', gamma_k);
    fprintf(fres, '& {%3.2e} ', time_counter);
    fprintf(fres, '...done\n');
    fclose(fres);
    fprintf('\nTime taken: %8.3e\n', time_counter);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ phiBar, phiBar_g ] = phiBar_fg( sigma, delta, D, g_parallel )
%phiBar_f evaluates the continous extension of
%phi(sigma) = 1/ ||v|| - 1/delta and its derivative.
%
%
% Copyright (2015): Johannes Brust, Jennifer Erway, and Roummel Marcia
%
% The technical report and software are available at
% www.wfu.edu/~erwayjb/publications.html
% www.wfu.edu/~erwayjb/software.html
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
%--------------------------------------------------------------------------

%p=zeros(size(D,1),size(D,2));
m        = size(g_parallel,1);
D        = D + sigma.*ones(m,1);   %vector
eps_tol  = 1e-10;
phiBar_g = 0;

%% case 1: numerator or denominator has a zero
if ( sum( abs(g_parallel) < eps_tol ) > 0 ) || (sum( abs((D)) < eps_tol ) > 0 )
    pnorm2 = 0;
    for i = 1:m
        if (abs(g_parallel(i)) > eps_tol) && (abs(D(i)) < eps_tol)
            phiBar   = -1/delta;
            phiBar_g = 1/eps_tol;
            return;
        elseif abs(g_parallel(i)) > eps_tol && abs(D(i)) > eps_tol
            pnorm2   = pnorm2   +  (g_parallel(i)/D(i))^2;
            phiBar_g = phiBar_g + ((g_parallel(i))^2)/((D(i))^3);
        end
    end
    normP    = sqrt(pnorm2);
    phiBar   = 1/normP - 1/delta;
    phiBar_g = phiBar_g/(normP^3);
    return;
end

%% case 2: numerators and denominators are all nonzero
p      = g_parallel./D;
normP  = norm(p);
phiBar = 1/normP - 1/delta;

phiBar_g = sum((g_parallel.^2)./(D.^3));
phiBar_g = phiBar_g/(normP^3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ x,k ] = Newton(x0,maxIter,tol,delta,Lambda,g_parallel)
%Newton finds a zero of phiBar.
%
%For further details, please see the following article:
% 
% On Solving L-SR1 Trust-Region Subproblems
%
% Copyright (2015): Johannes Brust, Jennifer Erway, and Roummel Marcia
%
% The technical report and software are available at
% www.wfu.edu/~erwayjb/publications.html
% www.wfu.edu/~erwayjb/software.html
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


x = x0;  %initialization
k = 0;   %counter

[f, g]  = phiBar_fg(x,delta,Lambda,g_parallel);
f0     = eps*norm(f);  %relative error tolerance

while (k < maxIter)
    x       = x - f / g;
    [f, g]  = phiBar_fg(x,delta,Lambda,g_parallel);
    k = k + 1;
    if (norm(f) <= (norm(f0)+tol))
        return
    end
end


