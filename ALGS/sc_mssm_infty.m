function [p_star,iExit,Bp_star]= sc_mssm_infty(g, S_or_Psi, Y_or_invM, gamIn,...
    delta, flag, show, varargin)
%sc_mssm_infty: Implementation of the shape-changing infinity norm
% trust-region subproblem solver for the MSS matrix
%
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
% 06/10/22, J.B., Updates for the MSS matrix
% 08/04/22, J.B., Preparation for the dense initialization
% 08/11/22, J.B., Improving the computation of B p_star with the
%                   dense initialization

tic
iExit = 0;
optL.SYM = true;    % Option for solving symmetric linear systems using
% the build-in function "linsolve"

%Check the input parameters for errors
if ((nargin<7) || (9<nargin)) || (delta<0) || ...
        ((flag~=0) && (flag~=1)) || ((show~=0) && (show~=1))
    
    fprintf('\nError in the inputs\n\n');
    iExit = 1;
    p_star = 0;
    return;
end

% If there are no stored updates, p=-g;
numMem=size(S_or_Psi,2);
if numMem == 0
    p_star = - g/norm(g);
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
        PsiPsi = Psi(:,1:numMem)'*Psi(:,1:numMem);
    else
        PsiPsi = varargin{1};
    end
elseif nargin == 9
    PsiPsi = varargin{1};
    invM = varargin{2};
end

%% Initialization of constants
tau = sqrt(eps);     %tolerance for eigenvalues
Psi_m = size(g,1);  %problem dimension

%% Memory allocation
% Names of variables correspond to the names in the manuscript.
% Many variables are allocated as function outputs, as suggested by
% "Code Analyzer" (MATLAB)
R = zeros(numMem,numMem);
w_star = zeros(numMem,1);

%% Algorithm 3 begins

%%Algorithm 1 in the paper: ComputeSpectral

[Ldl,Dldl,Pldl] = ldl(PsiPsi,'vector');  %ldl decomposition of PsiPsi'
Dldl = diag(Dldl);  %vectorize Dldl
maskldl = find(Dldl>(sqrt(eps)*max(Dldl)));  %generate mask to find lin indep columns of Psi

rankPsi = length(maskldl);     %compute rank of Psi
Pldl_fullRank = Pldl(maskldl);  %index of lin indep columns of Psi (used later)
R(1:rankPsi,:) = diag(sqrt(Dldl(maskldl)))*Ldl(:,maskldl)';

%eigenvalue and partial eigenbasis computation
invP(Pldl) = 1:numMem;

% Uncommenting next line computes the matrix RMR using the "backslash"
% operation. Currently this implementation uses the "linsolve"
%RMR = R(1:rankPsi,invP)* (invM(1:sr1_mem,1:sr1_mem)\R(1:rankPsi,invP)'); %form RP^TMPR

%RMR = R(1:rankPsi,invP)*(linsolve(invM(1:sr1_mem,1:sr1_mem),R(1:rankPsi,invP)',optL)); %form RP^TMPR
RMR = R(1:rankPsi,invP)*(invM(1:numMem,1:numMem)*R(1:rankPsi,invP)'); %form RP^TMPR
RMR = (RMR+(RMR)')./2;  %symmetrize RMR' to avoid complex
RMR(isnan(RMR)) = 0.0;
% Check on well definedness (un-commenting the next lines will include
% a safeguard against "NaN" values in RMR)
% if sum(sum(isnan(RMR))) > 0
%   p_star  = - g/norm(g);
%   return;
% end

%roundoff error
[U, D ]      = eig(RMR);

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

%form v_parallel
v_parallel      = zeros(rankPsi,1);
for i=1:rankPsi
    if abs(g_parallel(i))<delta*abs(Lambda(i)) && Lambda(i)>tau
        v_parallel(i) = -g_parallel(i)/Lambda(i);
    elseif  abs(g_parallel(i))<tau && abs(Lambda(i))<tau
        v_parallel(i) = -delta/2;
    elseif  abs(g_parallel(i))>tau && abs(Lambda(i))<tau
        v_parallel(i) = -sign(g_parallel(i))*delta;
    elseif  abs(g_parallel(i))<tau && (Lambda(i))<-tau
        v_parallel(i)=delta;
    else
        v_parallel(i)=-(delta/abs(g_parallel(i)))*g_parallel(i);
    end
end


%% Algorithm 2: ComputeW and computing p_star
% Also computes scalar beta
hasBeta = 0;
beta    = 0.0;
if (gamma_kp>0) && (a_kp2<=delta*gamma_kp)
    %w_star = -(1/gamma_k)*g;
    hasBeta   = 1;
    beta      = -(1/gamma_kp);
    w_star    = beta.*g;
elseif (gamma_kp<=0) && (a_kp2<tau)
    for i=1:numMem+1
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
            w_star = -(gamma_kp/(Pei))*ei;
            break;
        end
    end
else
    hasBeta = 1;
    beta    = -(delta/(a_kp2));
    %w_star = -(delta/(a_kp2))*g;
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


% For dense initializations compute Bp_star
% if hasDense == 1
%     
%     p_starP = U'*(R(1:rankPsi,maskldl)'\Psi(:,Pldl_fullRank)'*p_star);
%     p_starP = (Lambda-gamma_kp).*p_starP;
%     p_starP = Psi(:,Pldl_fullRank)*(R(1:rankPsi,maskldl)\...
%                 U*p_starP);
%     Bp_star = gamma_kp.*p_star + p_starP;
%     
%     Bp_star_ = gamma_kp.*w_star + V2P(:,2);
%     er2 = norm(Bp_star - Bp_star_);
%     
% else
%     Bp_star = 0;
% end

%% End of Algorithm 3

time_counter=toc;


%%% output for (P,infty)-norm

if show==1
    fprintf('\nVerbose mode on:');
    if flag==1
        fprintf('\n1. inputs were Psi and M^{-1}');
    else
        fprintf('\n1. inputs were S and Y; converted to compact form.');
    end
    fprintf('\n2. gamma =%5.3e', gamma_k);
    %print to a file what I need for paper
    fres = fopen('res1000out.text', 'a+');
    fprintf(fres, ' {%3.2e}', Psi_m);
    fprintf(fres, '& {%3.2e} ', gamma_k);
    fprintf(fres, '& {%3.2e} ', time_counter);
    fprintf(fres, '...done\n');
    fclose(fres);
    fprintf('\nTime taken: %8.3e\n', time_counter);
end


