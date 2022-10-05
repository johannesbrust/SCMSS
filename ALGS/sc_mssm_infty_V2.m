function [p_star,Bp,iExit]= sc_mssm_infty_V2(g, S_or_Psi, Y_or_invM, gamma_k,...
    delta, flag, show, varargin)
%sc_sr1_infty: Implementation of the shape-changing infinity norm
% trust-region subproblem solver from article (Alg. 3)
%
% SC-SR1: MATLAB Software for Limited-Memory SR1 Trust-Region Methods,
% (Brust,Burdakov,Erway,Marcia,2018-2021)                               [1]
%
% INPUTS:
%         g: the gradient vector and right-hand side of the first
%                optimality condition
%         (S, Y) or (Psi, invM): used to define the quasi-Newton matrix
%         gamma_k: scaling of initial Hessian matrix approximation B0,
%                i.e., B0 = gamma_k*I
%         delta: the trust-region radius
%         flag:  flag==0 means given (S,Y); flag=1 means (Psi,invM)
%         show:   0=Runs silently; 1=Verbose setting
%         varargin: Optional argument of precomputed
%                   Psi'*Psi (nargin=8) and/or
%                   invM (nargin=9)
%
% OUTPUTS: p_star: the optimal solution to the trust-region subproblem
%         iExit:  0=no errors; 1=error
%
%--------------------------------------------------------------------------
%
% Copyright (2021): Johannes J. Brust, Jennifer Erway, Roummel Marcia
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
sr1_mem=size(S_or_Psi,2);
if sr1_mem == 0
    p_star = - g/norm(g);
    return;
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
tau = sqrt(eps);     %tolerance for eigenvalues
Psi_m = size(g,1);  %problem dimension

%% Memory allocation
% Names of variables correspond to the names in the manuscript.
% Many variables are allocated as function outputs, as suggested by
% "Code Analyzer" (MATLAB)
R = zeros(sr1_mem,sr1_mem);
w_star = zeros(sr1_mem,1);

%% Algorithm 3 begins

%%Algorithm 1 in the paper: ComputeSpectral

[Ldl,Dldl,Pldl] = ldl(PsiPsi,'vector');  %ldl decomposition of PsiPsi'
Dldl = diag(Dldl);  %vectorize Dldl
maskldl = find(Dldl>(sqrt(eps)*max(Dldl)));  %generate mask to find lin indep columns of Psi

rankPsi = length(maskldl);     %compute rank of Psi
Pldl_fullRank = Pldl(maskldl);  %index of lin indep columns of Psi (used later)
R(1:rankPsi,:) = diag(sqrt(Dldl(maskldl)))*Ldl(:,maskldl)';

%eigenvalue and partial eigenbasis computation
invP(Pldl) = 1:sr1_mem;

% Uncommenting next line computes the matrix RMR using the "backslash"
% operation. Currently this implementation uses the "linsolve"
%RMR = R(1:rankPsi,invP)* (invM(1:sr1_mem,1:sr1_mem)\R(1:rankPsi,invP)'); %form RP^TMPR

%RMR = R(1:rankPsi,invP)*(linsolve(invM(1:sr1_mem,1:sr1_mem),R(1:rankPsi,invP)',optL)); %form RP^TMPR
RMR = R(1:rankPsi,invP)*(invM(1:sr1_mem,1:sr1_mem)*R(1:rankPsi,invP)'); %form RP^TMPR
RMR = (RMR+(RMR)')./2;  %symmetrize RMR' to avoid complex

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


%% Algorithm 2: ComputeW and computing p_star from  (Section 3.5)
% Also computes scalar beta
hasBeta = 0;
beta    = 0.0;
if (gamma_k>0) && (a_kp2<=delta*gamma_k)
    %w_star = -(1/gamma_k)*g;
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
    p_star_parallel = Psi(:,Pldl_fullRank)*(R(1:rankPsi,maskldl)\...
        (U*(v_parallel-p_star_w1)));
end

p_star = p_star_parallel + w_star;

Bp = gamma_k.*p_star + S_or_Psi * (invM*(S_or_Psi'*p_star));

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


