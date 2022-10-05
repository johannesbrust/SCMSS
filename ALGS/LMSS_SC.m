function [ xk, gk, fk, out ] = LMSS_SC( x, func, grad, pars )
%LMSS_SC: L-Multipoint-Symmetric-Secant trust-region implementation
% of the MSSM matrix
%
% Problem:
%
%   minimize f(x),
%
% where f:R^n -> R. Steps are computed using a L-MSS trust-region method.
% The L-MSS quasi-Newton compact representation with an initialization
%   that satisfies B0*Sk = gm*Sk (where gm=gamma a scalar) is:
%
% Bk    = B0    + Psi M Psi'
%       = B0    + [Sk Yk] [-gm*W - W (T + D + T') W, W; W, 0 ] [Sk Yk]',
%
% W     = inv(Sk'Sk)
%
% where Sk = [s_{k-1},...,s0]*Fk, Yk = [y_{k-1},...,y0] and
% s_{k-1} = xk - x_{k-1}, y_{k-1} = f'(xk) - f'(x_{k-1}). Here
% L = tril(Sk'Yk,-1)*Fk, T = triu(Sk'Yk,1)*Fk,  D = diag(Sk'Yk)*Fk, and
% Fk = diag(1/norm(s_{k-1}),...,1/norm(s0)). (Note that this implementation
%   normalizes the columns of Sk, i.e., uses a matrix like Fk)
% The trust-region subproblems are solved by
%
%   sc_mssm_infty.m          or
%   sc_mssm_2.m              or
%   obs.m                   or
%   lstrs.m                 or
%   st_truncated_CG_qn.m,
%
% of which the first 2 are shape-changing norms, and the remaining three
% are external solvers.
% This is a limited memory method, and the parameter m is used
% to limit the storage size of Psi to nx2m, which defines the
% sizes of the matrices that are derived from this quantity (i.e., M,
%       Psi'*Psi etc.)
%
% The trust-region updating strategy is used from
% J. Nocedal & S. J. Wright, Numerical Optimization, 1999, Algorithm 8.2.
%
% INPUTS:
% x     := Initial point
% func  := Objective function; f = func(x)
% grad  := Gradient function; g = grad(x)
% pars  := Struct with parameters
%   pars.tol    := Tolerance; Stop if norm(gk,'inf') < tol
%   pars.maxiter:= Maximum iterations
%   pars.print  := Flag to print iteration outputs
%   pars.m          := Memory parameter
%   pars.TRtol  := Trust-region subproblem tolerance
%   pars.TRmaxit:= Trust-region subproblem max. iterations
%   pars.c1     := TR acceptance. If pred/ared > c1 update xk
%   pars.c2     := Good step. If pred/ared > c2 'good step'
%   pars.c3     := Closeness to boundary for 'good step',
%                   if norm(sk) <= c3*Deltak keep radius, otherwise
%                   increase
%   pars.c4     := Increase radius by factor c4
%   pars.c5     := Lower bound for 'ok step', such that radius is not
%                   increased.
%   pars.c6     := Upper bound for 'ok step', such that radius is not
%                   increased.
%   pars.c7     := Trust-region shrinkage for 'not good step'
%   pars.whichInit     := Quasi-Newton initialization B_0 = gammak.*I
%       whichInit = 1 (all ones)
%       whichInit = 2 (constant specified value
%           whichInitCons = 1 (Constant initialization within algorithm.
%                           It uses gammak = min(y0'*y0/s0'*y0,1e4) and
%                                   gammak = max(gammak,1)
%       whichInit = 3 (gammak=yk'*yk / sk'*yk if positive, else no change)
%       whichInit = 4 (gammak=max(yk'*yk / sk'*yk,gams), where gams stores
%                       "q" previous "gammak" values
%           pars.q = (If whichInit = 4), stores the previous "q" gamma
%           values.
%   pars.whichSub := Which TR subproblem solver
%       whichSub = 1 (sc_mssm_infty)
%       whichSub = 2 (sc_mssm_2)
%       whichSub = 3 (l2_mssm)
%       whichSub = 4 (lstrs), "additional parameters in pars.LSTRS when
%                       this solver is selected"
%       whichSub = 5 (truncated CG), "additional parameters in pars.CG when
%                       this solver is selected"
%   pars.initDelta := Initial trust-region radius
%   pars.isSW      := Flag for representations
%                       isSW = 1 for [SkW,Yk]
%                       isSW = 0 for [Sk,Yk]
%   pars.whichDenseInit := Which dense initialization
%       whichDenseInit  := 1 (gammak_p = gammak (both parameters same))
%       whichDenseInit  := 2 (gammak_p = 0.5*gammak  0.5*max_k(gammak))
%       whichDenseInit  := 3 (gammak_p = 0.5*gammak  0.5*(gammakq))

% OUTPUTS:
% xk := Current iterate
% gk := Current gradient
% fk := Current function value
% out:= Struct with outputs
%   out.numiter     := Number of iterations
%   out.numf        := Number of function evaluations
%   out.numg        := Number of gradient evaluations
%   out.ng          := norm(gk,'inf')
%   out.ex          := Exit condition
%   out.Deltak      := Trust-region radius
%   out.numAccept   := Number step accepted
%   out.numTRInc    := Number trust-region radius increase
%   out.numTRDec    := Number trust-region radius decrease
%   out.ctime       := Computational time
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
% 07/16/20, J.B., Initial version
% 08/26/20, J.B., Modification to normalize "skip condition"
%                   abs(si'(yi-Bi*si)) > MSSMtol * norm(si) * norm(yi-Bi*si)
% 02/09/21, J.B., Including of deltaTol safeguard
% 02/10/21, J.B., Addition of LSTRS as subproblem solver,
%                 Addition of trucated CG as a subproblem solver
% 02/17/21, J.B., Storing Psi'*Psi when desired
% 02/24/21, J.B., Branching for initialization strategy and implementation
%                   of efficient updating strategies
% 02/25/21, J.B., Limited memory updating based on the branched
% initialization
% 03/02/21, J.B., Initialization strategies for gammak; 3,4
% 03/04/21, J.B., Constant initialization strategy
% 04/08/21, J.B., Modified skipping condition max(nsk*npsik,MSSMtolAdj)
% 04/27/21, J.B., Preparation for release
% 06/02/22, J.B., Preparation for a MSS method
% 06/07/22, J.B., Development of the compact representation by including
%                   vectors that store the norms of sk and yk
%                   Use an inverse memory index to represent the reverse
%                   order of vectors
% 06/08/22, J.B., Limited-memory updating
% 06/09/22, J.B., Further limited-memory updating
%                   Including option for [Sk*W Yk] representation,
%                       i.e, flag "isSY" added
% 06/10/22, J.B., Preparation for release
% 08/04/22, J.B., Including a dense initialization
% 08/18/22, J.B., Additional option for a dense initialization

% Parameter initializations
if isfield(pars,'tol')
    tol = pars.tol;
else
    tol = 1e-6;
end
if isfield(pars,'maxiter')
    maxiter = pars.maxiter;
else
    maxiter = 10000;
end
if isfield(pars,'print')
    print = pars.print;
else
    print = 0;
end
if isfield(pars,'m')
    m = pars.m;
else
    m = 8;
end
if isfield(pars,'c1')
    c1 = pars.c1;
else
    c1 = 1e-5;
end
if isfield(pars,'c2')
    c2 = pars.c2;
else
    c2 = 0.75;
end
if isfield(pars,'c3')
    c3 = pars.c3;
else
    c3 = 0.8;
end
if isfield(pars,'c4')
    c4 = pars.c4;
else
    c4 = 2;
end
if isfield(pars,'c5')
    c5 = pars.c5;
else
    c5 = 0.1;
end
if isfield(pars,'c6')
    c6 = pars.c6;
else
    c6 = 0.75;
end
if isfield(pars,'c7')
    c7 = pars.c7;
else
    c7 = 0.5;
end
if isfield(pars,'whichInit')
    gammak  = 1.0;
    whichInit = pars.whichInit;
    % Initialization 4 stores "q" previous gamma values
    if whichInit==4
        if isfield(pars,'q')
            q = pars.q;
        else
            q = m;
        end
        gams = zeros(q,1);
        numq = 0;
        mq = mod(numq,q);
        if mq==0; gams(1,1) = gammak; else gams(mq+1,1) = gammak; end;
        numq = numq+1;
    end
    if whichInit == 2
        if isfield(pars,'gammaInit')
            gammak = pars.gammaInit;
        else
            gammak = 1.0;
        end
        if isfield(pars,'whichInitCons')
            whichInitCons = pars.whichInitCons;
        else
            whichInitCons = 0;
        end
    end
else
    gammak = 1.0;
    whichInit = 1;
end
if isfield(pars,'whichDenseInit')
    whichDenseInit = pars.whichDenseInit;
    gammak_p       = gammak;
else
    whichDenseInit = 1;
    gammak_p       = gammak;
end

% Line search parameters
if isfield(pars,'LSftol')
    LSftol = pars.LSftol;
else
    LSftol = 1e-4; % Sufficient decrease f(xk+alpha.sk) <= fk + ftol.*alpha.*gk'*sk
end
if isfield(pars,'LSgtol')
    LSgtol = pars.LSgtol;
else
    LSgtol = 0.9; % Curvature abs(gk1'*sk) <= gtol.*abs(gk'*sk)
end
if isfield(pars,'LSxtol')
    LSxtol = pars.LSxtol;
else
    LSxtol = 1e-6; % Parameter for length of line search uncertainty interval
end
if isfield(pars,'LSstpmin')
    LSstpmin = pars.LSstpmin;
else
    LSstpmin = 0;
end
if isfield(pars,'LSstpmax')
    LSstpmax = pars.LSstpmax;
else
    LSstpmax = 5.0;
end
if isfield(pars,'LSmaxfev')
    LSmaxfev = pars.LSmaxfev;
else
    LSmaxfev = 3000;
end
% Subproblem solver
if isfield(pars,'whichSub')
    whichSub = pars.whichSub;
else
    whichSub = 1;
end
% Trust-region radius safeguard
if isfield(pars,'deltaTol')
    deltaTol = pars.deltaTol;
else
    deltaTol = 1e-18;
end
% Parameters for LSTRS
if whichSub == 4
    hpar = struct([]);
    Hfunc = @MultiplyBvec_lstrs_mssm;
    EigFunc = @eigensolver_lstrs;
    if isfield(pars,'LSTRS')
        if isfield(pars.LSTRS,'eigensolverpar')
            eigensolverpar = pars.LSTRS.eigensolverpar;
        else
            eigensolverpar.K = 2;
        end
        if isfield(pars.LSTRS,'epsilon')
            epsilon = pars.LSTRS.epsilon;
            if ~isfield(epsilon,'Delta')
                epsilon.Delta = 1e-4;
            end
        else
            epsilon.Delta = 1e-4;
        end
        if isfield(pars.LSTRS,'lopts')
            lopts = pars.LSTRS.lopts;
            if ~isfield(lopts,'message_level')
                lopts.message_level = 0;
            end
        else
            lopts.message_level = 0;
        end
    else
        eigensolverpar.K = 2;
        epsilon.Delta = 1e-4;
        lopts.message_level = 0;
    end
end
% Parameters for truncated CG
if whichSub == 5
    if isfield(pars,'CG')
        parsCG = pars.CG;
        if isfield(parsCG,'tolCG')
            tolCG = pars.CG.tolCG;
        else
            tolCG = 1e-4;
        end
        if isfield(parsCG,'maxitCG')
            maxitCG = pars.CG.maxitCG;
        else
            maxitCG = 100;
        end
        if isfield(parsCG,'show')
            showCG = pars.CG.show;
        else
            showCG = 0;
        end
    else
        tolCG = 1e-4;
        maxitCG = 100;
        showCG = 0;
    end
end
if isfield(pars,'initDelta')
    Deltak = pars.initDelta;
else
    Deltak = 1e5;
end
if isfield(pars,'isSW')
    isSW = pars.isSW;
else
    isSW = 0;
end

ts      = tic; % Measure time

% Preallocating counters
n       = length(x);
k       = 0;
numf    = 0;
numg    = 0;
numskip = 0;
numTRit = 0;

numAccept = 0;
numTRInc = 0;
numTRDec = 0;
numMem = 0;
linInd = 0;

% L-MSS initializations
% This implementation updates Sk, Yk, SkSk, SkYk, DSk, DYk, PsiPsi
% and computes M, W, Tk
% An option is included to use a representation with [Sk*W Yk] in
% which case the updates to M and PsiPsi change

m2   = 2*m;
M    = zeros(m2,m2);
Im   = eye(m);

% Intermediate buffers
b       = zeros(m2,1);

Sk = zeros(n,m);
Yk = zeros(n,m);
Dk = zeros(m,1);
Lk = zeros(m,m);
Tk = zeros(m,m);
SkSk = zeros(m,m);
YkYk = zeros(m,m);
SkYk = zeros(m,m);
W    = zeros(m,m);
% Intermediate buffers
b1 = zeros(m,1);
b2 = zeros(m,1);

DSk = zeros(m,1);
DYk = zeros(m,1);

PsiPsi = zeros(m2,m2);

midx    = 1:m; % Memory index
imidx   = zeros(m,1); % Inverse memory index

xk      = x;
fk      = func(xk);
gk      = grad(xk);
numf    = numf + 1;
numg    = numg + 1;
ng      = norm(gk,'inf');

if print == 1
    
    fprintf('----------- Running Algorithm ----------- \n'); % 41 chars
    fprintf(' L-MSS Trust-Region                      \n');
    fprintf('Iter \t fk      \t norm(gk)  \t TR     \t numf  \n');
    fprintf('%i \t %.4e \t %.4e \t %.4e \t %i  \n',k,fk,ng,Deltak,numf);
    
end

if ng < tol
    
    out.numiter = k;
    out.numf    = numf;
    out.numg    = numg;
    out.ng      = ng;
    out.ex      = 1;
    out.Deltak  = Deltak;
    out.numTRit = numTRit;
    
    out.numAccept = numAccept;
    out.numTRInc = numTRInc;
    out.numTRDec = numTRDec;
    out.ctime    = toc(ts);
    
    % Initial point satisfactory
    return;
    
end

% Line search in the steepest descent direction.
% Enables computing, sk and yk at k=0.
sk = -gk;
alpha0 = 1;

[xk1,fk1,gk1,~,info,nfev] = cvsrch_INTF(func,grad,n,xk,fk,...
    gk,sk,alpha0,LSftol,LSgtol,LSxtol,LSstpmin,LSstpmax,LSmaxfev);

numf = numf + nfev;
numg = numg + nfev;

% Computations for the L-MSS matrix
sk  = xk1-xk;
nsk = norm(sk);
sk1 = sk/nsk;
yk  = gk1 - gk;
%ss  = sk'*sk;
yy  = yk'*yk;
sy  = sk'*yk;
s1y = sk1'*yk;

% Setting-up of values for the initialization
switch whichInit
    case 2
        % Constant initialization based on initial s and y
        switch whichInitCons
            case 1
                gammak = min(yy/sy,1e4);
                gammak = max(gammak,1);
        end
    case 3
        if sy > 0; gammak = yy/sy; end
    case 4
        gammak = max(max(yy/sy,gams(1:min(numq,q),1)));
        mq = mod(numq,q);
        if mq==0; gams(1,1) = gammak; else gams(mq+1,1) = gammak; end;
        numq = numq+1;
end

% Dense initialization parameter
switch whichDenseInit
    case 1
        gammak_p = gammak;
    case 2
        maxgam   = gammak;
        gammak_p = 0.5*maxgam + 0.5*gammak;
    case 4
        gammak_p = 1;
end

% Container to store initializations
if whichDenseInit > 1
    gamIn = [gammak,gammak_p];
else
    gamIn = gammak;
end

% Updates to L-MSS compact representation

Sk(:,k+1) = sk1;
Yk(:,k+1) = yk;

% Scalings
Dk(k+1) = s1y/nsk;
Tk(k+1,k+1) = s1y/nsk;

SkYk(k+1,k+1) = s1y;
SkSk(k+1,k+1) = 1;%ss;
YkYk(k+1,k+1) = yy;

% Updates to 2mx2m block matrices
PsiPsi(k+1,k+1) = 1;
PsiPsi(m+1,m+1) = yy;
PsiPsi(k+1,m+1) = s1y;
PsiPsi(m+1,k+1) = s1y;

W(k+1,k+1) = 1;

M(k+1,k+1) = -gammak - Dk(k+1);
M(k+1,m+1) = 1/nsk;
M(m+1,k+1) = 1/nsk;
M(m+1,m+1) = 0;

DSk(k+1)    = nsk;
DYk(k+1)    = sqrt(yy);

numMem = numMem + 1;
linInd = linInd + 1;
linIndX= 1:linInd;

imidx(1:numMem) = midx(numMem:-1:1);

xk      = xk1;
gk      = gk1;
Deltak  = 2*nsk;
fk      = fk1;

k       = k + 1;

ng      = norm(gk,'inf');
if print == 1
    fprintf('%i \t %.4e \t %.4e \t %.4e \t %i  \n',k,fk,ng,Deltak,numf);
end


% Compute initial secant error
% Used to check whether B is correctly computed to this point
%---------------------------- Compute errors ------------------------------
% linInd = 1;
% idxSelC = [1:linInd, (m+1):(m+linInd)];
% PsiS = [Sk(:,midx(linIndX)),Yk(:,midx(linIndX))];
%
%
% B = gammak*eye(n) + PsiS*(M(idxSelC,idxSelC)*PsiS');
%
% Bs = B*sk;
%
% err0 = norm(yk-Bs);
% %------------------------- End Compute errors ---------------------------

% Main loop
while((tol < ng) && (k <= maxiter-1) ...
        && info ~= 6)
    
    % Trust-region step direction using the methods in [1]
    %   whichSub == 1 : Shape-changning infinity norm
    %   whichSub == 2 : Shape-changing 2 norm
    %   whichSub == 3 : Classical 2 norm
    %   whichSub == 4 : LSTRS
    %   whichSub == 5 : Truncated CG
    switch whichSub
        case 1
            
            % Implementation enables two representations
            % isSW == 0 then [Sk Yk]
            % isSW ~= 0 then [Sk*W Yk]
            if isSW == 0
                
                idxSel      = [1:linInd, m + (1:linInd)];
                [sk,iExit,Bsd]  = sc_mssm_infty( gk, [Sk(:,imidx(linIndX)),...
                    Yk(:,imidx(linIndX))],M(idxSel,idxSel),...
                    gamIn, Deltak, 1, 0,...
                    PsiPsi(idxSel,idxSel));
            else
                
                idxSel      = [1:linInd, m + (1:numMem)];
                [sk,iExit,Bsd]  = sc_mssm_infty( gk, [Sk(:,imidx(linIndX))*...
                    W(1:linInd,1:linInd),Yk(:,imidx(1:numMem))],...
                    M(idxSel,idxSel),gamIn, Deltak, 1, 0,PsiPsi(idxSel,idxSel));
            end
            if iExit == 1; info = 6; end;
            
        case 2
            if isSW == 0
                
                idxSel  = [1:linInd, m + (1:linInd)];
                [sk,~,Bsd]  = sc_mssm_2( gk, [Sk(:,imidx(linIndX)),...
                    Yk(:,imidx(linIndX))],M(idxSel,idxSel),...
                    gamIn, Deltak, 1, 0,...
                    PsiPsi(idxSel,idxSel));
            else
                
                idxSel  = [1:linInd, m + (1:numMem)];
                [sk,~,Bsd]  = sc_mssm_2( gk, [Sk(:,imidx(linIndX))*...
                    W(1:linInd,1:linInd),Yk(:,imidx(1:numMem))],...
                    M(idxSel,idxSel),gamIn, Deltak, 1, 0,PsiPsi(idxSel,idxSel));
            end
            %if iExit == 1; info = 6; end;
        case 3
            iExit = 0;
            try
                if isSW == 0
                    
                    idxSel       = [1:linInd, m + (1:linInd)];
                    [~,sk,~,~,~] = obs_mssm( gk, [], [], Deltak, gammak,...
                        [Sk(:,imidx(linIndX)),Yk(:,imidx(linIndX))],...
                        M(idxSel,idxSel),PsiPsi(idxSel,idxSel));
                else
                    
                    idxSel       = [1:linInd, m + (1:numMem)];
                    [~,sk,~,~,~] = obs_mssm( gk, [], [], Deltak, gammak,...
                        [Sk(:,imidx(linIndX))*W(1:linInd,1:linInd),...
                        Yk(:,imidx(1:numMem))],M(idxSel,idxSel),...
                        PsiPsi(idxSel,idxSel));
                    
                end
            catch ME %#ok<NASGU>
                iExit = 1;
            end
            if iExit == 1; info = 6; end;
        case 4
            iExit = 0;
            try
                hpar(1).f = gammak;
                if isSW == 0
                    idxSel       = [1:linInd, m + (1:linInd)];
                    hpar(2).f    = [Sk(:,imidx(linIndX)),Yk(:,imidx(linIndX))];
                    %nsel         = 2*linInd;
                else
                    idxSel       = [1:linInd, m + (1:numMem)];
                    hpar(2).f    = [Sk(:,imidx(linIndX))*W(1:linInd,1:linInd),...
                                        Yk(:,imidx(1:numMem))];                    
                    %nsel         = linInd + numMem;                
                end
                
                hpar(3).f = M(idxSel,idxSel); %\Im(1:nsel,1:nsel);
                
                [sk,~,infoLS,~] = lstrs(Hfunc,gk,Deltak,epsilon,EigFunc,...
                    lopts,hpar,eigensolverpar);
            catch ME %#ok<NASGU>
                iExit = 1;
                infoLS = -4;
            end
            if iExit == 1; info = 6; end;
            if infoLS== -4; info = 6; sk = zeros(n,1); end
        case 5
            if isSW == 0
                idxSel      = [1:linInd, m + (1:linInd)];
                nsel        = 2*linInd;
                [sk,~,~,~]  = st_truncated_CG_qn_mssm(gk,[Sk(:,imidx(linIndX)),...
                    Yk(:,imidx(linIndX))],M(idxSel,idxSel),gammak,Deltak,nsel,tolCG,...
                    maxitCG,showCG);
            else
                idxSel       = [1:linInd, m + (1:numMem)];
                nsel         = linInd + numMem;
                [sk,~,~,~]  = st_truncated_CG_qn_mssm(gk,[Sk(:,imidx(linIndX))*W(1:linInd,1:linInd),...
                    Yk(:,imidx(1:numMem))],M(idxSel,idxSel),gammak,Deltak,nsel,tolCG,...
                    maxitCG,showCG);
            end
        otherwise
            idxSel          = [1:linInd, m + (1:linInd)];
                [sk,iExit]  = sc_mssm_infty( gk, [Sk(:,imidx(linIndX)),...
                    Yk(:,imidx(linIndX))],M(idxSel,idxSel),...
                    gammak, Deltak, 1, 0,...
                    PsiPsi(idxSel,idxSel));
            if iExit == 1; info = 6; end;
    end
    
    % Trial step
    xk1t        = xk+sk;
    
    gk1t        = grad(xk1t);
    numg        = numg + 1;
    
    yk          = gk1t - gk;
    
    fk1t        = func(xk1t);
    numf        = numf + 1;
    
    ared        = fk - fk1t;
    
    % Intermediate "buffers" that are reused later for
    % updating matrices
    b1(1:numMem,1) = Yk(:,midx(1:numMem))'*sk;
    b2(1:numMem,1) = Sk(:,midx(1:numMem))'*sk;
    
    % Computing a buffer "b" for the quadratic model
    % based on which option was selected 
    % "b" can be used to compute 
    % s*B*s = gammak*ss + b(1:bmem,1)'*(M(idxSel,idxSel)*b(1:bmem,1));
    if isSW == 0
        
        bmem = 2*linInd;
        idxSel          = [1:linInd, m + (1:linInd)];
        
        if numMem < m
            idxSelb = imidx(linIndX);
        else
            idxSelb = imidxL(linIndX);
        end
        
        b(1:bmem,1) = [b2(idxSelb,1);b1(idxSelb,1)];
        
    else
        
        bmem = linInd+numMem;
        idxSel = [(1:linInd), (m + (1:numMem))];
        
        if numMem < m
            idxSelb2 = imidx(linIndX);
            idxSelb1 = imidx(1:numMem);
        else
            idxSelb2 = imidxL(linIndX);
            idxSelb1 = imidxL(1:numMem);
        end
        
        b(1:bmem,1) = [W(1:linInd,1:linInd)*b2(idxSelb2,1);b1(idxSelb1,1)];
        
    end
    
    ss = sk'*sk;
    
    yy = yk'*yk;
    
    sy = sk'*yk;
    
    nsk = sqrt(ss);
    sk1 = sk/nsk;
    s1y = sk1'*yk;
    
    % If no dense initialization is used then s'Bs can be computed
    % efficiently
    if whichDenseInit > 1
        sBs = sk'*Bsd;        
    else
        sBs = gammak*ss + b(1:bmem,1)'*(M(idxSel,idxSel)*b(1:bmem,1));
    end
    
    gs = gk'*sk;
    
    switch whichInit
        case 3
            if sy > 0; gammak = yy/sy; end
        case 4
            gammak = max(max(yy/sy,gams(1:min(numq,q),1)));
            mq = mod(numq,q);
            %if mq==0; gams(1,1) = gammak; else gams(mq+1,1) = gammak; end;
            if mq==0; gams(1,1) = yy/sy; else gams(mq+1,1) = yy/sy; end;
            numq = numq+1;
    end
    
    
    % Dense initialization parameter
    switch whichDenseInit
        case 1
            gammak_p = gammak;
        case 2
            if gammak > maxgam   
                maxgam = gammak;
            end
            gammak_p = 0.5*maxgam + 0.5*gammak;
        case 3
            if sy > 0; maxgam = yy/sy; else maxgam = gammak; end;
            gammak_p = 0.5*maxgam + 0.5*gammak;
        case 4
            if sy > 0; gammak_p = yy/sy; end;
            %gammak_p = 0.5*maxgam + 0.5*gammak;
    end
    
    % Including a dense initialization
    % For input in the TR subproblem solvers
    if whichDenseInit > 1
        gamIn(:) = [gammak,gammak_p];
    else
        gamIn    = gammak;
    end
    
    pred        = -(gs + 0.5*sBs);
    
    rho         = ared/pred;
    
    % Trust-region radius updates
    if c1 < rho
        xk1 = xk1t;
        gk1 = gk1t;
        fk1 = fk1t;
        
        numAccept = numAccept + 1;
    else
        xk1 = xk;
    end
    
    nsk = sqrt(ss);
    
    if c2 < rho
        if nsk <= c3*Deltak
            Deltak = 1*Deltak;
        else
            Deltak = c4*Deltak;
            numTRInc = numTRInc + 1;
        end
    elseif (c5 <= rho) && (rho <= c6)
        Deltak = 1*Deltak;
    else
        Deltak = c7*Deltak;
        numTRDec = numTRDec + 1;
    end
    
    % Update limited memory quasi-Newton vectors/matrices
    if(numMem < m)
        
        numMem = numMem + 1;
        
        imidx(1:numMem) = midx(numMem:-1:1);
        
        if numMem == m
            
            % Full inverse index
            % Because of memory access patterns
            imidxL          = numMem:-1:1;
            
        end
        
        numMemM1 = numMem-1;
        Sk(:,midx(numMem)) = sk1;
        Yk(:,midx(numMem)) = yk;
        Dk(numMem) = s1y/nsk;
        
        
        DSk(numMem)    = nsk;
        DYk(numMem)    = sqrt(yy);
        
        
        Lk(numMem,1:numMemM1) = (b1(1:numMemM1,1)/nsk); %./DSk(imidx(1:numMemM1));
        
        Skyk = Sk(:,midx(1:numMem))'*yk;
        
        %Tk(1:numMem,numMem) = Skyk/nsk;
        
        SkSk(1:numMemM1,numMem) = b2(1:numMemM1,1)/nsk;
        SkSk(numMem,numMem) = 1;%ss;
        SkSk(numMem,1:numMemM1) = SkSk(1:numMemM1,numMem); %b2(1:numMemM1,1)/nsk;
        YkYk(1:numMem,numMem) = Yk(:,midx(1:numMem))'*yk;
        YkYk(numMem,1:numMemM1) = YkYk(1:numMemM1,numMem);
        
        SkYk(1:numMem,numMem) = Skyk;
        SkYk(numMem,1:numMemM1) = b1(1:numMemM1)/nsk;
        
        %Tk(numMem,1:numMemM1) = (b1(1:numMemM1,1)/nsk)./DSk(1:numMemM1); % SkYk(numMem,1:numMem) /nsk; %(Skyk/nsk)./DSk(1:numMem);
        
        Tk(1:numMem,1:numMem) = tril(SkYk(1:numMem,1:numMem)*diag(1./DSk(1:numMem)));
        
        % LDL' factorization of Sk'Sk
        
        [Ldl,dR,Pldl] = ldl(SkSk(imidx(1:numMem),imidx(1:numMem)),'vector');  %ldl decomposition of PsiPsi'
        dR = diag(dR);  %vectorize Dldl
        linIndX1 = find(dR>(sqrt(eps)*max(dR)));  %generate mask to find lin indep columns of Psi
        linIndX = Pldl(linIndX1);
        linInd = length(linIndX);
        
        R(1:linInd,1:linInd) = diag(sqrt(dR(linIndX1)))*Ldl(linIndX1,linIndX1)';
        
        W(1:linInd,1:linInd) =  ...
            R(1:linInd,1:linInd)\(R(1:linInd,1:linInd)'\Im(1:linInd,1:linInd));
        
        
        % Computing the blocks in M
        % Depending on the representation
        if isSW == 0
            
            M(1:linInd,1:linInd) = -gammak*W(1:linInd,1:linInd) - ...
                W(1:linInd,1:linInd) * (( Tk(imidx(linIndX),imidx(linIndX)) + ...
                Tk(imidx(linIndX),imidx(linIndX))' - diag(Dk(imidx(linIndX))) )*W(1:linInd,1:linInd));
            
            M(1:linInd,(m+1):(m+linInd)) = W(1:linInd,1:linInd)*...
                diag(1./DSk(imidx(linIndX)));
            M((m+1):(m+linInd),1:linInd) = diag(1./DSk(imidx(linIndX)))*...
                W(1:linInd,1:linInd);
            
            PsiPsi(1:linInd,1:linInd) = SkSk(imidx(linIndX),imidx(linIndX));
            PsiPsi(1:linInd,(m+1):(m+linInd)) = SkYk(imidx(linIndX),imidx(linIndX));
            PsiPsi((m+1):(m+linInd),1:linInd) = (SkYk(imidx(linIndX),imidx(linIndX)))';
            PsiPsi((m+1):(m+linInd),(m+1):(m+linInd)) = (YkYk(imidx(linIndX),imidx(linIndX)));
            
        else
            
            M(1:linInd,1:linInd) = -gammak*SkSk(imidx(linIndX),imidx(linIndX)) - ...
                ( Tk(imidx(linIndX),imidx(linIndX)) + ...
                Tk(imidx(linIndX),imidx(linIndX))' - diag(Dk(imidx(linIndX))) );
            
            ddsk = diag(1./DSk(imidx(1:numMem)));
            
            M(1:linInd,(m+1):(m+numMem)) = ddsk(linIndX,1:numMem);
            
            M((m+1):(m+numMem),1:linInd) = ddsk(1:numMem,linIndX);
            
            PsiPsi(1:linInd,1:linInd) = W(1:linInd,1:linInd);
            PsiPsi(1:linInd,(m+1):(m+numMem)) = W(1:linInd,1:linInd)*SkYk(imidx(linIndX),imidx(1:numMem));
            PsiPsi((m+1):(m+numMem),1:linInd) = PsiPsi(1:linInd,(m+1):(m+numMem))'; %(SkYk(imidx(linIndX),imidx(linIndX)))';
            PsiPsi((m+1):(m+numMem),(m+1):(m+numMem)) = (YkYk(imidx(1:numMem),imidx(1:numMem)));
            
        end
        
%---------------------------- Compute errors ------------------------------
% This is for testing purposes only and intended for verifying
% the secant condition and that updates like Sk'*Sk are valid
        
        %             %imidx(1:numMem) = midx(numMem:-1:1);
        %
        %             %idxSelC = [imidx(1:linInd), m+imidx(1:linInd)];
        %             if isSW == 0
        %                 idxSelC = [1:linInd, m+(1:linInd)];
        %                 idxSelY = imidx(linIndX);
        %                  PsiS = [Sk(:,imidx(linIndX)),Yk(:,idxSelY)];
        %             else
        %
        %                 idxSelC = [1:linInd, m+(1:numMem)];
        %                 idxSelY = imidx(1:numMem);
        %                  PsiS = [Sk(:,imidx(linIndX))*W(1:linInd,1:linInd),Yk(:,idxSelY)];
        %             end
        %
        %
        %     %             idxSelC = [1:linInd, (m+1):(m+linInd)];
        %     %             PsiS = [Sk(:,midx(linIndX)),Yk(:,midx(linIndX))];
        %
        %                 PSPS = PsiS'*PsiS;
        %
        %                 SY = Sk(:,imidx(1:numMem))'*Yk(:,imidx(1:numMem));
        %                 SS = Sk(:,imidx(1:numMem))'*Sk(:,imidx(1:numMem));
        %                 YY = Yk(:,imidx(1:numMem))'*Yk(:,imidx(1:numMem));
        %
        %
        %                 % Unscaled quantities/Scaled
        %     %             Fk = diag(1./DSk(imidx(linIndX)));
        %     %             Fki = diag(DSk(imidx(linIndX)));
        %     %             Fki = eye(linInd);
        %     %
        %     %             SkU = Sk(:,imidx(midx(linIndX)))*Fki;
        %     %             PsiU = [SkU,Yk(:,imidx(midx(linIndX)))];
        %     %
        %     %             SYU = SkU'*Yk(:,imidx(midx(1:numMem)));
        %     %             SYU = SkU'*(Yk(:,imidx(midx(1:numMem)))*Fk);
        %     %
        %     %             SSU = SkU'*SkU;
        %     %             WU  = SSU\eye(numMem);
        %     %             TU = triu(SYU);
        %     %             DU = diag(TU);
        %     %             MU = [-gammak*WU - WU*(TU+TU'-diag(DU))*WU, WU*Fk;...
        %     %                 Fk*WU, zeros(numMem,numMem)];
        %     %
        %     %             % Direct (unscaled) formula of B
        %     %
        %     %             BU = gammak*eye(n) + PsiU*(MU*PsiU');
        %     %
        %     %             Bsu = BU*sk;
        %     %
        %     %             err0U = norm(yk-Bsu);
        %
        %                 B = gammak*eye(n) + PsiS*(M(idxSelC,idxSelC)*PsiS');
        %
        %                 Bs = B*sk;
        %
        %                 err0 = norm(yk-Bs);
        %
        %                 err1 = norm(SkSk(imidx(1:numMem),imidx(1:numMem))-SS,'fro');
        %                 err2 = norm(SkYk(imidx(1:numMem),imidx(1:numMem))-SY,'fro');
        %                 err3 = norm(YkYk(imidx(1:numMem),imidx(1:numMem))-YY,'fro');
        %
        %                 err4 = norm(PSPS-PsiPsi(idxSelC,idxSelC),'fro');
        %
        %                 err5 = norm(SkSk(imidx(linIndX),imidx(linIndX))*W(1:linInd,1:linInd)-...
        %                     Im(1:linInd,1:linInd),'fro');
        %
        %                 dbg = 1;
%
%---------------------------- End Compute errors --------------------------
        
    elseif(numMem == m)
        
        numMemM1 = numMem-1;
        midxs = midx(1);
        
        Sk(:,midxs) = sk1;
        Yk(:,midxs) = yk;
        
        midx(1:(end-1)) = midx(2:end);
        midx(end) = midxs;
        
        imidx(1:numMem) = midx(numMem:-1:1);
        
        % Copy previous (small) arrays
        Dk(1:numMemM1) = Dk(2:numMem);
        Lk(1:numMemM1,1:numMemM1) = Lk(2:numMem,2:numMem);
        Tk(1:numMemM1,1:numMemM1) = Tk(2:numMem,2:numMem);
        SkSk(1:numMemM1,1:numMemM1) = SkSk(2:numMem,2:numMem);
        YkYk(1:numMemM1,1:numMemM1) = YkYk(2:numMem,2:numMem);
        
        SkYk(1:numMemM1,1:numMemM1) = SkYk(2:numMem,2:numMem);
        
        DSk(1:numMemM1) = DSk(2:numMem);
        DYk(1:numMemM1) = DYk(2:numMem);
        
        % Updates
        % First temporary "b" buffers
        b1(1:numMemM1) = b1(2:numMem);
        b1(numMem) = s1y;
        b2(1:numMemM1) = b2(2:numMem);
        b2(numMem) = ss/nsk;
        
        Dk(numMem) = s1y/nsk;
        
        DSk(numMem) = nsk;
        DYk(numMem) = sqrt(yy);
        
        Lk(numMem,1:numMemM1) = b1(1:numMemM1,1);
        %Tk(1:numMem,numMem) = Sk(:,midx(1:numMem))'*yk;
        
        Skyk = Sk(:,midx(1:numMem))'*yk;
        
        SkSk(1:numMemM1,numMem) = b2(1:numMemM1,1)/nsk;
        SkSk(numMem,numMem) = 1;
        SkSk(numMem,1:numMemM1) = SkSk(1:numMemM1,numMem);
        
        YkYk(1:numMem,numMem) = Yk(:,midx(1:numMem))'*yk;
        %YkYk(numMem,numMem) = yy;
        YkYk(numMem,1:numMemM1) = YkYk(1:numMemM1,numMem);
        
        SkYk(1:numMem,numMem) = Skyk;
        SkYk(numMem,1:numMemM1) = b1(1:numMemM1)/nsk;
        
        Tk(1:numMem,1:numMem) = tril(SkYk(1:numMem,1:numMem)*diag(1./DSk(1:numMem)));
        
        %% Update SS and M
        %[Ldl,dR,Pldl] = ldl(SkSk(imidx(1:numMem),imidx(1:numMem)),'vector');  %ldl decomposition of PsiPsi'
        [Ldl,dR,Pldl] = ldl(SkSk(imidxL(1:numMem),imidxL(1:numMem)),'vector');  %ldl decomposition of PsiPsi'
        dR = diag(dR);  %vectorize Dldl
        linIndX1 = find(dR>(sqrt(eps)*max(dR)));  %generate mask to find lin indep columns of Psi
        linIndX = Pldl(linIndX1);
        
        linInd = length(linIndX);
        
        R(1:linInd,1:linInd) = diag(sqrt(dR(linIndX1)))*Ldl(linIndX1,linIndX1)';
        
        W(1:linInd,1:linInd) =  ...
            R(1:linInd,1:linInd)\(R(1:linInd,1:linInd)'\Im(1:linInd,1:linInd));
        
        % Computing the blocks in M
        if isSW == 0
            
            M(1:linInd,1:linInd) = -gammak*W(1:linInd,1:linInd) - ...
                W(1:linInd,1:linInd) * (( Tk(imidxL(linIndX),imidxL(linIndX)) + ...
                Tk(imidxL(linIndX),imidxL(linIndX))' - diag(Dk(imidxL(linIndX))) )*W(1:linInd,1:linInd));
            
            M(1:linInd,(m+1):(m+linInd)) = W(1:linInd,1:linInd)*...
                diag(1./DSk(imidxL(linIndX)));
            M((m+1):(m+linInd),1:linInd) = diag(1./DSk(imidxL(linIndX)))*...
                W(1:linInd,1:linInd);
            
            PsiPsi(1:linInd,1:linInd) = SkSk(imidxL(linIndX),imidxL(linIndX));
            PsiPsi(1:linInd,(m+1):(m+linInd)) = SkYk(imidxL(linIndX),imidxL(linIndX));
            PsiPsi((m+1):(m+linInd),1:linInd) = (SkYk(imidxL(linIndX),imidxL(linIndX)))';
            PsiPsi((m+1):(m+linInd),(m+1):(m+linInd)) = (YkYk(imidxL(linIndX),imidxL(linIndX)));
            
        else
            
            M(1:linInd,1:linInd) = -gammak*SkSk(imidxL(linIndX),imidxL(linIndX)) - ...
                ( Tk(imidxL(linIndX),imidxL(linIndX)) + ...
                Tk(imidxL(linIndX),imidxL(linIndX))' - diag(Dk(imidxL(linIndX))) );
            
            ddsk = diag(1./DSk(imidxL(1:numMem)));
            
            M(1:linInd,(m+1):(m+numMem)) = ddsk(linIndX,1:numMem);
            
            M((m+1):(m+numMem),1:linInd) = ddsk(1:numMem,linIndX);
            
            PsiPsi(1:linInd,1:linInd) = W(1:linInd,1:linInd);
            PsiPsi(1:linInd,(m+1):(m+numMem)) = W(1:linInd,1:linInd)*SkYk(imidxL(linIndX),imidxL(1:numMem));
            PsiPsi((m+1):(m+numMem),1:linInd) = PsiPsi(1:linInd,(m+1):(m+numMem))'; %(SkYk(imidx(linIndX),imidx(linIndX)))';
            PsiPsi((m+1):(m+numMem),(m+1):(m+numMem)) = (YkYk(imidxL(1:numMem),imidxL(1:numMem)));
            
        end
        
%---------------------------- Compute errors ------------------------------
% Checking errors when the limited memory strategy is used

        %             %imidx(1:numMem) = midx(numMem:-1:1);
        %
        %             %idxSelC = [imidx(1:linInd), m+imidx(1:linInd)];
        %             if isSW == 0
        %                 idxSelC = [1:linInd, m+(1:linInd)];
        %                 PsiS = [Sk(:,imidx(linIndX)),Yk(:,imidx(linIndX))];
        %             else
        %                 idxSelC = [1:linInd, m+(1:numMem)];
        %                 PsiS = [Sk(:,imidx(linIndX))*W(1:linInd,1:linInd),...
        %                     Yk(:,imidx(1:numMem))];
        %             end
        %
        % %             idxSelC = [1:linInd, (m+1):(m+linInd)];
        % %             PsiS = [Sk(:,midx(linIndX)),Yk(:,midx(linIndX))];
        %
        %             PSPS = PsiS'*PsiS;
        %
        %             SY = Sk(:,imidx(1:numMem))'*Yk(:,imidx(1:numMem));
        %             SS = Sk(:,imidx(1:numMem))'*Sk(:,imidx(1:numMem));
        %             YY = Yk(:,imidx(1:numMem))'*Yk(:,imidx(1:numMem));
        %
        %
        %             B = gammak*eye(n) + PsiS*(M(idxSelC,idxSelC)*PsiS');
        %
        %             Bs = B*sk;
        %
        %             %sk11 = Sk(:,midx(1));
        %             %yk11 = Yk(:,midx(1));
        %
        %             %err00= norm(yk11-B*(sk11*DSk(1)));
        %
        %             err0 = norm(yk-Bs);
        %
        %             err1 = norm(SkSk(imidxL(1:numMem),imidxL(1:numMem))-SS,'fro');
        %             err2 = norm(SkYk(imidxL(1:numMem),imidxL(1:numMem))-SY,'fro');
        %             err3 = norm(YkYk(imidxL(1:numMem),imidxL(1:numMem))-YY,'fro');
        %
        %             err4 = norm(PSPS-PsiPsi(idxSelC,idxSelC),'fro');
        %
        %             err5 = norm(SkSk(imidxL(linIndX),imidxL(linIndX))*W(1:linInd,1:linInd)-...
        %                 Im(1:linInd,1:linInd),'fro');
        %
        %             dbg = 1;
        
%---------------------------- Compute errors ------------------------------        

    end
    
    k = k + 1;
    
    % Updated iterates
    xk       = xk1;
    gk       = gk1;
    fk       = fk1;
    
    ng       = norm(gk,'inf');
    
    if print == 1
        fprintf('%i \t %.4e \t %.4e \t %.4e \t %i  \n',k,fk,ng,Deltak,numf);
    end
    
    % Safeguard for vanishing radius
    if Deltak < deltaTol
        info = 6;
    end
    
end

out.numiter = k;
out.numf    = numf;
out.numg    = numg;
out.ng      = ng;
if k == maxiter && tol < ng; info = 6; end;
out.ex      = info;
out.numskip = numskip;
out.numTRit = numTRit;

out.numAccept = numAccept;
out.numTRInc = numTRInc;
out.numTRDec = numTRDec;
out.ctime    = toc(ts);
% out.numinfo2 = numinfo2;
% out.numinfo3 = numinfo3;
% out.numinfo4 = numinfo4;
% out.numinfo5 = numinfo5;

return

end

