function [ xk, gk, fk, out ] = LSR1_SC( x, func, grad, pars )
%LSR1_SC: L-SR1 trust-region implementation using the subproblem solvers
% from article SC-SR1: MATLAB Software for Limited-Memory SR1
% Trust-Region Methods, (Brust,Burdakov,Erway,Marcia,2018-2021)         [1]
%
% The trust-region updating strategy is used from
% J. Nocedal & S. J. Wright, Numerical Optimization, 1999, Algorithm 8.2.
%
% Problem:
%
%   minimize f(x), 
%
% where f:R^n -> R. Steps are computed using a L-SR1 trust-region method. 
% The L-SR1 quasi-Newton compact representation is:
%
% Bk    = B0    + Psi inv(M) Psi'
%       = B0    + [Yk-B0*Sk] inv(L + L' + D - Sk'B0Sk) [Yk-B0*Sk]',
%       
% where Sk = [s0,...,s_{k-1}], Yk = [y0,...,y_{k-1}] and
% s_{k-1} = xk - x_{k-1}, y_{k-1} = f'(xk) - f'(x_{k-1}). Here
% L = tril(Sk'Yk,-1), D = diag(Sk'Yk). In this implementation, B0 =
% gamma*I. The trust-region subproblems are solved by 
%
%   sc_sr1_infty.m          or
%   sc_sr1_2.m              or
%   obs.m                   or
%   lstrs.m                 or
%   st_truncated_CG_qn.m,
%
% of which the first 2 are developed in [1], and the remaining three
% are external solvers.
% This is a limited memory method, and the parameter m is used
% to limit the storage size of Psi to nxm, which defines the
% sizes of the matrices that are derived from this quantity (i.e., L, D, etc)
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
%   pars.SR1tol := SR1 update tolerance: If norm(sk'*(ak)) > SR1tol then
%                   update QN matrix.
%   pars.SR1tolAdj := Adjustment to SR1 update condition:
%                       norm(sk'*(ak)) > SR1tol*max(norm(sk)*norm(ak),SR1tolAdj)
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
%       whichSub = 1 (sc_sr1_infty)
%       whichSub = 2 (sc_sr1_2)
%       whichSub = 3 (l2_sr1)
%       whichSub = 4 (lstrs), "additional parameters in pars.LSTRS when
%                       this solver is selected"
%       whichSub = 5 (truncated CG), "additional parameters in pars.CG when
%                       this solver is selected"
%   pars.initDelta := Initial trust-region radius

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
% 07/16/20, J.B., Initial version
% 08/26/20, J.B., Modification to normalize "skip condition" 
%                   abs(si'(yi-Bi*si)) > SR1tol * norm(si) * norm(yi-Bi*si)
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
% 04/08/21, J.B., Modified skipping condition max(nsk*npsik,SR1tolAdj)
% 04/27/21, J.B., Preparation for release

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
if isfield(pars,'SR1tol') 
    SR1tol = pars.SR1tol;
else
    SR1tol = 1e-8; 
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
    deltaTol = 1e-20; 
end
% Parameters for LSTRS
if whichSub == 4
    hpar = struct([]);
    Hfunc = @MultiplyBvec_lstrs;
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
% Flag to store matrix Psi'*Psi
if isfield(pars,'storePsiPsi') 
    storePsiPsi = pars.storePsiPsi;
else
    storePsiPsi = 0; 
end
if isfield(pars,'initDelta')
    Deltak = pars.initDelta;
else
    Deltak = 1e5;
end
if isfield(pars,'SR1tolAdj')
    SR1tolAdj = pars.SR1tolAdj;
else
    SR1tolAdj = 1e-2;
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

%SR1tolAdj = 1e-12;

% L-SR1 initializations
% This implementation updates invM and depending on input choices
% updates Psi (or Sk, Yk, Dk, Lk, SkSk, YkYk). The product Psi'*Psi may also 
% chosen to be updated.
invM    = zeros(m,m);

% Intermediate buffers
b       = zeros(m,1);

% Constant initialization
if whichInit == 1 || whichInit == 2
    Psi = zeros(n,m);  
else
    % Nonconstant initialization
    Sk = zeros(n,m);
    Yk = zeros(n,m);
    Dk = zeros(m,1);
    Lk = zeros(m,m);
    Tk = zeros(m,m);
    SkSk = zeros(m,m);
    YkYk = zeros(m,m);    
    % Intermediate buffers
    b1 = zeros(m,1);
    b2 = zeros(m,1);
end

if storePsiPsi == 1
    PsiPsi = zeros(m,m);
end

midx    = 1:m; % Memory index

xk      = x;
fk      = func(xk); 
gk      = grad(xk);
numf    = numf + 1;
numg    = numg + 1;
ng      = norm(gk,'inf');

if print == 1

    fprintf('----------- Running Algorithm ----------- \n'); % 41 chars
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

% Computations for the L-SR1 matrix
sk = xk1-xk;
yk = gk1 - gk;
ss = sk'*sk;
yy = yk'*yk;
sy = sk'*yk;

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
   
psik    = yk - gammak.*sk;
spsi    = sy - gammak.*ss;
%spsi    = sk'*psik;

%nsk     = norm(sk);
nsk     = sqrt(ss);
%psipsi  = psik'*psik;
psipsi  =  yy - 2*gammak*sy + gammak*gammak*ss;
npsik   = sqrt(psipsi);
%npsik   = norm(psik);

% Updates if L-SR1 matrix is well defined
if SR1tol < (abs(spsi)*nsk*npsik)    
    % Branched updating depending on initialization strategy
    if whichInit == 1 || whichInit == 2
        Psi(:,k+1)  = psik;
        invM(k+1,k+1) = spsi;    
        if storePsiPsi == 1
            PsiPsi(k+1,k+1) = psipsi;
        end
    else
        Sk(:,k+1) = sk;
        Yk(:,k+1) = yk;
        Dk(k+1) = sy;
        Tk(k+1,k+1) = sy;
        SkSk(k+1,k+1) = ss;
        YkYk(k+1,k+1) = yy;    
        invM(k+1,k+1) = spsi;
        if storePsiPsi == 1
            PsiPsi(k+1,k+1) = psipsi;
        end
    end    
    numMem = numMem + 1;
end

xk      = xk1;
gk      = gk1;
Deltak  = 2*nsk;
fk      = fk1;

k       = k + 1;

ng      = norm(gk,'inf');
if print == 1                
        fprintf('%i \t %.4e \t %.4e \t %.4e \t %i  \n',k,fk,ng,Deltak,numf);
end

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
            if storePsiPsi == 1
                if whichInit~=1&&whichInit~=2
                    [sk,iExit] = sc_sr1_infty( gk, Sk(:,midx(1:numMem)),...
                        Yk(:,midx(1:numMem)), gammak, Deltak, 1, 0,...
                        PsiPsi(1:numMem,1:numMem),invM(1:numMem,1:numMem));
                else
                    [sk,iExit] = sc_sr1_infty( gk, Psi(:,midx(1:numMem)),...
                        invM(1:numMem,1:numMem), gammak, Deltak, 1, 0, PsiPsi(1:numMem,1:numMem));
                end
            else                
                [sk,iExit] = sc_sr1_infty( gk, Psi(:,midx(1:numMem)),...
                    invM(1:numMem,1:numMem), gammak, Deltak, 1, 0);
            end
            if iExit == 1; info = 6; end;
        case 2
            if storePsiPsi == 1
                if whichInit~=1&&whichInit~=2
                    [sk,~] = sc_sr1_2( gk, Sk(:,midx(1:numMem)),...
                        Yk(:,midx(1:numMem)), gammak, Deltak, 1, 0,...
                        PsiPsi(1:numMem,1:numMem),invM(1:numMem,1:numMem));
                else
                    [sk,~] = sc_sr1_2( gk, Psi(:,midx(1:numMem)),...
                        invM(1:numMem,1:numMem), gammak, Deltak, 1, 0, PsiPsi(1:numMem,1:numMem));
                end
            else                
                [sk,~] = sc_sr1_2( gk, Psi(:,midx(1:numMem)),...
                    invM(1:numMem,1:numMem), gammak, Deltak, 1, 0);
            end
            %if iExit == 1; info = 6; end;
        case 3
            iExit = 0;
            try
                if storePsiPsi == 1
                    if whichInit~=1&&whichInit~=2
                        [~,sk,~,~,~] = obs( gk, Sk(:,midx(1:numMem)),...
                            Yk(:,midx(1:numMem)), Deltak, gammak,...
                        [],invM(1:numMem,1:numMem),PsiPsi(1:numMem,1:numMem));
                    else
                        [~,sk,~,~,~] = obs( gk, [], [], Deltak, gammak,...
                        Psi(:,midx(1:numMem)),invM(1:numMem,1:numMem));
                    end
                else
                    [~,sk,~,~,~] = obs( gk, [], [], Deltak, gammak,...
                        Psi(:,midx(1:numMem)),invM(1:numMem,1:numMem));
                end
            catch ME %#ok<NASGU>
                iExit = 1;
            end
            if iExit == 1; info = 6; end;
        case 4
            iExit = 0;
            try
                if whichInit~=1&&whichInit~=2
                    hpar(1).f = gammak;
                    hpar(2).f = Sk(:,midx(1:numMem));
                    hpar(3).f = Yk(:,midx(1:numMem));
                    hpar(4).f = invM(1:numMem,1:numMem);
                else
                    hpar(1).f = gammak;
                    hpar(2).f = Psi(:,midx(1:numMem));
                    hpar(3).f = invM(1:numMem,1:numMem);
                end                
                [sk,~,infoLS,~] = lstrs(Hfunc,gk,Deltak,epsilon,EigFunc,...
                    lopts,hpar,eigensolverpar);
            catch ME %#ok<NASGU>
                iExit = 1;
                infoLS = -4;
            end
            if iExit == 1; info = 6; end;
            if infoLS== -4; info = 6; sk = zeros(n,1); end
        case 5
            if whichInit~=1&&whichInit~=2
                [sk,~,~,~] = st_truncated_CG_qn(gk,[],...
                    invM(1:numMem,1:numMem),gammak,Deltak,numMem,tolCG,...
                    maxitCG,showCG,Sk(:,midx(1:numMem)),Yk(:,midx(1:numMem)));
            else
                [sk,~,~,~] = st_truncated_CG_qn(gk,Psi(:,midx(1:numMem)),...
                    invM(1:numMem,1:numMem),gammak,Deltak,numMem,tolCG,maxitCG,showCG);
            end    
        otherwise
            [sk,iExit] = sc_sr1_infty( gk, Psi(:,midx(1:numMem)),...
                invM(1:numMem,1:numMem), gammak, Deltak, 1, 0);
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
    
    if whichInit == 1 || whichInit == 2    
        b(1:numMem,1) = Psi(:,midx(1:numMem))'*sk;
    else
        b1(1:numMem,1) = Yk(:,midx(1:numMem))'*sk;
        b2(1:numMem,1) = Sk(:,midx(1:numMem))'*sk;
        b(1:numMem,1) = b1(1:numMem,1) - gammak.*b2(1:numMem,1);
    end
    
    ss = sk'*sk;
    yy = yk'*yk;
    sy = sk'*yk;
    sBs = gammak*ss + b(1:numMem,1)'*(invM(1:numMem,1:numMem)\b(1:numMem,1));    
    gs = gk'*sk;
    
    switch whichInit
        case 3
            if sy > 0; gammak = yy/sy; end
        case 4    
            gammak = max(max(yy/sy,gams(1:min(numq,q),1)));      
            %gammak = min(gammak,1e13);
            mq = mod(numq,q);
            %if mq==0; gams(1,1) = gammak; else gams(mq+1,1) = gammak; end;        
            if mq==0; gams(1,1) = yy/sy; else gams(mq+1,1) = yy/sy; end;
            numq = numq+1;            
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
    
    psik = yk-gammak.*sk;
    %spsi = sk'*psik;
    spsi = sy - gammak.*ss;
       
    %psipsi = psik'*psik;
    psipsi = yy - 2*gammak*sy + gammak*gammak*ss;
    npsik = sqrt(psipsi);
    %npsik = norm(psik);
    
    % Update limited memory quasi-Newton vectors/matrices
    if(numMem < m) && (SR1tol*max(nsk*npsik,SR1tolAdj) < abs(spsi))
          %(SR1tol*nsk*npsik < abs(spsi))
        
        numMem = numMem + 1;
        
        % Constant gammak updating
        if whichInit == 1 || whichInit == 2
            Psi(:,midx(numMem)) = psik;
            b(numMem,1)         = spsi;
            invM(1:numMem,numMem) = b(1:numMem,1);
            invM(numMem,1:numMem) = invM(1:numMem,numMem);
            if storePsiPsi == 1
                numMemM1 = numMem-1;
                PsiPsi(1:numMemM1,numMem) = Psi(:,midx(1:numMemM1))'*psik;
                PsiPsi(numMem,numMem) = psipsi; %npsik*npsik;
                PsiPsi(numMem,1:numMemM1) = PsiPsi(1:numMemM1,numMem);
            end
        else % gammak not constant
            numMemM1 = numMem-1;
            Sk(:,midx(numMem)) = sk;
            Yk(:,midx(numMem)) = yk;
            Dk(numMem) = sy;
            Lk(numMem,1:numMemM1) = b1(1:numMemM1,1);
            Tk(1:numMem,numMem) = Sk(:,midx(1:numMem))'*yk;
            SkSk(1:numMemM1,numMem) = b2(1:numMemM1,1);
            SkSk(numMem,numMem) = ss;
            SkSk(numMem,1:numMemM1) = b2(1:numMemM1,1);
            YkYk(1:numMem,numMem) = Yk(:,midx(1:numMem))'*yk;
            YkYk(numMem,1:numMemM1) = YkYk(1:numMemM1,numMem);
            invM(1:numMem,1:numMem) = diag(Dk(1:numMem)) + Lk(1:numMem,1:numMem) +...
                Lk(1:numMem,1:numMem)' - gammak.*SkSk(1:numMem,1:numMem);
            % Storing the product Psi'*Psi without explicitly recomputing
            % it
            if storePsiPsi == 1
                PsiPsi(1:numMem,1:numMem) = YkYk(1:numMem,1:numMem)-...
                    gammak.*(Tk(1:numMem,1:numMem) + ...
                    Tk(1:numMem,1:numMem)' + ...
                    Lk(1:numMem,1:numMem) + ...
                    Lk(1:numMem,1:numMem)') + gammak.*gammak*SkSk(1:numMem,1:numMem);
            end
        end
        
    elseif(numMem == m) && (SR1tol*max(nsk*npsik,SR1tolAdj) < abs(spsi)) % Limited memory 
        %(SR1tol*nsk*npsik < abs(spsi)) % Limited memory 
        
        numMemM1 = numMem-1;
        midxs = midx(1);        
        
        if whichInit == 1 || whichInit == 2
            Psi(:,midxs) = psik;
        
            midx(1:(end-1)) = midx(2:end);        
            midx(end) = midxs;
            
            b(1:numMemM1,1) = b(2:numMem,1);
            b(numMem,1)     = spsi;
        
            invM(1:numMemM1,1:numMemM1) = invM(2:numMem,2:numMem);
            invM(1:numMem,numMem) = b(1:numMem,1);        
            invM(numMem,1:numMem) = invM(1:numMem,numMem);
        
            if storePsiPsi == 1
                            
                PsiPsi(1:numMemM1,1:numMemM1) = PsiPsi(2:numMem,2:numMem);
                PsiPsi(1:numMemM1,numMem) = Psi(:,midx(1:numMemM1))'*psik;
                PsiPsi(numMem,numMem) = psipsi; %npsik*npsik;
                PsiPsi(numMem,1:numMemM1) = PsiPsi(1:numMemM1,numMem);
            
            end
        else
            Sk(:,midxs) = sk;
            Yk(:,midxs) = yk;
            midx(1:(end-1)) = midx(2:end);        
            midx(end) = midxs;
            % Copy previous (small) arrays
            Dk(1:numMemM1) = Dk(2:numMem);
            Lk(1:numMemM1,1:numMemM1) = Lk(2:numMem,2:numMem);
            Tk(1:numMemM1,1:numMemM1) = Tk(2:numMem,2:numMem);
            SkSk(1:numMemM1,1:numMemM1) = SkSk(2:numMem,2:numMem);
            YkYk(1:numMemM1,1:numMemM1) = YkYk(2:numMem,2:numMem);
            
            % Updates
            % First temporary "b" buffers
            b1(1:numMemM1) = b1(2:numMem);
            b1(numMem) = sy;
            b2(1:numMemM1) = b2(2:numMem);
            b2(numMem) = ss;
            
            Dk(numMem) = sy;
            Lk(numMem,1:numMemM1) = b1(1:numMemM1,1);
            Tk(1:numMem,numMem) = Sk(:,midx(1:numMem))'*yk;
            SkSk(1:numMemM1,numMem) = b2(1:numMemM1,1);
            SkSk(numMem,numMem) = ss;
            SkSk(numMem,1:numMemM1) = b2(1:numMemM1,1);
            %YkYk(1:numMem,numMem) = Yk(:,midx(1:numMem))'*yk;
            YkYk(1:numMemM1,numMem) = Yk(:,midx(1:numMemM1))'*yk;
            YkYk(numMem,numMem) = yy;
            YkYk(numMem,1:numMemM1) = YkYk(1:numMemM1,numMem);
            
            invM(1:numMem,1:numMem) = diag(Dk(1:numMem)) + Lk(1:numMem,1:numMem) +...
                Lk(1:numMem,1:numMem)' - gammak.*SkSk(1:numMem,1:numMem);
            
            if storePsiPsi == 1                
                PsiPsi(1:numMem,1:numMem) = YkYk(1:numMem,1:numMem)-...
                    gammak.*(Tk(1:numMem,1:numMem) + ...
                    Tk(1:numMem,1:numMem)' + ...
                    Lk(1:numMem,1:numMem) + ...
                    Lk(1:numMem,1:numMem)') + gammak.*gammak*SkSk(1:numMem,1:numMem);
            end
        end
    end
    

    %----------------------- Error computations --------------------------%
%      In = eye(n);        
%      if whichInit ~= 1 && whichInit ~= 2
%          Psi = Yk(:,1:numMem) - gammak.*Sk(:,1:numMem);
%          
%          errSS = norm(SkSk(1:numMem,1:numMem)-...
%              Sk(:,midx(1:numMem))'*Sk(:,midx(1:numMem)),'fro');         
%          errYY = norm(YkYk(1:numMem,1:numMem)-...
%              Yk(:,midx(1:numMem))'*Yk(:,midx(1:numMem)),'fro');         
%          SY = Sk(:,midx(1:numMem))'*Yk(:,midx(1:numMem));         
%          errT = norm(Tk(1:numMem,1:numMem)-triu(SY(1:numMem,1:numMem)),'fro');
%          errL = norm(Lk(1:numMem,1:numMem)-tril(SY(1:numMem,1:numMem),-1),'fro');
%          errD = norm(Dk(1:numMem)-diag(SY(1:numMem,1:numMem)),'fro');
%      end
%      Bk = gammak*In + Psi(:,midx(1:numMem))*...
%         (invM(1:numMem,1:numMem)\Psi(:,midx(1:numMem))');
%      err3 = norm(Bk*sk-yk); % Secant condition
%      if storePsiPsi == 1
%          err4 = norm(PsiPsi(1:numMem,1:numMem)-...
%              Psi(:,midx(1:numMem))'*Psi(:,midx(1:numMem)),'fro');
%      end
    %-------------------- End error computations -------------------------%
    
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

