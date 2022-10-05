%% Function from CUTEST EXPERIMENT to retrieve outputs of running an algorithm.
function [ex,numf,numg,numit,tcpu,tract,numrst,numskip,var1,var2,var3,var4]=runAlgorithm(alg,...
        obj,grad,x0,params,numRuns)
tcpu=zeros(numRuns,1);

for i=1:numRuns
    %[~,~,outinfo]=alg(obj,cons,x0,params);
    [~,~,~,outinfo]=alg(x0,obj,grad,params);
        %[~,~,outinfo]=alg(@cuter_fun,x0,params);
    tcpu(i)=outinfo.ctime;
end

numf=outinfo.numf;
numg=outinfo.numg;
ex=outinfo.ex;
numit=outinfo.numiter;   
tract=outinfo.numTRit;
if nargout==7
    numrst=outinfo.numAccept;
elseif nargout==8
    numrst=outinfo.numAccept;
    numskip=outinfo.numTRInc;
elseif nargout == 9
    numrst=outinfo.numAccept;
    numskip=outinfo.numTRInc;
    var1    = outinfo;
elseif nargout == 10
    numrst=outinfo.numAccept;
    numskip=outinfo.numTRInc;
    var1    = outinfo;
    var2    = outinfo.varout3;
elseif nargout == 11
    numrst=outinfo.numAccept;
    numskip=outinfo.numTRInc;
    var1    = outinfo;
    var2    = outinfo.varout3;
    var3    = outinfo.varout4;
elseif nargout == 12 
    numrst=outinfo.numAccept;
    numskip=outinfo.numTRInc;
    var1    = outinfo;
    var2    = outinfo.varout3;
    var3    = outinfo.varout4;
    var4    = outinfo.varout5;
end
