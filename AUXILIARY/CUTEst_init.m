function CUTEst_init()

if ispc
  disp('Running on Windows');
  disp('No installation available!');
end

if isunix
  disp('Running on Unix');
  
  % Please adapt the next line for your CUTEr/(st) installation folder  
  %addpath('~/Documents/CUTEST/cutest/src/matlab');
  %addpath('~/Documents/CUTEST_201/cutest/src/matlab');
  addpath('/Users/jjbrust/Documents/CUTEST_201/cutest/src/matlab');
end