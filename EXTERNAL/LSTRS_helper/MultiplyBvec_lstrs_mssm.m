%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [ Bx ] = MultiplyBvec_lstrs_mssm( x, Hpar )
%MultiplyBVec multiplies a vector with the limited memory matrix B
%using the compact representation of B:
% Bx = gammaIx + Psi M * Psi'x
%
% Output : 
% Bx := nxn vector.
%
% Dervied from previous L-SR1 methods. 
% Please see the following technical report:
%
% "OBS: MATLAB Solver for L-SR1 Trust-Region Subproblems"
% by Johannes Burst, Jennifer Erway, and Roummel Marcia
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 02/10/21, J.B., Minor modification to use "\" in computing the product
% 03/02/21, J.B., Instead of Psi, the components Yk,Sk were passed 
%                   (i.e., the length of number of inputs differs)
% 06/10/22, J.B., Update for the MSS matrix


Bx = (Hpar(1).f).*x + (Hpar(2).f)*((Hpar(3).f)*((Hpar(2).f)'*x));

