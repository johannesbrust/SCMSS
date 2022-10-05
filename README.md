# LMSS_SC
###########################################################################
                                                                       
 LMSS_SC: Shape-changing limited-memory Multipoint Symmetric Secant<br/> 
 trust-region algorithms for large-scale optimization problems                                     
                                                                       
                            minimize f(x)                                   
 
 The methods seek a point "norm(grad(x)) = 0", where grad(x) is<br/>
 the gradient of the objective f(x).

 J.J Brust, J.B. Erway, R.F. Marcia (2022)                                                                        
###########################################################################

## Setup
Files are organized in 6 folders:

	ALGS/ 		(main algorithms)
	AUXILIARY/ 	(support routines)
	DATA/ 		(data for figures and tables)
	EXPERIMENTS/ 	(examples and experiments)
	EXTERNAL/ 	(third party software)
	FIGS/		(figures)

To rerun any experiment from the folder "EX_COMP_LMSSM_LSR1_m_q5" the 
CUTEst (https://github.com/ralna/CUTEst/wiki) test set has to be installed. 
If CUTEst is installed, then the relative path to the CUTEst installation
is to be updated in 

	AUXILIARY/CUTEst_init.m

Without a CUTEst installation, the figures can still be plotted 
from within the "AUXILIARY/" folder by calling 
fig*.m, *={1,2,3,4,5,6,7,8}

Detailed descriptions of the signatures of the algorithms in "ALGS/"
are found in the beginning comments of the respective files

The codes have been tested on:

Matlab R2016a, macOS Catalina, {CUTEst n/a, July 2017}

## Examples
To try-out three examples navigate inside folder "EXPERIMENTS/".
From within the folder you can run the examples by

	>> example_1

	>> example_2

	>> example_3

