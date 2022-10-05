# SCMSS
                                                                       
 SCMSS: Shape-changing limited-memory Multipoint Symmetric Secant 
 trust-region algorithms for large-scale optimization problems                                     
                                                                       
                            minimize f(x)                                   
 
 The methods seek a point "norm(grad(x)) = 0", where grad(x) is
 the gradient of the objective f(x).

 J.J Brust, J.B. Erway, R.F. Marcia (2022), [[Article](https://arxiv.org/abs/2209.12057 "Technical Report")]                                                                        


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

The outcomes of running example_1 look like:
```
>> example_1

Example 1 ##########################################
Rosenbrock objective: f(x)                         
                      n = 100                       
                     f0 = 8.1089e+05                    
               norm(g0) = 1.0807e+05                    

Algorithm: LMSS_SC                                  
####################################################

----------- Running Algorithm ----------- 
 L-MSS Trust-Region                      
Iter 	 fk      	 norm(gk)  	 TR     	 numf  
0 	 8.1089e+05 	 1.0806e+05 	 1.0000e+05 	 1  
1 	 6.4527e+01 	 2.5728e+01 	 6.3640e+01 	 7  
2 	 6.4266e+01 	 2.5418e+01 	 6.3640e+01 	 8  
3 	 6.4011e+01 	 2.5115e+01 	 6.3640e+01 	 9  
4 	 6.3760e+01 	 2.4818e+01 	 6.3640e+01 	 10  
5 	 6.3514e+01 	 2.4527e+01 	 6.3640e+01 	 11  
6 	 6.3272e+01 	 2.4242e+01 	 6.3640e+01 	 12  
7 	 6.3035e+01 	 2.3962e+01 	 6.3640e+01 	 13  
8 	 4.6265e+01 	 4.8674e+00 	 6.3640e+01 	 14  
9 	 4.1142e+01 	 1.9100e+00 	 6.3640e+01 	 15  
10 	 3.0268e+01 	 1.5262e+00 	 6.3640e+01 	 16  
11 	 3.0268e+01 	 1.5262e+00 	 3.1820e+01 	 17  
12 	 3.0268e+01 	 1.5262e+00 	 1.5910e+01 	 18  
13 	 3.0268e+01 	 1.5262e+00 	 7.9550e+00 	 19  
14 	 3.0268e+01 	 1.5262e+00 	 3.9775e+00 	 20  
15 	 2.0653e+01 	 1.6094e+00 	 7.9550e+00 	 21  
16 	 1.5106e+01 	 1.9218e+00 	 7.9550e+00 	 22  
17 	 1.5106e+01 	 1.9218e+00 	 3.9775e+00 	 23  
18 	 6.9360e+00 	 2.4369e+00 	 3.9775e+00 	 24  
19 	 6.9360e+00 	 2.4369e+00 	 1.9888e+00 	 25  
20 	 6.1816e+00 	 8.8037e+00 	 9.9438e-01 	 26  
21 	 2.8697e+00 	 2.3486e+00 	 9.9438e-01 	 27  
22 	 2.2602e+00 	 1.2387e+00 	 4.9719e-01 	 28  
23 	 5.2984e-01 	 1.3380e+00 	 4.9719e-01 	 29  
24 	 2.5109e-01 	 1.3279e+00 	 4.9719e-01 	 30  
25 	 2.5109e-01 	 1.3279e+00 	 2.4859e-01 	 31  
26 	 1.1164e-01 	 7.2574e-01 	 2.4859e-01 	 32  
27 	 3.8678e-02 	 1.1388e-01 	 4.9719e-01 	 33  
28 	 1.8057e-03 	 4.3420e-02 	 9.9438e-01 	 34  
29 	 4.5221e-05 	 8.9849e-03 	 9.9438e-01 	 35  
30 	 1.2118e-06 	 2.8453e-03 	 9.9438e-01 	 36  
31 	 2.4660e-07 	 2.6741e-04 	 9.9438e-01 	 37  
32 	 3.3913e-09 	 2.3060e-04 	 9.9438e-01 	 38  
33 	 1.0920e-11 	 1.0197e-05 	 9.9438e-01 	 39  
34 	 1.0324e-17 	 3.2252e-09 	 9.9438e-01 	 40
```

## Cite
You can cite this work as (bibtex)

```
@TechReport{scmss22,
  author      = {Johannes J. Brust, Jennifer B. Erway and Roummel F. Marcia},
  title       = {Shape-Changing Trust-Region Methods Using Multipoint Symmetric Secant Matrices},
  institution = {Mathematics Department, University of California, San Diego, CA},
  type        = {Technical Report},
  year        = {2022},
  url         = {https://doi.org/10.48550/arXiv.2209.12057}
}
```
