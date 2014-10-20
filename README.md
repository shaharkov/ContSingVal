Controlling Singular Values with Semidefinite Programming
====

A Matlab implementation of the paper ["Controlling Singular Values with Semidefinite Programming"](http://www.wisdom.weizmann.ac.il/~shaharko/ContSingVal.html).

----
This package includes three examples:
- `example_optimizeSingleMatrix.m` demonstrates simple single matrix optimization problems with constraints on singular values. Algorithm 1 of the paper is applied to a few example optimization problems, aiming to illustrate the implementation and usage of the theory presented in the paper (section 4 and 5). For example, solving problems of the form
```
minimize        f(A) 
subject to      s_min(A)>=k
                s_max(A)<=K
				det(A)>=0
```
over the matrix variable `A`, where `s_min(A)` and `s_max(A)` are the minimal and maximal singular values of the matrix `A`.

- `example_BarDeformation.m` demonstrates volumetric mesh optimization problems with constraints on singular values. Algorithm 2 of the paper is applied to some example optimization problems, generating some of the example deformations of Figure 2 in the paper. See section 6.1 in the paper for additional details. Either run the entire code for generating all examples (might take a while), or the relevant code block for your desired example.

- `example_ExtremalQuasiConformal.m` implements an algorithm for computing extremal quasiconformal mappings of volumetric meshes (i.e., minimizing maximal conformal distortion). The code reproduces the example presented in Figure 1. See section 6.1 for additional details.


**Getting started:** Before running the code, make sure you have YALMIP and MOSEK installed. Update the paths in `initialize.m` accordingly. The code was tested with Matlab (2013a), YALMIP (20140915) and MOSEK (7.0.0.92).

**Disclaimer:**
The code is provided as-is for academic use only and without any guarantees. Please contact the authors to report any bugs.
Written by [Shahar Kovalsky](http://www.wisdom.weizmann.ac.il/~shaharko/) and [Noam Aigerman](http://www.wisdom.weizmann.ac.il/~noamaig/)