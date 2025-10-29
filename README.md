This code can be used to estimate the net growth rate of a birth and death process as described in the paper ["Estimating the growth rate of a birth and death process using data from a small sample"](https://arxiv.org/abs/2508.16110) by Carola Heinzel and Jason Schweinsberg. It uses parts of the R package [cloneRate](https://github.com/bdj34/cloneRate?tab=readme-ov-file).


### Overview of the Files and Folders

The following files can be used to apply our method to new data: <br>

* New_Estimator.R: Implementation of \hat r_{MSE}, \hat r_{Bias} and \hat r_{Inv}. This file includes an example how to use it to new data. <br> 
* CPP_infinity.R: Determination of the prefactors for \hat r_{MSE} and simulates CPPs for T = \infty. This can be used to calculate \hat r_{MSE} for other samples sizes n than in the paper. Basically, this script can also be used to calculate c_{Bias}. The calculation of c_{Inv}(n) is analytical. <br>

The following folder contains code to reproduce the results from the paper. <br>
* Comparison_methods: Comparison of the different estimation methods, including Phylofit and the new ones. This also includes to apply \hat r_{MSE} and \hat r_{MCMC} to real data. <br> 


### Packages and Dependencies

We tested the code with:

- **R version:** 4.5.1  
- **Platform:** x86_64-w64-mingw32  
- **R packages:**
  - cloneRate 0.2.3
  - ape 5.8.1
  - gtools 3.9.5
  - pracma 2.4.4
  - Rmpfr 1.1.0
  - reshape2 1.4.4


### Funding Acknowledgement

Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) – Project-ID 499552394 – SFB 1597.
