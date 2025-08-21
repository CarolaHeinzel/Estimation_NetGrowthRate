This code can be used to estimate the net growth rate as described in the paper ["Estimating the growth rate of a birth and death process using data from a small sample"]() by Carola Heinzel and Jason Schweinsberg. It uses parts of the R package [cloneRate](https://github.com/bdj34/cloneRate?tab=readme-ov-file).


### Overview of the Files and Folders

The following files might be important to apply our method to new data: <br>

New_Estimator.R: Implementation of \hat r_{MSE}, \hat r_{Bias} and \hat r_{Inv}. This file includes an example how to use it to new data. <br> 
CPP_infinity.R: Determination of the prefactors for \hat r_{MSE} and simulates CPPs for T = \infty. This can be used to calculate \hat r_{MSE} for other samples sizes n than in the paper. Basically, this script can also be used to calculate c_{Bias}. The calculation of c_{Inv}(n) is analytical. <br>

The following folders contain code to reproduce the results from the paper. <br>
Results_CPP: Simulated values for the internal branch lengths of a coalescent point process (called S(n) in the paper). <br>
Comparison_methods: Comparison of the different estimation methods, including Phylofit and the new ones. <br> 
