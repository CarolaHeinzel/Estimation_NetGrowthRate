This code can be used to estimate the net growth rate as described in the paper [""]() by Carola Heinzel and Jason Schweinsberg. It uses parts of the R package [cloneRate](https://github.com/bdj34/cloneRate?tab=readme-ov-file)


### Overview of the Files and Folders

The following files might be important to apply our method to new data: <br>

new_Estimator.R: Implementation of \hat r_{new}. <br> 
Calc_quantile.R: Calculation of the quantiles. <br>
CPP_infinity.R: Determination of the prefactors for \hat r_{MSE} and simulates CPPs for T = \infty. <br>

The following folders contain code to reproduce the results from the paper. <br>
Data: Simulated values for the internal branch lengths of a coalescent point process (called S(n) in the paper). <br>
Comparison_methods: Comparison of the different estimation methods, including Phylofit and the new ones. <br> 
Plots: Contains the code to create the figures in our paper. <br>

### Usage

To use the code to calculate the net growth rate with the new methods, just apply the function calc_esimator() in new_estimation.R.
