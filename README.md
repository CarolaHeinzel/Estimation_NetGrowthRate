This code can be used to estimate the net growth rate. It uses parts of the R package [cloneRate](https://github.com/bdj34/cloneRate?tab=readme-ov-file)

### Overview of the Files

Estimation_Durrett.R: Estimation of the net growth rate with the formula of Durrett. <br>
calc_prefactor.R: Determination of the prefactors for \hat r_{JS, MSE}. <br>
simulation_CPP.R: Simulation of the CPP, for which we need some new outputs compared to cloneRate. However, the code is basically the same as in cloneRate. <br>
main.R: Comparison of the different estimation methods, including Phylofit and the new ones. <br> 
Create_Plot.R: Plot the results. <br>

### Usage

To use the code to calculate the net growth rate with the new methods, just apply the function calc_esimator() in new_estimation.R.
