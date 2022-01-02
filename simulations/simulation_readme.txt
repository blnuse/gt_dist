##  Guide to scripts in the 'simulations' folder of the 'gt_dist' repository accompanying the Nuse et al.'s 'Alternatives to the g(0) = 1 assumption in distance sampling of populations where an individual covariate influences detection'.
##  
##  
##  Dec 2021
#########################################################################

System requirements:
1) program R  (https://cran.r-project.org)
2) JAGS  (https://sourceforge.net/projects/mcmc-jags)
3) R pacakges: rjags, parallel, xtable, shape, viridis



None of the script files provided in folder 'simulations/scripts' is intended to be run with Rscript; they should all be run interactively, some have particulr instructions and variables to set before proceeding.

The simulated datasets and model fits summarized in the manuscript were produced using the script 'simulations/scripts/distance_simulation_automated_oct2019.r'.

Plots relying on those model fits were made with 'simulations/scripts/distance_simulation_collect_results_jul2019.r'.

A newer version of the burrow survey simulator is also included, and was used to make the example survey summary figures in the manuscript.  See 'simulations/scripts/guide_new_sim.txt' and 'simulations/scripts/make_figures.r'.


