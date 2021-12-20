##  Basic steps for using scripts and models in the 'case_study' folder of the 'gt_dist' repository accompanying the Nuse et al.'s 'Alternatives to the g(0) = 1 assumption in distance sampling of populations where an individual covariate influences detection'.
##  
##  
##  Dec 2021
#########################################################################

System requirements:
1) program R  (https://cran.r-project.org)
2) JAGS  (https://sourceforge.net/projects/mcmc-jags)
3) R pacakges: rjags, parallel, xtable

System note:
Parallel processing of MCMC chains is carried out with function mclapply() in R package {parallel}.  This function uses process forking, which does not work on Microsoft Windows.  See notes in the file 'case_study/functions/ich_funs.r' if you're using Windows.  


Steps: interactive version.
----------------------------
1)  Start an interactive R session; set working directory to the top directory of the repository (i.e., the one containing sub-directories 'case_study' and 'simulations').
2)  Open file 'case_study/scripts/ich_script.r'
3)  Load data and functions at the top of the script; run particular models as desired.
4)  Open file 'case_study/scripts/ich_collect_results_nov2019.r' and gather results and summarize; produce plots.


Steps: script version.
-----------------------
1)  Open a shell terminal, set working directory to the top directory of the repository (i.e., the one containing sub-directories 'case_study' and 'simulations').
2)  run, from the command line:
Rscript ./case_study/scripts/ich_script.r
(This will take a while.  Output files will be deposited in 'case_study/output'.)
3)  run:
Rscript ./case_study/scripts/ich_script.r
(Figure files will be deposited in 'case_study/figures'.)


