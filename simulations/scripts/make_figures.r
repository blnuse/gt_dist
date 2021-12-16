##  Make figures for Nuse et al. manuscript 'Alternatives to the g(0) = 1 assumption...' 
##
##
##  16 dec 2021
##########################################3

##  Set working directory to top level of the 'gt_dist' git repository (where folders 'simulations' and 'case_study' live.)

# Fig. 3 -- simulated survey includes stochastic processes, so it will never look the same twice.
# ------------------------------------
  source('simulations/scripts/new_sim_09jul2019')
  source('simulations/scripts/distance_simulation_plots_jul2019.txt')
  
  s <- burrSim(method="stable-age", S.juv = 0.7, S.subad = 0.85, S.ad = 0.95, F = 4, lambda = 0.94, mult = 20, sigma = 30, intrvl = 100, xi.brk = 300, xi.50 = 0.1)
  map.burrSim(s)
  plot.p.burrSim(s)

  
  
