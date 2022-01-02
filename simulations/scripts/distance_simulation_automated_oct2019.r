#  Simulating DISTANCE analysis data.
#  Parallelized, and automated (to a degree)! 
#
#
#  13 May 2015; Oct 2016
##############################################

####  SET THESE:
#  iter.start, the replicate number on which to start
#  iter.stop, the replicate number on which to stop
#  xi.switch : TRUE if xi is size-dependent in the simulated data; FALSE if sigma is size dependent
###  e.g., iter.start <- 1; iter.stop <- 20; xi.switch <- TRUE; doll <- FALSE

####  THEN TO CALL THIS SCRIPT:
# source('simulations/scripts/distance_simulation_automated_oct2019.r')



     library(parallel)
     library(rjags); load.module("glm"); load.module("lecuyer"); load.module("dic")
     library(random)
     
####  Set up the cores for the chains:     
     registerDoParallel(cl = 3, cores = 3)
     
     
####  Stuff that applies to all iterations of the following loops.

      
     
####  Applies to all pop's:
 edge <- 1000

####  Functions Bay:  ############################
           trans.y <- function (x, trans.num = 1, theta, interval, y.int) {
            theta.rad <- theta*pi/180
            tan(theta.rad)*x + (y.int + trans.num*interval/cos(theta.rad))
          }

          trans.x <- function (y, trans.num = 1, theta, interval, y.int) {
            theta.rad <- theta*pi/180
            (y - (y.int + trans.num*interval/cos(theta.rad)))/tan(theta.rad)
            }
##################################################



######### OUTER LOOP WILL START HERE #######################

# desired pop sizes:

pop.sz <- c(800, 400, 200, 100, 50)

int.vals <- c(100, 200)

max.diam <- 400

for (h in int.vals) {

 
#----------------------
# Lay out transects:

 theta.choice <- 50
 theta.ch.rad <- theta.choice*pi/180
 int.choice <- h

if (xi.switch) {
  if (doll) {
    setwd(paste("simulations/output/sim_modfits/nested_doll/interval", int.choice, sep = "_"))
    } else {setwd(paste("simulations/output/sim_modfits/xi_size/interval", int.choice, sep = "_"))}
  } else {
    setwd(paste("simulations/output/sim_modfits/sigma_size/interval", int.choice, sep = "_"))
    }

     
######### MIDDLE LOOP WILL START HERE #######################

     for (k in pop.sz) {

          fnm.base <- paste("dist_result", k, gsub("-", "_", Sys.Date(), fixed = T), sep = "_")


     ###  Pop setup:
       n.burr <- k
      
      
       ###  Find y-intercept for lowest line:
      opt <- optimize (
        f = function(y) {abs(trans.y(x = 1000, theta = theta.choice, interval = int.choice, y.int = y))},
        interval = c(-10*edge, 10*edge)
      )

      y.int.min <- opt$minimum + runif(1, min = 0.2, max = 0.8)*sin(theta.ch.rad)*int.choice
      
      
      ###  Find number of necessary transect lines to cover extent:
      n.lines <- floor((edge -  y.int.min)/(int.choice/cos(theta.ch.rad)))

      
     ###  OK, here I'll want to select whether or not:
               ## -- sigma is size-dependent
               ## -- xi is size-dependent
               ## -- both


      ###  Size-dependent sigma.  Hold relevant scalar values in a list:
      vl <- list()
      
      vl$xi.50 <- 0.3
      vl$xi.break <- 200
      vl$xi.beta <- (1 - vl$xi.50)/(vl$xi.break - 50)
      # vl$xi.beta.sq <- NA #(1 - vl$xi.int)/vl$xi.break^2
      
      ###  Fixed sigma, for now
      vl$sigma.int <- 1
      vl$sigma.beta <- 0.005     ## log(20) - 1
      
      

      
      
      
      #########################  BEGIN INNER LOOP #########################

     for (j in iter.start:iter.stop) {
           

            bw <- runif(n.burr, 50, max.diam)
           ###  beta:     
            # bw <- 5+35*rbeta(300, shape1 = 1.8, shape2 = 1.3)
           ###  Try a truncated normal:
           # require(truncnorm)
           # bw <- rtruncnorm(300, a = 5, b = 40, mean = 33, sd = 15)

           loc.x <- runif(n.burr, 0, edge)
           loc.y <- runif(n.burr, 0, edge)
           

           ##  find distance of each burrow to nearest transect, and which that is!
           ##  d = abs(y)*cos(theta.ch.rad)
           diff.mat <- matrix(NA, ncol = n.lines, nrow = n.burr)
           for (i in 1:n.lines) { 
             diff.mat[,i] <- abs(loc.y - trans.y(x = loc.x, trans.num = i, theta = theta.choice, interval = int.choice, y.int = y.int.min)) * cos(theta.ch.rad)
           }
           
           near.trans <- apply(diff.mat, MAR = 1, which.min)
           min.dist <- apply(diff.mat, MAR = 1, min)

            ###  Put all together:
           df <- data.frame(loc.x = round(loc.x, 1), loc.y = round(loc.y, 1), diam = round(bw, 1), trans = near.trans, dist = round(min.dist, 1))
           
                                 
           ### calculate sigma's
           # 
           
           ### calculate xi's
           df$xi <- vl$xi.50 + vl$xi.beta*(df$diam-50)
           df$xi <- ifelse(df$xi > 1, 1, df$xi)
           
           ### calculate p 
           if (xi.switch) {
             vl$sigma <- 20
             df$p <- df$xi*exp(-(df$dist^2)/(2*vl$sigma^2))
             } else {
             df$sigma <- exp( vl$sigma.int + vl$sigma.beta*df$diam )
             df$p <- exp(-(df$dist^2)/(2*df$sigma^2))
             }

           ###  Now do sampling!
           df$found <- rbinom(n = df$p, size = 1, prob = df$p)

           obs <- df[df$found == 1, c("found",  "trans", "dist", "diam")]

           aug <- 4*dim(obs)[1]
           
           ### append augmented burrows:
           app <- as.data.frame(matrix(rep(c(0, rep(NA, dim(obs)[2]-1)), times = aug), ncol = dim(obs)[2], byrow = T))
           names(app) <- names(obs)
#           app$diam <- rep_len(obs$diam, aug)
#           app$diam <- round(50+rbeta(aug, 0.2, 0.9)*350,1)
       
           full <- rbind(obs, app)

           M <- dim(full)[1]
           
           ###  Strip width, for observed burrows:
           Bx <- max(full$dist, na.rm = T)
           

           ###  Find the area actually sampled, given Bx
           ###  This is approximately the total length minus, for each transect, Bx^2*(tan(theta.ch.rad) + 1/tan(theta.ch.rad))
           

           len <- NA
           hits.side <- NA
           overest <- NA
           
           for (i in 1:n.lines) { 
             beg.x <- min ( max( c(0, trans.x(y = 0, theta = theta.choice, interval = int.choice, y.int =  y.int.min, trans.num = i)) ), 1000)
             beg.y <- min( max( c(0, trans.y(x = 0, theta = theta.choice, interval = int.choice, y.int =  y.int.min, trans.num = i)) ), 1000)
             end.x <- min ( max( c(0, trans.x(y = 1000, theta = theta.choice, interval = int.choice, y.int =  y.int.min, trans.num = i)) ), 1000)
             end.y <- min( max( c(0, trans.y(x = 1000, theta = theta.choice, interval = int.choice, y.int =  y.int.min, trans.num = i)) ), 1000)
           
             len[i] <- sqrt((beg.x - end.x)^2 + (beg.y - end.y)^2)
             if (beg.x == 0) {
                 hits.side[i] <- ifelse(end.x != 1000, "both", "side-side")
                 } else {
                 hits.side[i] <- ifelse(end.x != 1000, "top-bottom", "both")
                 }

             #####  This NEEDS TO BE CHECKED!!! --> Incorporate hits.side !!!! ---> ok, done: it's still not perfect, missing some extra area when a transect is near a corner!
             if (hits.side[i] == "both") {
                 overest[i] <- 0.5*Bx^2*(tan(theta.ch.rad) + 1/tan(theta.ch.rad))
                 }
             if (hits.side[i] == "top-bottom") {
                 overest[i] <- Bx^2*(1/tan(theta.ch.rad))
                 }       
             if (hits.side[i] == "side-side") {
                 overest[i] <- Bx^2*(tan(theta.ch.rad))
                 }       

            }       
            
            
               area <- 2*Bx*sum(len) - sum(overest) 
           
           
           ###  Should be ready now to fit the model....

           
          ###  Number of burrows observed:
          N.obs <- dim(obs)[1]

          ###  Number of burrows in the true population that are within Bx:
          N.Bx <- dim(df[df$dist <= max(full$dist, na.rm = T),])[1]


           ########################### ########################### ###########################
            ########################### ########################### ###########################
               ########################### ########################### ###########################  Now do fitting with JAGS....
               

               
               j.dat <- with(full, list(Y = found,
                                           diam = diam, 
                                           trans = trans,
                                           x = dist, 
                                           M = M,
                                           Bx = Bx,
                                           A.samp = area,
                                           A.tot = edge^2,
                                           bin.floor = seq(50,max.diam,50),
                                           max.diam = max.diam))
               
               j.init <- list(w = full$found)

           
               
               ######   PARALLELIZE !! #################################################
               
               init.fun <- function() {
                  ### all the other params ###
                 list(.RNG.name = "lecuyer::RngStream",
                       .RNG.seed = sample(x = 1e+06, size = 1),          # as.numeric(randomNumbers(n = 1, min = 1, max = 1e+06,col=1)),
                       w = unlist(j.init$w)               
                       )
               }
               


           
          #     ----- Xi !! --------
          #####  For single-site xi model:   ######################################################
          time.xi <- system.time(
               samp.xi <- foreach(i=1:getDoParWorkers(), .combine = "c", .final = mcmc.list, .verbose = T) %dopar% {  ##.export = c("randomNumbers", "jags.model", "coda.samples", "load.module"), 
          #       load.module("lecuyer")
          #       load.module("dic")
                 
                 jags.mod <- jags.model(file = "simulations/jags/R_&_D_JAGS_broken_stick_single.txt",
                                                data = j.dat,
                                                inits=init.fun,
                                                n.chain=1, n.adapt=1000)
                 update(jags.mod, n.iter = 5000)                                 
                 result <- coda.samples(jags.mod, 
                                              variable.names = c("sigma", "xi.50", "xi.break", "N", "N.tot", "deviance"),
                                              n.iter = 5000)
                 return(result)
               }
               ) 
               s.xi <- summary(samp.xi)
               conv.xi <- gelman.diag(samp.xi, mult = F)
               
          #     ----- Sigma constant!! --------
          #####  For single-site sigma-only model:   ######################################################
          time.sigma <- system.time(
               samp.sigma <- foreach(i=1:getDoParWorkers(), .export = c("randomNumbers", "jags.model", "coda.samples", "load.module"), .combine = "c", .final = mcmc.list, .verbose = T) %dopar% {
                 load.module("lecuyer")
                 load.module("dic")
                 
                 jags.mod <- jags.model(file = "simulations/jags/R_&_D_JAGS_just_sigma.txt",
                                                data = j.dat,
                                                inits=init.fun,
                                                n.chain=1, n.adapt=1000)
                 update(jags.mod, n.iter = 5000)                                 
                 result <- coda.samples(jags.mod, 
                                              variable.names = c("sigma", "N", "N.tot", "deviance"),
                                              n.iter = 5000)
                 return(result)
               }
               ) 
               s.sigma <- summary(samp.sigma)
               conv.sigma <- gelman.diag(samp.sigma, mult = F)
               
          #     ----- Sigma size-dependent!! --------
          #####  For single-site sigma-only model:   ######################################################     
          time.sigsize <- system.time(
               samp.sigsize <- foreach(i=1:getDoParWorkers(), .combine = "c", .final = mcmc.list, .verbose = T) %dopar% {  ##.export = c("randomNumbers", "jags.model", "coda.samples", "load.module"), 
                 load.module("lecuyer")
                 load.module("dic")
                 
                 jags.mod <- jags.model(file = "simulations/jags/R_&_D_JAGS_sigsize.txt",
                                                data = j.dat,
                                                inits=init.fun,
                                                n.chain=1, n.adapt=1000)
                 update(jags.mod, n.iter = 5000)                                 
                 result <- coda.samples(jags.mod, 
                                              variable.names = c("sigma.int", "sigma.beta", "N", "N.tot", "deviance"),
                                              n.iter = 5000)
                 return(result)
               }
               ) 
               s.sigsize <- summary(samp.sigsize)
               conv.sigsize <- gelman.diag(samp.sigsize, mult = F)

          #     ----- Xi vs. Sigma Decision !! --------
          #####  For single-site model with switch to determine whether xi or sigma is size-dependent:   ############
          time.decide <- system.time(
               samp.decide <- foreach(i=1:getDoParWorkers(), .combine = "c", .final = mcmc.list, .verbose = T) %dopar% { ##.export = c("randomNumbers", "jags.model", "coda.samples", "load.module"), .combine = "c", .final = mcmc.list, .verbose = T) %dopar% {
                 load.module("lecuyer")
                 load.module("dic")

                 jags.mod <- jags.model(file = "simulations/jags/R_&_D_JAGS_xi_sigma_decision.txt",
                                                data = j.dat,
                                                inits=init.fun,
                                                n.chain=1, n.adapt=1000)
                 update(jags.mod, n.iter = 5000)                                 
                 result <- coda.samples(jags.mod, 
                                              variable.names = c("sig.size", "sigma.int", "sigma.beta", "sigma.const", "xi.50", "xi.break", "N", "N.tot", "deviance"),
                                              n.iter = 5000)
                 return(result)
               }
               ) 
               s.decide <- summary(samp.decide)
               conv.decide <- gelman.diag(samp.decide, mult = F)
               
               
               
          #     ----- Nested Doll !! --------
          #####  For Nested Doll single-site model:   ############
          time.RD <- system.time(
               samp.RD <- foreach(i=1:getDoParWorkers(), .combine = "c", .final = mcmc.list, .verbose = T) %dopar% { ##.export = c("randomNumbers", "jags.model", "coda.samples", "load.module"), 
                 load.module("lecuyer")
                 load.module("dic")
                 
                 jags.mod <- jags.model(file = "simulations/jags/R_&_D_JAGS_Nested_Doll.txt",
                                                data = j.dat,
                                                inits=init.fun,
                                                n.chain=1, n.adapt=1000)
                 update(jags.mod, n.iter = 5000)                                 
                 result <- coda.samples(jags.mod, 
                                              variable.names = c("sigma", "xi.50", "xi.break", "N", "N.tot", "deviance"),
                                              n.iter = 5000)
                 return(result)
               }
               ) 
               s.RD <- summary(samp.RD)
               conv.RD <- gelman.diag(samp.RD, mult = F)
               
               

          #######################################  Collect results, summaries, etc.



          ###  Save results:
               simul.list <- list(extent = list(value = edge^2/1000^2, units = "km2"),
                                      N.tot = list(value = n.burr, descr = "total number of burrows in the sample frame"),
                                      size.distr = list(value = "uniform", descr = "size distribution of burrows"),
                                      n.transects = list(value = n.lines, descr = "number of transects"),
                                      trans.len = list(value = sum(len)/1000, units = "km", descr = "total length of all transects"),
                                      Bx.area = list(value = area/1000000, units = "km2", descr = "total area within the strip"),
                                      sigma = list(value = vl$sigma, descr = "fixed value of sigma"),
                                      sigma.int = list(value = vl$sigma.int, descr = "intercept of sigma's linear predictor"),
                                      sigma.beta = list(value = vl$sigma.int, descr = "coefficient on diam, in sigma's linear predictor"),
                                      xi.int = list(value = vl$xi.int, descr = "intercept of xi's linear predictor"),
                                      xi.beta = list(value = vl$xi.beta, descr = "coefficient on diam, in xi's linear predictor"),
                                      Bx = list(value = Bx, descr = "strip width"),
                                      N.Bx = list(value = N.Bx, descr = "true count of burrows within the strip"),
                                      N.obs = list(value = N.obs, descr = "number of burrows actually observed"),
                                      theta = list(value = theta.choice, units = "degrees", descr = "angle of the transects"),
                                      interval = list(value = int.choice, units = "meters"))

               samp.list <- list(xi.size = samp.xi,
                                     sigma.size = samp.sigsize,
                                     sigma.const = samp.sigma,
                                     decide = samp.decide,
                                     nested.doll = samp.RD)
               
               summary.list <- list(xi.size = s.xi,
                                         sigma.size = s.sigsize,
                                         sigma.const = s.sigma,
                                         decide = s.decide,
                                         nested.doll = s.RD)
               
               convergence.list <- list(xi.size = conv.xi,
                                         sigma.size = conv.sigsize,
                                         sigma.const = conv.sigma,
                                         decide = conv.decide,
                                         nested.doll = conv.RD)
                                   
               time.list <- list(xi.size = time.xi,
                                         sigma.size = time.sigsize,
                                         sigma.const = time.sigma,
                                         decide = time.decide,
                                         nested.doll = time.RD)
               
          #########
          # Summarize stuff in a df

          # mod <- names(summary.list)

          # minutes <- sapply(time.list, FUN = "[", 3)

          # columns <- c("xi.int", "xi.break", "sigma", "sigma.int", "sigma.beta", "N", "N.tot")     

          # cols.get <- function (x, cols = columns) {
            # require(plyr)
            # stopifnot(class(x) ==  "summary.mcmc")
            # tmp <- x$quantiles
            # tmp.df <- join(x = data.frame(names = cols),
                               # y = data.frame(names = row.names(tmp), tmp),
                               # by = "names")[,c(1:2,4,6)]
            # return(tmp.df)                     
            # }

            
          # sig.size.get <- function (x) {
               # "sig.size",      
               
               
               
               
               save(simul.list, samp.list, summary.list, convergence.list, time.list, j.dat,
                     file = paste(fnm.base, "_iter_", j, ".R", sep = ""))

}


}}
