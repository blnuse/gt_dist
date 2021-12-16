##  Function to simulate a line transect distance survey.
##
##
##  15 Dec 2017
########################################

###  source('/home/nuse/Projects/GT_stuff/distance/simulations/new_dec2017/new_sim_15dec2017' )


####  Antecedents:  ###################################
rad <- function (x) {x*pi/180}

`%[]%` <- function (x, y) {stopifnot(length(y) == 2); x >= y[1] && x <= y[2] }
`%()%` <- function (x, y) {stopifnot(length(y) == 2); x > y[1] && x < y[2] }
`%[)%` <- function (x, y) {stopifnot(length(y) == 2); x >= y[1] && x < y[2] }
`%(]%` <- function (x, y) {stopifnot(length(y) == 2); x >= y[1] && x <= y[2] }

trans.y <- function (x, trans.num = 1, theta, spacing, y.int) {
  tan(rad(theta))*x + (y.int + trans.num*spacing/cos(rad(theta)))
}

trans.x <- function (y, trans.num = 1, theta, spacing, y.int) {
  (y - (y.int + trans.num*spacing/cos(rad(theta))))/tan(rad(theta))
}


#S.juv <- 0.7
#S.subad <- 0.8
#S.ad <- 0.98
#F <- 5
#                    max.age <- 80
#                    lambda <- 1
#                    A.f <- 306      ## female
#                    k.f <- 0.056
#                    A.m <- 289      ## male
#                    k.m <- 0.061
#                    prop.f <- 0.5
                    

###  generate a very realistic size-class histogram:
stable.age <- function (S.juv, S.subad, S.ad, F,
                    max.age = 80,
                    lambda = 1,
                    A.f = 306,      ## female
                    k.f = 0.056,
                    A.m = 289,      ## male
                    k.m = 0.061,
                    prop.f = 0.5) {


#####   VITAL RATES
    S1 <- NA                ##  Hatchling (& egg) survival; S1.calc will actually be used in the projection matrix, A
    S2 <- S.juv       ##  Juvenile (soft-shell) survival (until 15cm)
    S3 <- S.subad   ## Subadult suvival.
#    E3 <- ifelse(E3.cap < 0, -runif(1, 0, abs(E3.cap)), runif(1, 0, E3.cap))   ## Subadult dispersal.
#    A3 <- S3 - E3         ## Apparent subadult survival.
    A3 <- S3
    S4 <- S.ad     ## Adult survival.

#####   MATURITY
    maturity.f <- round(-log(1-((230-50)/(A.f-50)))/k.f)
    subadult.f <- round(-log(1-((150-50)/(A.f-50)))/k.f)
    maturity.m <- round(-log(1-((220-50)/(A.m-50)))/k.m)
    subadult.m <- round(-log(1-((150-50)/(A.m-50)))/k.m)
 
 
#####   BUILD MATRIX

    S.vect.f <- c(S1, rep(S2, subadult.f-2), rep(A3, maturity.f-subadult.f), rep(S4, max.age-maturity.f))
    S.vect.m <- c(S1, rep(S2, subadult.m-2), rep(A3, maturity.m-subadult.m), rep(S4, max.age-maturity.m))
    S.vect <- S.vect.f*prop.f + S.vect.m*(1-prop.f)
    
    F.vect <- c(rep(0, maturity.f-1), rep(F, max.age-maturity.f+1))


    ex <- seq(-1, -max.age)
    lambda.exp <- lambda^ex

    ind.list <- list()
    for (i in 2:(max.age-1)) {
      ind.list[[i-1]] <- seq(2, i)
    }

    S.func <- function(x) { prod(S.vect[x]) }
    S.prod <- sapply(ind.list, FUN = S.func, simplify = T)

    numerator <- c(F.vect[1:2], S.prod * F.vect[3:max.age])

    S1.calc <- as.numeric(1/(numerator %*% lambda.exp))

    if (S1.calc > 1) { return(paste0("Whoops, S1.calc = ", round(S1.calc,2), "!")); break }
    
    S.vect.calc <- c(S1.calc, S.vect[-1])
    A.bottom <- cbind(diag(S.vect.calc), 0)

    A <- rbind(F.vect, A.bottom)
    
    eigen.vect <- as.numeric(eigen(A)$vectors[,1])
    stable.age <- eigen.vect/eigen.vect[1]        ###  so "mult" value now represents the number of HATCHLINGS!
    
    out.lst <- mget(c("stable.age", "A.f", "k.f", "A.m", "k.m", "prop.f", "maturity.f", "subadult.f", "maturity.m", "subadult.m"))
    class(out.lst) <- "age.lst"
    
    return(out.lst)
}    

#stable.size <- function(S.jv, S.sa, S.ad, F, mult, bin.width, n.burr.sa = 3, n.burr.ad = 4, ...) {
stable.size <- function(mult, bin.width, n.burr.sa = 3, n.burr.ad = 4, ...) {
  
  lst <- stable.age(...)
  stopifnot(class(lst) == "age.lst")
#  if (class(lst) != "age.lst") {return(lst)}
  
  size.key.f <- with(lst, 50 + (A.f-50)*(1-exp(-k.f*0:(length(stable.age)-1))))
  size.key.m <- with(lst, 50 + (A.m-50)*(1-exp(-k.m*0:(length(stable.age)-1))))
  
  zeros <- rep(0, length(lst$stable.age))

  tween <- function(x, y, burr.mult) {
    ifelse(x %[)% y, burr.mult, 1)
  }
  
  sa.vect.f <- with(lst, sapply(1:length(stable.age), FUN = tween, y = c(subadult.f, maturity.f), burr.mult = n.burr.sa))
  ad.vect.f <- with(lst, sapply(1:length(stable.age), FUN = tween, y = c(maturity.f, Inf), burr.mult = n.burr.ad))
  all.vect.f <- sa.vect.f *  ad.vect.f
 
  sa.vect.m <- with(lst, sapply(1:length(stable.age), FUN = tween, y = c(subadult.m, maturity.m), burr.mult = n.burr.sa))
  ad.vect.m <- with(lst, sapply(1:length(stable.age), FUN = tween, y = c(maturity.m, Inf), burr.mult = n.burr.ad))
  all.vect.m <- sa.vect.m *  ad.vect.m
 
  
  age.f <- with(lst, mult*stable.age*prop.f)
  age.m <- with(lst, mult*stable.age*(1-prop.f))

  burrs.f <- age.f * all.vect.f
  burrs.m <- age.m * all.vect.m
  
  bin.num.f <- ceiling(size.key.f/bin.width)
  bin.num.m <- ceiling(size.key.m/bin.width)  
  
  agg.f <- aggregate(burrs.f, by = list(bin.num.f), FUN = sum)
  agg.m <- aggregate(burrs.m, by = list(bin.num.m), FUN = sum)
  
  add.mat <- matrix(0, ncol = 2, nrow = max(agg.f$Group.1, agg.m$Group.1))
  
  add.mat[agg.f$Group.1,1] <- agg.f$x
  add.mat[agg.m$Group.1,2] <- agg.m$x
  
  hist.all <- rowSums(add.mat)
  
  return(hist.all)
}

### "mult" value represents the number of tortoises in first occupied bin, if method = "well-known"; but HATCHLINGS if using "stable-age"!
ind.sz.lst <- function(mult, method = "well-known", bin.width = 50, rel.abund = c(0,1,0.3,0.25,0.2,0.15,0.3,0.4,0.1), ...) {
   stopifnot(method %in% c("well-known", "stable-age"))
   
   if (method == "stable-age") {
     abund <- round(stable.size(bin.width = bin.width, mult = mult, ...))
     ### stable.size() puts the hatchlings in the 0-50 bin, which is not what I want!
     abund[1:2] <- c(0,sum(abund[1:2]))
   } else {
     abund <- round(rel.abund*mult)
   }
   h.base <- 0:length(abund)*bin.width
   h.setup <- cbind(abund, h.base[-length(h.base)], h.base[-1])     
   inds <- round(unlist(mapply(n = h.setup[,1], min = h.setup[,2], max = h.setup[,3], FUN=runif)))     
   return(inds)
}   
   
###  Examples:  by default, subadults keep 3 burrows each, and adults 4.
#ind.sz.lst(method = "stable-age", mult = 100, S.juv = 0.7, S.subad = 0.8, S.ad = 0.98, F = 5, lambda = 0.98)

#ind.sz.lst(mult = 100)


###  Questions, from Craig Guyer's data:
#1. What is the (size-dependent) rate of occupancy, in burrows?  --> actually, may not be able to answer this from his data....
#2. What is the (size-dependent) relationship of burrow width to carapace length?

      ###  calculate xi and tau: x.shift will always be 50 mm, here.
        broken.stick <- function (x, int, brk, slope, sq=FALSE, x.shift = 50) {
         if (is.nan(slope) || slope == 0 || is.infinite(slope)) {
           y <- 1
         } else {
           if (sq) {
             y <- ifelse(x < brk, int + slope*(x-x.shift)^2, 1)
           } else {
             y <- ifelse(x < brk, int + slope*(x-x.shift), 1)
           }
         }
         return(y)
        }

########################################################               
########################################################               



#----------------------
# Lay out transects:

#edge = 1000
#               n.burr = NULL        ##  true number of burrows
#               intrvl = 100         ##  distance between transects
#               theta = 50       ##  angle of transects
#               sigma = 30
#               xi.50 = 1
#               xi.brk = 50        
#               xi.sq = FALSE
#               tau.50 = 1
#               tau.brk = 50               
#               tau.sq = FALSE




burrSim <- function (
               intrvl,         ##  distance between transects
               sigma,
               mult,
               n.burr = NULL,         ##  true number of burrows
               edge = 1000,
               theta = 50,       ##  angle of transects
               xi.50 = 1,       ##  the intercept of the detection curve at the minimum size = 50mm !!
               xi.brk = 50,        
               xi.sq = FALSE,
               tau.50 = 1,      ##  the intercept of the detection curve depression factor, at the minimum size = 50mm !!
               tau.brk = 50,               
               tau.sq = FALSE,
               plot = FALSE,
               ...) {            ##  you must provide arguments to ind.sz.lst() here !!
               
               
  stopifnot(sigma > 0)
  stopifnot(xi.brk >= 50)               
  stopifnot(xi.50 %[]% c(0,1))               
  stopifnot(tau.brk >= 50)               
  stopifnot(tau.50 %[]% c(0,1))

###  Find y-intercept for lowest line:
  opt <- optimize (
         f = function(y) {abs(trans.y(x = edge, theta = theta, spacing = intrvl, y.int = y))},
         interval = c((-10)*edge, 10*edge)
        )

        y.int.min <- opt$minimum + runif(1, min = 0.2, max = 0.8)*sin(rad(theta))*intrvl
        
        
        ###  Find number of necessary transect lines to cover extent:
        n.lines <- floor((edge -  y.int.min)/(intrvl/cos(rad(theta))))


      ###  find parameters of xi and tau:
        if (xi.sq) {
         xi.beta <- (1 - xi.50)/(xi.brk - 50)^2
        } else {
         xi.beta <- (1 - xi.50)/(xi.brk - 50)
        }
        
        if (tau.sq) { 
         tau.beta <- (1 - tau.50)/(tau.brk - 50)^2
        } else {
         tau.beta <- (1 - tau.50)/(tau.brk - 50)
        }



##################

        ###  get the individuals, who'll be in the population:
        diam.v <- ind.sz.lst(mult, ...)
        stopifnot(all(diam.v >= 50))
        if (is.null(n.burr)) {
          n.burr <- length(diam.v)
          } else {
          if (n.burr > length(diam.v)) {
            warning("Argument n.burr = ", n.burr, " was greater than the number of burrows generated (", length(diam.v), ")!  Using the latter value!")
            n.burr <- length(diam.v)
          } else {
            diam.v <- sample(diam.v, n.burr)
          }
        }
        ###  diam.v <- ind.sz.lst(mult = 100)

        ###  no clumping here!
        loc.x <- runif(n.burr, 0, edge)
        loc.y <- runif(n.burr, 0, edge)
        
        xi <- broken.stick(diam.v, brk = xi.brk, int = xi.50, slope = xi.beta, sq =  xi.sq)
        tau <- broken.stick(diam.v, brk = tau.brk, int = tau.50, slope = tau.beta, sq =  tau.sq)
        
        ##  find distance of each burrow to nearest transect, and which that is!
        ##  d = abs(y)*cos(theta.ch.rad)
        diff.mat <- matrix(NA, ncol = n.lines, nrow = n.burr)
        for (i in 1:n.lines) {
         diff.mat[,i] <- abs(loc.y - trans.y(x = loc.x, trans.num = i, theta = theta, spacing = intrvl, y.int = y.int.min)) * cos(rad(theta))
        }
               
        near.trans <- apply(diff.mat, MAR = 1, which.min)
        min.dist <- apply(diff.mat, MAR = 1, min)
        
        p <- xi*exp(-(min.dist^2/(2*tau*sigma^2)))
        
        ###  Now do sampling!
        found <- rbinom(n = p, size = 1, prob = p)
        
        ###  Find the area actually sampled, given Bx
        Bx <- min(max(min.dist[found == 1]), 0.5*intrvl)

               ###  This is approximately the total length minus, for each transect, Bx^2*(tan(theta.ch.rad) + 1/tan(theta.ch.rad))
        
        ####  Are any of the corners within a transects' strip?
        #  d = abs(y)*cos(theta.ch.rad)
        corn.mat <- matrix(NA, ncol = n.lines, nrow = 4)
        
        ### Keep track of corners thus: 1 -- Lower Left, 2 -- Upper Left, 3 -- Upper Right, 4 -- Lower Right.
        corn.x <- c(0, 0, edge, edge)
        corn.y <- c(0, edge, edge, 0)
        for (i in 1:n.lines) {
         corn.mat[,i] <- abs(corn.y - trans.y(x = corn.x, trans.num = i, theta = theta, spacing = intrvl, y.int = y.int.min)) * cos(rad(theta))
        }
        
        under.Bx <- corn.mat < Bx
        which.corn <- apply(under.Bx, MAR=2, FUN=which)
        ####  LEFT OFF HERE !!!!
        ####  LEFT OFF HERE !!!!
                ####  LEFT OFF HERE !!!!
                        ####  LEFT OFF HERE !!!!
                                ####  LEFT OFF HERE !!!!
        dist.corn <- apply(corn.mat, MAR = 1, FUN = `<`, Bx)
        

               
               beg.x <- beg.y <- end.x <- end.y <- len <- hits.side <- len <- hits.side <- overest <- NA
               for (i in 1:n.lines) { 
                beg.x[i] <- min ( max( c(0, trans.x(y = 0, theta = theta, spacing = intrvl, y.int =  y.int.min, trans.num = i)) ), edge)
                beg.y[i] <- min( max( c(0, trans.y(x = 0, theta = theta, spacing = intrvl, y.int =  y.int.min, trans.num = i)) ), edge)
                end.x[i] <- min ( max( c(0, trans.x(y = edge, theta = theta, spacing = intrvl, y.int =  y.int.min, trans.num = i)) ), edge)
                end.y[i] <- min( max( c(0, trans.y(x = edge, theta = theta, spacing = intrvl, y.int =  y.int.min, trans.num = i)) ), edge)
               
                len[i] <- sqrt((beg.x[i] - end.x[i])^2 + (beg.y[i] - end.y[i])^2)
                if (beg.x[i] == 0) {
                      hits.side[i] <- ifelse(end.x[i] != edge, "both", "side-side")
                      } else {
                      hits.side[i] <- ifelse(end.x[i] != edge, "top-bottom", "both")
                      }
                
                
                #####  This NEEDS TO BE CHECKED!!! --> Incorporate hits.side !!!! ---> ok, done: it's still not perfect, missing some extra area when a transect is near a corner!
                if (hits.side[i] == "both") {
                      overest[i] <- 0.5*Bx^2*(tan(rad(theta)) + 1/tan(rad(theta)))
                      ###  LEFT OFF HERE  !!!!!!   if (length(which.corn[[i]] 
                      }
                if (hits.side[i] == "top-bottom") {
                      overest[i] <- Bx^2*(1/tan(rad(theta)))
                      }        
                if (hits.side[i] == "side-side") {
                      overest[i] <- Bx^2*(tan(rad(theta)))
                      }        

               }        
        
        total.len <- sum(len)
        area <- 2*Bx*total.len - sum(overest) 
        
        ###  Note well that the following df of burrows is filtered by Bx !!!
        all.burr <- as.data.frame(cbind(diam.v, loc.x, loc.y, near.trans, min.dist, p, found)[order(near.trans, diam.v),])
        names(all.burr) <- c("diam", "loc.x", "loc.y", "trans", "dist", "p", "found")
        all.burr$strip <- ifelse(all.burr$dist < Bx, 1, 0)
        
        ###  Info to allow plotting survey lines:
        survey <- mget(c("edge", "n.lines", "theta", "intrvl", "y.int.min", "total.len", "Bx", "area"))
        
        N <- dim(all.burr)[1]
        N.strip <- sum(all.burr$strip)
        N.obs <- sum(found)
        
        ###  Estimated parameters:
        par <- mget(c("N", "N.strip", "N.obs", "sigma", "xi.50", "xi.brk", "xi.sq", "xi.beta", "tau.50", "tau.brk", "tau.sq", "tau.beta"))
        
        want.out <- c("par", "all.burr", "survey")  
        out.lst <- mget(want.out)
        class(out.lst) <- "burrSim"
        return(out.lst)
}  ## end sim()

####left off testing sim(), with method="stable-age":
#tst <- do.call(sim, c(method="stable-age", mult = 100, sigma = 20, intrvl = 100, s))
#s <- list(S.juv = 0.8, S.subad = 0.85, S.ad = 0.95, F = 4, n.burr = 200)

