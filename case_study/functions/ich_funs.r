##  Function to run JAGS models on Ichauway pre- and post-burn data.
##
##
##  Oct 2019, Dec 2021
##############################################


##  Notes:
##  ----------
##  1)  This function is quite idiosyncratic.  It uses 4 cores to run MCMC chains in parallel; to change this, alter the function definition here and reload.
##  2)  Likewise for file paths.
##  3)  Requires: 1) JAGS to be installed on the system; 2) these R packages: {rjags}, {parallel}. 
##  4)  The function will not work on Windows systems as written, because it uses the function mclapply() (in package {parallel}) which relies on forking of processes.  To use on Windows, replace mclapply() with parLapply() (which takes considerably more setup; see: https://stackoverflow.com/questions/17196261/understanding-the-differences-between-mclapply-and-parlapply-in-r); or, set n.cores = 1 when calling function runfun() (you'll only be able to run one MCMC chain this way).




runfun <- function(modl,                ##   Character string: name of the model, for use in filenames, etc.
                   datl,                ##   List: necessary vectors and matrices to be passed to JAGS.
                   initl,               ##   List: intital values for JAGS.
                   iter.samp,           ##   Scalar: number of iterations over which to sample posterior distributions of interest (i.e., 'monitor' in JAGS parlance).
                   vars,                ##   Character vector: nodes to monitor.
                   pre = T,             ##   Logical: should the pre-burn survey data be included?
                   post = F,            ##   Logical: should the post-burn survey data be included?
                   name.suffix = "",    ##   Character string: suffix to be added to filenames, etc.
                   n.cores = 4          ##   Number of cores to use, for parallel processing.
                   ) { 

      stopifnot(pre || post)

      require(parallel); options(mc.cores = n.cores)

      require(rjags); load.module("glm"); load.module("dic"); load.module("lecuyer")

  jags.fun <- function (l) {
    init.fun <- function() {
      list(.RNG.state = l$.RNG.state,
           .RNG.name = l$.RNG.name,
           w = unlist(initl$w))
    }
  
    jm <- jags.model(file=l$file.nm, inits=init.fun(), data=l$dat, n.chains=1) #inits=init.lst,
    update(jm, n.iter=l$iter.burn)
    js <- coda.samples(jm, var=l$var, n.iter=l$iter.samp)
    return(js)
  }  ##  END jags.fun() definition.


  rnd <- parallel.seeds("lecuyer::RngStream", n.cores)

  jl.unit <- list(iter.burn = iter.samp,
                  iter.samp = iter.samp,
                  var = vars,
                  dat = datl,
                  file.nm = paste0("case_study/jags/", modl))
  
  jl <- lapply(rnd, FUN=c, jl.unit)                 
                  
time <- system.time(
  out <- mclapply(jl, FUN=jags.fun)  ## mclapply() parallelizes this!  Note this function is not available on Windows.
)
  
  out <- mcmc.list(sapply(out, FUN = `[[`, 1, simplify=F))

#  return(out)
#}
  s <- summary(out)

  save(datl, time, out, s, file = paste0("case_study/jags/output", modl, name.suffix, ".R"))

  if (pre && !post) { p.pre <- s$quantiles[grep(x = rownames(s$quantiles), pattern = "p.pred", fixed = T), c(1,3,5)] }
  if (post && !pre) { p.post <- s$quantiles[grep(x = rownames(s$quantiles), pattern = "p.pred", fixed = T), c(1,3,5)] }
  if (pre && post) {
    p.pre <- s$quantiles[grep(x = rownames(s$quantiles), pattern = "p\\.pred\\[[0-9]{1,2},1\\]", fixed = F), c(1,3,5)]
    p.post <- s$quantiles[grep(x = rownames(s$quantiles), pattern = "p\\.pred\\[[0-9]{1,2},2\\]", fixed = F), c(1,3,5)]
  }

png(file=paste0("case_study/figs/", modl, name.suffix, ".png"), res=300, height = 7, width = 10, units="in")
par(mfrow = c(2,2), oma = c(3.5,3.5,1,1), mar = c(0,0,0,0), xaxs="i", yaxs = "i")

  lg <- rgb(t(col2rgb("lightgrey")/255), alpha = 0.5)
  lc <- rgb(t(col2rgb("lightcyan2")/255), alpha = 0.5)
    
for (i in 1:4) { 

with(datl,
  { tst <- p.vals[,1] == 10*i

    if (pre) { poly.pre <- cbind(c(p.vals[tst,2], rev(p.vals[tst,2])), c(p.pre[tst,1], rev(p.pre[tst,3]))) } 
    if (post) { poly.post <- cbind(c(p.vals[tst,2], rev(p.vals[tst,2])), c(p.post[tst,1], rev(p.post[tst,3]))) }

    if (pre) { plot(p.vals[tst,2], p.pre[tst,2], ylim = c(0,1), type="n", main = "", ylab = "", xlab = "", xaxt = "n", yaxt = "n")
    } else {
      plot(p.vals[tst,2], p.post[tst,2], ylim = c(0,1), type="n", main = "", ylab = "", xlab = "", xaxt = "n", yaxt = "n")
    }
    
    if (pre) { polygon(poly.pre, col = lg, border = NA) }
    if (post) { polygon(poly.post, col = lc, border = NA) }
  
    if (pre) { lines(p.vals[tst,2], p.pre[tst,2], col = "black") }
    if (post) { lines(p.vals[tst,2], p.post[tst,2], col = "blue") }

    u <- par("usr")
    text(x = u[1] + 0.85*(u[2]-u[1]), y = u[3] + 0.85*(u[4]-u[3]), labels = paste(10*i, "cm"), cex = 1.3, font = 2)
    axis(side = 1, at = if (i %in% c(1,3)) {seq(0,20,5)} else {seq(5,25,5)}, labels = ifelse(i %in% 1:2, F, T), tcl = 0.5, padj = -1, lwd = 0, lwd.ticks = 1)
    axis(side = 2, at = if (i %in% c(1,2)) {seq(0.2,1,0.2)} else {seq(0,0.8,0.2)}, labels = ifelse(i %in% c(1,3), T, F), tcl = 0.5, padj = 1, lwd = 0, lwd.ticks = 1)
    axis(side = 3, labels = F, tcl = 0.5, lwd = 0, lwd.ticks = 1)
    axis(side = 4, labels = F, tcl = 0.5, lwd = 0, lwd.ticks = 1)
  })


      
  if (i == 4) {
    l <- legend(x = 3, y = 0.1, yjust = 0, lty = 1, legend = c("Pre-burn", "Post-burn"), plot = F)
    if (pre) { leg.poly.pre <- cbind(l$text$x[1] + c(-0.3, -0.1, -0.1, -0.3)*l$rect$w,
                          l$text$y[1] + c(0.1, 0.1, -0.1, -0.1)*l$rect$h)
               polygon(leg.poly.pre, col = lg, border = NA) }
    if (post) { leg.poly.post <- cbind(l$text$x[2] + c(-0.3, -0.1, -0.1, -0.3)*l$rect$w, l$text$y[2] + c(0.1, 0.1, -0.1, -0.1)*l$rect$h)
                polygon(leg.poly.post, col = lc, border = NA) }

    legend(x = 3, y = 0.1, yjust = 0, lty = 1, col = c("black", "blue"), legend = c("Pre-burn", "Post-burn"), bg = NULL)
  }
}

mtext(side = 1, line = 2.5, text = "Distance (m)", outer = T)
mtext(side = 2, line = 2.5, text = "Detection probability", outer = T)


dev.off()

return(mget(c("time", "s")))
}

