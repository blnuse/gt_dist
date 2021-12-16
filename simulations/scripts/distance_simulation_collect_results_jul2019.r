#  Collecting, inspecting results of model-fitting to distance simulations
#  
#
#
#  originally produced only figures 18 May 2015 (see'nov2016' version); altered to produce a table July 2019 (at bottom)
##############################################


## OK, I see that the summary list 'simul.res' was saved as "simul_res_2016-11-04.R" -- so I just need to load that...
load("/home/nuse/Projects/GT_stuff/distance/simulations/work_Jul2019/sim_modfits_from_Obelisk/simul_res_2016-11-04.R")

####  The base file path is:
#base.path <- "/home/nuse/Projects/GT_stuff/distance/simulations/work_Jul2019/sim_modfits_from_Obelisk"

####  The folders containing the two types of simulations are:
#fold.1 <- c("sigma_size", "xi_size", "nested_doll")

####  Within each of those are:

#fold.2 <- paste("interval", int.vals, sep = "_")

# fold.3 <- "2015_07_29"
gen.mods <- gsub(x = names(simul.res), pattern = "_", replace = ".", fixed = T)
int.vals <- gsub(x = names(simul.res[[1]]), pattern = ".*([1|2]00)", replace = "\\1")
pop.sz <- gsub(x = names(simul.res[[1]][[1]]), pattern = ".*([84215]0+)", replace = "\\1")

fit.old <- names(simul.res[[1]][[1]][[1]][1:5])
fit.old <- fit.old[!(fit.old == "decide")]
fit.new <- tolower(gsub(gsub(x = fit.old, pattern = "_", replace = ".", fixed = T), pattern = "Russian", replace = "nested", fixed = T))


res.arr <- array(NA, dim = c(gen.mods = length(gen.mods), int.vals = length(int.vals), pop.sz = length(pop.sz), fit.mods = length(fit.old), stat = 5), dimnames = list(gen.mods = gen.mods, interval = int.vals, pop.size = pop.sz, fit.mods = fit.new, stat = c("bias.med", "bias.sd", "bias.5th", "bias.95th", "coverage")))
pop.int <- as.integer(pop.sz)
for (i in 1:dim(res.arr)[1]) {
  for (j in 1:dim(res.arr)[2]) {
    for (k in 1:dim(res.arr)[3]) {
      for (l in 1:dim(res.arr)[4]) {
        tmp <- simul.res[[i]][[j]][[k]][[fit.old[l]]]
        res.arr[i,j,k,l,"bias.med"] <- median(100*tmp$N.tot.diff/pop.int[k])
        res.arr[i,j,k,l,"bias.sd"] <- sd(100*tmp$N.tot.diff/pop.int[k])
        res.arr[i,j,k,l,"bias.5th"] <- quantile(100*tmp$N.tot.diff/pop.int[k], prob = 0.05)
        res.arr[i,j,k,l,"bias.95th"] <- quantile(100*tmp$N.tot.diff/pop.int[k], prob = 0.95)
        res.arr[i,j,k,l,"coverage"] <- sum(tmp$N.tot.inCrI)
      }
    }
  }
}

library(rgl)
library(shape)
library(viridis)

vp <- viridis_pal(option = "C")
col.vec <- vp(4)
#vp.alpha <- viridis_pal(alpha = 0.5, option = "D")
#alpha.vec <- vp.alpha(4)

bg.gray <- "gray75"
fg.gray <- "gray90"

x.rng <- c(0,100)
y.rng <- 1.1*range(res.arr[,,,,"bias.med"])
x.at <- pretty(x.rng)
y.at <- pretty(y.rng)


###  Bias med:
#x11(width = 15, height = 10)

png(file = "/home/nuse/Projects/GT_stuff/distance/simulations/work_Jul2019/figs/simul_sum_v4.png", units = "in", width = 10, height = 6, res = 400)

    ###  Version 1:
#par(mfrow = c(2,3), mar = c(0,0,0,0), oma = c(6,6,4,6.5), cex.axis = 1.3)
    ###  Version 2:
par(mfrow = c(2,3), mar = c(0.5,0.5,0.5,0.5), oma = c(4.5,5,3,6), cex.axis = 1.3)

ct <- 1 
for (i in 1:dim(res.arr)["int.vals"]) {
for (g in 1:dim(res.arr)["gen.mods"]) {
    ###  Version 1:
#  plot(x = res.arr[g,i,,"xi.size","coverage"], y = res.arr[g,i,,"xi.size","bias.med"], type = "n", xlim = x.rng, ylim = y.rng, xaxt = "n", yaxt = "n", xaxs="i", yaxs = "i", bty = "n")
    ###  Version 2:
  plot(x = res.arr[g,i,,"xi.size","coverage"], y = res.arr[g,i,,"xi.size","bias.med"], type = "n", xaxt = "n", yaxt = "n", xlim = x.rng, ylim = y.rng, bty = "n")
  u <- par("usr")
    ###  Version 4:
  rect(u[1], u[3], u[2], u[4], col = bg.gray)
    ###  Version 1:
#  sapply(x.at[-c(1,length(x.at))], FUN = function (s) abline(v = s, col = "gray90"))
    ###  Version 2:
  sapply(x.at, FUN = function (s) abline(v = s, col = fg.gray))
  sapply(y.at, FUN = function (s) abline(h = s, col = fg.gray))
#  abline(h = 0, col = "gray70")
#  abline(v = 100, col = "gray70")
  points(x = 100, y = 0, pch = 8, col = "white", cex = 1.5)
    ###  Version 1:
#  axis(side = 1, at = {if (ct %% 3 == 1) {tmp <- x.at[-length(x.at)]}; if (ct %% 3 == 2) {tmp <- x.at[-c(1,length(x.at))]}; if (ct %% 3 == 0) {tmp <- x.at[-1]}; tmp}, labels = ifelse(ct > 3, T, F), tcl = 0.5, padj = -1, lwd = 0, lwd.ticks = 0)
#  axis(side = 2, at = {if (ct <= 3) {tmp <- y.at[-1]}; if (ct > 3) {tmp <- y.at[-length(y.at)]}; tmp}, labels = ifelse(ct %% 3 == 1, T, F), tcl = 0.5, padj = 1, lwd = 0, lwd.ticks = 0)
#  axis(side = 3, at = {if (ct %% 3 == 1) {tmp <- x.at[-length(x.at)]}; if (ct %% 3 == 2) {tmp <- x.at[-c(1,length(x.at))]}; if (ct %% 3 == 0) {tmp <- x.at[-1]}; tmp}, labels = F, tcl = 0.5, lwd = 0, lwd.ticks = 0)
#  axis(side = 4, at = {if (ct <= 3) {tmp <- y.at[-1]}; if (ct > 3) {tmp <- y.at[-length(y.at)]}; tmp}, labels = F, tcl = 0.5, lwd = 0, lwd.ticks = 0)
    ###  Version 2:
  axis(side = 1, at = x.at, labels = ifelse(ct > 3, T, F), tcl = 0.5, padj = -1, lwd = 0, lwd.ticks = 0)
  axis(side = 2, at = {if (ct <= 3) {tmp <- y.at[-1]}; if (ct > 3) {tmp <- y.at[-length(y.at)]}; tmp}, labels = ifelse(ct %% 3 == 1, T, F), tcl = 0.5, padj = 1, lwd = 0, lwd.ticks = 0)
  axis(side = 3, at = x.at, labels = F, tcl = 0.5, lwd = 0, lwd.ticks = 0)
  axis(side = 4, at = {if (ct <= 3) {tmp <- y.at[-1]}; if (ct > 3) {tmp <- y.at[-length(y.at)]}; tmp}, labels = F, tcl = 0.5, lwd = 0, lwd.ticks = 0)

    ###  Version 1:
#  box(bty = {if (ct == 1) {tmp <- "c"}; if (ct == 2) {tmp <- "o"}; if (ct == 3) {tmp <- "]"}; if (ct == 4) {tmp <- "l"}; if (ct == 5) {tmp <- "l"}; if (ct == 6) {tmp <- "u"}; tmp})
    ###  Version 2:
  box()
  if (ct == 6) {
#    u <- par("usr")
    l <- legend(x = "topleft", yjust = 0, lty = 1, pch = 15, legend = c("Intercept-varying", "Scale-varying", "Conventional DA", "Nested doll"), bg = bg.gray, col = col.vec, cex = 1.3, plot = F)
    rect(xleft = l$rect$left,
         xright = l$rect$left + l$rect$w,
         ybottom = l$rect$top - l$rect$h,
         ytop = l$rect$top,
         col = "white", border = NA)
    rect(xleft = l$rect$left,
         xright = l$rect$left + 0.23*l$rect$w,
         ybottom = l$rect$top - l$rect$h,
         ytop = l$rect$top,
         col = bg.gray, border = NA)
legend(x = "topleft", yjust = 0, lty = 1, pch = 15, legend = c("Intercept-varying", "Scale-varying", "Conventional DA", "Nested doll"), bg = "transparent", col = col.vec, cex = 1.3)
  }  
  
  for (m in 1:dim(res.arr)["fit.mods"]) {
    for (p in dim(res.arr)["pop.sz"]:2) {
    ###  Version 2:
      Arrows(x0 = res.arr[g,i,p, m,"coverage"],
             x1 = res.arr[g,i,p-1, m,"coverage"],
             y0 = res.arr[g,i,p, m,"bias.med"],
             y1 = res.arr[g,i,p-1, m,"bias.med"],
             col = col.vec[m], arr.type = "triangle",
             arr.adj = 0,
             arr.length = ifelse(p == 2, 0.3, 0), lwd = 1.5)
             
    ###  Version 1:
#      Arrows(x0 = res.arr[g,i,p, m,"coverage"],
#             x1 = res.arr[g,i,p-1, m,"coverage"],
#             y0 = res.arr[g,i,p, m,"bias.med"],
#             y1 = res.arr[g,i,p-1, m,"bias.med"], col = m, arr.type = "triangle", arr.adj = 1, arr.length = 0.3)
      }
    points(x = res.arr[g,i,,m,"coverage"], y = res.arr[g,i,,m,"bias.med"], pch = ifelse(dimnames(res.arr)$gen.mods[g] == dimnames(res.arr)$fit.mods[m], 0, 15), col = col.vec[m], bg = bg.gray, cex = 1.5)   ## 22
    }
    ct <- ct + 1
  }  
}

mtext(side = 1, text = "Coverage (%)", line = 2.5, outer = T, cex = 1.3, font = 2) 
mtext(side = 2, text = "Median deviation (%)", line = 2.4, outer = T, font = 2, cex = 1.3)
mtext(side = 3, text = c("Scale-varying", "Intercept-varying", "Nested doll"), at = c(1,3,5)/6, line = 1, outer = T, font = 3, cex = 1.3)
mtext(side = 4, text = c("100 m", "200 m"), at = c(1,3)/4, line = 0.5, outer = T, font = 3, cex = 1.3, las = 1)

dev.off()




####  Bias med and sd:
#plot3d(x = matrix(res.arr["xi.size", "100",i, m,c("coverage","bias.med","bias.sd")], nrow = 1), type = "n", xlim = range(res.arr["xi.size","100",,,"coverage"]), ylim = range(res.arr["xi.size","100",,,"bias.med"]), zlim = range(res.arr["xi.size","100",,,"bias.sd"]), xlab = "Coverage", ylab = "Median bias", zlab = "Sd(bias)")

#for (m in 1:length(fit.mods)) {
#  for (i in length(pop.sz):2) {
#    arrow3d(p0 = res.arr["xi.size", "100",i, m,c("coverage","bias.med","bias.sd")],
#           p1 = res.arr["xi.size", "100",i-1, m,c("coverage","bias.med","bias.sd")],
#           barblen = c(10,10,10),
#           col = m)

#  }
#}  



######  TABLE-MAKE!
## with units like:

### OK, do two tables: one with transects spaced 100m, the other with 200m.
			Gen mod
			
Fit mod	50 ha^-1	Coverage	Bias [med (lcr, ucr)]
	100
	200
	400
	800
------------------------

###  100-m spacing
cov.100 <- matrix(aperm(res.arr[,"100",,,"coverage"], perm = c(2,3,1)), ncol = dim(res.arr)["gen.mods"])
cov.200 <- matrix(aperm(res.arr[,"200",,,"coverage"], perm = c(2,3,1)), ncol = dim(res.arr)["gen.mods"])

bias.100 <- matrix(format(aperm(res.arr[,"100",,,"bias.med"], perm = c(2,3,1)), digits = 1, nsmall = 1, trim = T), ncol = dim(res.arr)["gen.mods"])
bias.int.100 <- matrix(paste0("[",format(aperm(res.arr[,"100",,,"bias.5th"], perm = c(2,3,1)), digits = 1, nsmall = 1, trim = T), ",",
                          format(aperm(res.arr[,"100",,,"bias.95th"], perm = c(2,3,1)), digits = 0, nsmall = 1, trim = T), "]"),
                   ncol = dim(res.arr)["gen.mods"])

bias.200 <- matrix(format(aperm(res.arr[,"200",,,"bias.med"], perm = c(2,3,1)), digits = 1, nsmall = 1, trim = T), ncol = dim(res.arr)["gen.mods"])
bias.int.200 <- matrix(paste0("[",format(aperm(res.arr[,"200",,,"bias.5th"], perm = c(2,3,1)), digits = 1, nsmall = 1, trim = T), ",",
                          format(aperm(res.arr[,"200",,,"bias.95th"], perm = c(2,3,1)), digits = 0, nsmall = 1, trim = T), "]"),
                   ncol = dim(res.arr)["gen.mods"])

tab.100 <- matrix(rbind(cov.100, bias.100, bias.int.100), nrow = nrow(cov.100))
colnames(tab.100) <- rep(c("Coverage", "Bias", "Credible Interval"), length.out = ncol(tab.100))

d <- dimnames(res.arr)
pop.lab <- sapply(as.integer(d$pop.size)/100, FUN = function (s) { ns <- ifelse(s < 1, 1, 0); format(s, digits = 0, nsmall = ns)})
fit.mod.col <- c(t(matrix(c(d$fit.mods, rep("", length(d$fit.mods) * (length(d$pop.size)-1))), nrow = length(d$fit.mods))))

df.100 <- data.frame(fit.mod = fit.mod.col, pop.dens = rep(pop.lab, times = length(d$fit.mods)), as.data.frame(tab.100))

tab.200 <- matrix(rbind(cov.200, bias.200, bias.int.200), nrow = nrow(cov.200))
colnames(tab.200) <- rep(c("Coverage", "Bias", "Credible Interval"), length.out = ncol(tab.200))
df.200 <- data.frame(fit.mod = fit.mod.col, pop.dens = rep(pop.lab, times = length(d$fit.mods)), as.data.frame(tab.200))

library(xtable)

x.100 <- xtable(df.100)
x.200 <- xtable(df.200)

print.xtable(x.100, include.rownames = F, booktabs = T)
print.xtable(x.200, include.rownames = F, booktabs = T)


#rect(xleft = 0.5, xright = 29.5, ybottom = 93, ytop = 105, col = "grey95", border = F, xpd = T)

#text(x = seq(3,27,6), y = 98.5, labels = paste(pop.sz/100, "/ ha"))
#abline(v = seq(6,24,6))

#  
#axis(side = 2, las = 2)
#mtext(side = 2, text = expression(paste("deviation from ", italic(N), " (%)")), line = 2.5)

#leg.place <- c(1, 5, 9, 13, 17)
#points(x = leg.place, , y = rep(-112,5), bg = cols, pch = 22, cex = 2.2, xpd = TRUE)
#y.adj <- c(-1,0,0,-1,0)*0.4
#text(x = leg.place, y = rep(-125,5)+y.adj, pch = 19, cex = 1.1, xpd = TRUE, labels = c("Intercept-var.",               #expression(paste(xi(italic(d)))),
#                                                                                       "Scale-var.",                     # expression(paste(sigma(italic(d)))),
#                                                                                       "Conv. DA",                # expression(paste(sigma)),
#                                                                                       "Scale/Intercept",           # expression(paste(xi(italic(d)), " or ", sigma(italic(d)))),
#                                                                                       "Nested doll"))

#text(x = 22, y = -115, labels = "generating model = Scale-var.", xpd = T, pos = 4)          ## expression(paste("generating model = ", sigma(italic(d))))
#text(x = 22, y = -125, labels = "transect interval = 100m", xpd = T, pos = 4)

#dev.off()




## ####  Now sigma_size, interval_200
## x11(width = 11, height = 5)
#png(file = "/media/nuse/E/GT_project/Analysis/Distance/simulation/figures/nov2016/sig_size_int200.png", width = 11, height = 5, res = 300, units = "in")
#par(mar = c(3,4,2,1))
#plot.new()     
#plot.window(xlim = c(1,30), ylim = c(-95,95))


## abline(h = 0)
#segments(x0 = 0, x1 = 30, y0 = 0)
#rect(xleft = 0.5, xright = 29.5, ybottom = -102, ytop = -90, col = "grey95", border = F, xpd = T)

#k <- 1
#for (i in 1:length(pop.sz)) {
#  for (j in 1:length(mod.names)) {
#    tmp.y <- 100*simul.res$sigma_size$interval_200[[i]][[j]]$N.tot.diff/pop.sz[i]
#     tmp.sum <- sum(simul.res$sigma_size$interval_200[[i]][[j]]$N.tot.inCrI)
#    # points(x = rep(k, 100), y = tmp.y, col = cols[j], pch = 19, cex = 1.3)
#     boxplot(x = tmp.y, at = k, add = T, col = cols[j], axes = F, yaxt = "n", ylim = c(-50,50))
#     text(x = k, y = -95, labels = tmp.sum)
#    k <- k + 1
#  }
#  k <- k + 1
#  }

#rect(xleft = 0.5, xright = 29.5, ybottom = 93, ytop = 105, col = "grey95", border = F, xpd = T)

#text(x = seq(3,27,6), y = 98.5, labels = paste(pop.sz/100, "/ ha"))
#abline(v = seq(6,24,6))

#  
#axis(side = 2, las = 2)
#mtext(side = 2, text = expression(paste("deviation from ", italic(N), " (%)")), line = 2.5)

#leg.place <- c(1, 5, 9, 13, 17)
#points(x = leg.place, , y = rep(-112,5), bg = cols, pch = 22, cex = 2.2, xpd = TRUE)
#y.adj <- c(-1,0,0,-1,0)*0.4
#text(x = leg.place, y = rep(-125,5)+y.adj, pch = 19, cex = 1.1, xpd = TRUE, labels = c("Intercept-var.",               #expression(paste(xi(italic(d)))),
#                                                                                 "Scale-var.",                     # expression(paste(sigma(italic(d)))),
#                                                                                                     "Conv. DA",                # expression(paste(sigma)),
#                                                                                                     "Scale/Intercept",           # expression(paste(xi(italic(d)), " or ", sigma(italic(d)))),
#                                                                                                     "Nested doll"))

#text(x = 22, y = -115, labels = "generating model = Scale-var.", xpd = T, pos = 4)          ## expression(paste("generating model = ", sigma(italic(d))))
#text(x = 22, y = -125, labels = "transect interval = 200m", xpd = T, pos = 4)

#dev.off()



#####  And now xi_size, interval_100
## x11(width = 11, height = 5)
#png(file = "/media/nuse/E/GT_project/Analysis/Distance/simulation/figures/nov2016/xi_size_int100.png", width = 11, height = 5, res = 300, units = "in")
#par(mar = c(3,4,2,1))
#plot.new()     
#plot.window(xlim = c(1,30), ylim = c(-95,95))


## abline(h = 0)
#segments(x0 = 0, x1 = 30, y0 = 0)
#rect(xleft = 0.5, xright = 29.5, ybottom = -102, ytop = -90, col = "grey95", border = F, xpd = T)

#k <- 1
#for (i in 1:length(pop.sz)) {
#  for (j in 1:length(mod.names)) {
#    tmp.y <- 100*simul.res$xi_size$interval_100[[i]][[j]]$N.tot.diff/pop.sz[i]
#     tmp.sum <- sum(simul.res$xi_size$interval_100[[i]][[j]]$N.tot.inCrI)
#    # points(x = rep(k, 100), y = tmp.y, col = cols[j], pch = 19, cex = 1.3)
#     boxplot(x = tmp.y, at = k, add = T, col = cols[j], axes = F, yaxt = "n", ylim = c(-50,50))
#     text(x = k, y = -95, labels = tmp.sum)
#    k <- k + 1
#  }
#  k <- k + 1
#  }

#rect(xleft = 0.5, xright = 29.5, ybottom = 93, ytop = 105, col = "grey95", border = F, xpd = T)

#text(x = seq(3,27,6), y = 98.5, labels = paste(pop.sz/100, "/ ha"))
#abline(v = seq(6,24,6))

#  
#axis(side = 2, las = 2)
#mtext(side = 2, text = expression(paste("deviation from ", italic(N), " (%)")), line = 2.5)

#leg.place <- c(1, 5, 9, 13, 17)
#points(x = leg.place, , y = rep(-112,5), bg = cols, pch = 22, cex = 2.2, xpd = TRUE)
#y.adj <- c(-1,0,0,-1,0)*0.4
#text(x = leg.place, y = rep(-125,5)+y.adj, pch = 19, cex = 1.1, xpd = TRUE, labels = c("Intercept-var.",               #expression(paste(xi(italic(d)))),
#                                                                                 "Scale-var.",                     # expression(paste(sigma(italic(d)))),
#                                                                                                     "Conv. DA",                # expression(paste(sigma)),
#                                                                                                     "Scale/Intercept",           # expression(paste(xi(italic(d)), " or ", sigma(italic(d)))),
#                                                                                                     "Nested doll"))

#text(x = 22, y = -115, labels = "generating model = Intercept-var.", xpd = T, pos = 4)          ##  expression(paste("generating model = ", xi(italic(d))))
#text(x = 22, y = -125, labels = "transect interval = 100m", xpd = T, pos = 4)

#dev.off()


#####  Finally, xi_size, interval_200
## x11(width = 11, height = 5)
#png(file = "/media/nuse/E/GT_project/Analysis/Distance/simulation/figures/nov2016/xi_size_int200.png", width = 11, height = 5, res = 300, units = "in")
#par(mar = c(3,4,2,1))
#plot.new()     
#plot.window(xlim = c(1,30), ylim = c(-95,95))


## abline(h = 0)
#segments(x0 = 0, x1 = 30, y0 = 0)
#rect(xleft = 0.5, xright = 29.5, ybottom = -102, ytop = -90, col = "grey95", border = F, xpd = T)

#k <- 1
#for (i in 1:length(pop.sz)) {
#  for (j in 1:length(mod.names)) {
#    tmp.y <- 100*simul.res$xi_size$interval_200[[i]][[j]]$N.tot.diff/pop.sz[i]
#     tmp.sum <- sum(simul.res$xi_size$interval_200[[i]][[j]]$N.tot.inCrI)
#    # points(x = rep(k, 100), y = tmp.y, col = cols[j], pch = 19, cex = 1.3)
#     boxplot(x = tmp.y, at = k, add = T, col = cols[j], axes = F, yaxt = "n", ylim = c(-50,50))
#     text(x = k, y = -95, labels = tmp.sum)
#    k <- k + 1
#  }
#  k <- k + 1
#  }

#rect(xleft = 0.5, xright = 29.5, ybottom = 93, ytop = 105, col = "grey95", border = F, xpd = T)

#text(x = seq(3,27,6), y = 98.5, labels = paste(pop.sz/100, "/ ha"))
#abline(v = seq(6,24,6))

#  
#axis(side = 2, las = 2)
#mtext(side = 2, text = expression(paste("deviation from ", italic(N), " (%)")), line = 2.5)

#leg.place <- c(1, 5, 9, 13, 17)
#points(x = leg.place, , y = rep(-112,5), bg = cols, pch = 22, cex = 2.2, xpd = TRUE)
#y.adj <- c(-1,0,0,-1,0)*0.4
#text(x = leg.place, y = rep(-125,5)+y.adj, pch = 19, cex = 1.1, xpd = TRUE, labels = c("Intercept-var.",               #expression(paste(xi(italic(d)))),
#                                                                                 "Scale-var.",                     # expression(paste(sigma(italic(d)))),
#                                                                                                     "Conv. DA",                # expression(paste(sigma)),
#                                                                                                     "Scale/Intercept",           # expression(paste(xi(italic(d)), " or ", sigma(italic(d)))),
#                                                                                                     "Nested doll"))

#text(x = 22, y = -115, labels = "generating model = Intercept-var.", xpd = T, pos = 4)          ##  expression(paste("generating model = ", xi(italic(d))))
#text(x = 22, y = -125, labels = "transect interval = 200m", xpd = T, pos = 4)

#dev.off()


#####  And now nested_doll, interval_100
## x11(width = 11, height = 5)
#png(file = "/media/nuse/E/GT_project/Analysis/Distance/simulation/figures/nov2016/nested_doll_int100.png", width = 11, height = 5, res = 300, units = "in")
#par(mar = c(3,4,2,1))
#plot.new()     
#plot.window(xlim = c(1,30), ylim = c(-95,95))


## abline(h = 0)
#segments(x0 = 0, x1 = 30, y0 = 0)
#rect(xleft = 0.5, xright = 29.5, ybottom = -102, ytop = -90, col = "grey95", border = F, xpd = T)

#k <- 1
#for (i in 1:length(pop.sz)) {
#  for (j in 1:length(mod.names)) {
#    tmp.y <- 100*simul.res$nested_doll$interval_100[[i]][[j]]$N.tot.diff/pop.sz[i]
#    tmp.sum <- sum(simul.res$nested_doll$interval_100[[i]][[j]]$N.tot.inCrI)
#    # points(x = rep(k, 100), y = tmp.y, col = cols[j], pch = 19, cex = 1.3)
#    boxplot(x = tmp.y, at = k, add = T, col = cols[j], axes = F, yaxt = "n", ylim = c(-50,50))
#    text(x = k, y = -95, labels = tmp.sum)
#    k <- k + 1
#  }
#  k <- k + 1
#  }

#rect(xleft = 0.5, xright = 29.5, ybottom = 93, ytop = 105, col = "grey95", border = F, xpd = T)

#text(x = seq(3,27,6), y = 98.5, labels = paste(pop.sz/100, "/ ha"))
#abline(v = seq(6,24,6))

#  
#axis(side = 2, las = 2)
#mtext(side = 2, text = expression(paste("deviation from ", italic(N), " (%)")), line = 2.5)

#leg.place <- c(1, 5, 9, 13, 17)
#points(x = leg.place, , y = rep(-112,5), bg = cols, pch = 22, cex = 2.2, xpd = TRUE)
#y.adj <- c(-1,0,0,-1,0)*0.4
#text(x = leg.place, y = rep(-125,5)+y.adj, pch = 19, cex = 1.1, xpd = TRUE, labels = c("Intercept-var.",               #expression(paste(xi(italic(d)))),
#                                                                                       "Scale-var.",                     # expression(paste(sigma(italic(d)))),
#                                                                                       "Conv. DA",                # expression(paste(sigma)),
#                                                                                       "Scale/Intercept",           # expression(paste(xi(italic(d)), " or ", sigma(italic(d)))),
#                                                                                       "Nested doll"))

#text(x = 22, y = -115, labels = "generating model = Nested doll", xpd = T, pos = 4)          ##  expression(paste("generating model = ", xi(italic(d))))
#text(x = 22, y = -125, labels = "transect interval = 100m", xpd = T, pos = 4)

#dev.off()




#####  And now nested_doll, interval_200
##x11(width = 11, height = 5)
#png(file = "/media/nuse/E/GT_project/Analysis/Distance/simulation/figures/nov2016/nested_doll_int200.png", width = 11, height = 5, res = 300, units = "in")
#par(mar = c(3,4,2,1))
#plot.new()     
#plot.window(xlim = c(1,30), ylim = c(-95,95))

#segments(x0 = 0, x1 = 30, y0 = 0)
#rect(xleft = 0.5, xright = 29.5, ybottom = -102, ytop = -90, col = "grey95", border = F, xpd = T)

#k <- 1
#for (i in 1:length(pop.sz)) {
#  for (j in 1:length(mod.names)) {
#    tmp.y <- 100*simul.res$nested_doll$interval_200[[i]][[j]]$N.tot.diff/pop.sz[i]
#    tmp.sum <- sum(simul.res$nested_doll$interval_200[[i]][[j]]$N.tot.inCrI)
#    boxplot(x = tmp.y, at = k, add = T, col = cols[j], axes = F, yaxt = "n", ylim = c(-50,50))
#    text(x = k, y = -95, labels = tmp.sum)
#    k <- k + 1
#  }
#  k <- k + 1
#  }

#rect(xleft = 0.5, xright = 29.5, ybottom = 93, ytop = 105, col = "grey95", border = F, xpd = T)

#text(x = seq(3,27,6), y = 98.5, labels = paste(pop.sz/100, "/ ha"))
#abline(v = seq(6,24,6))

#axis(side = 2, las = 2)
#mtext(side = 2, text = expression(paste("deviation from ", italic(N), " (%)")), line = 2.5)

#leg.place <- c(1, 5, 9, 13, 17)
#points(x = leg.place, , y = rep(-112,5), bg = cols, pch = 22, cex = 2.2, xpd = TRUE)
#y.adj <- c(-1,0,0,-1,0)*0.4
#text(x = leg.place, y = rep(-125,5)+y.adj, pch = 19, cex = 1.1, xpd = TRUE, labels = c("Intercept-var.",               #expression(paste(xi(italic(d)))),
#                                                                                       "Scale-var.",                     # expression(paste(sigma(italic(d)))),
#                                                                                       "Conv. DA",                # expression(paste(sigma)),
#                                                                                       "Scale/Intercept",           # expression(paste(xi(italic(d)), " or ", sigma(italic(d)))),
#                                                                                       "Nested doll"))

#text(x = 22, y = -115, labels = "generating model = Nested doll", xpd = T, pos = 4)          ##  expression(paste("generating model = ", xi(italic(d))))
#text(x = 22, y = -125, labels = "transect interval = 200m", xpd = T, pos = 4)

#dev.off()

