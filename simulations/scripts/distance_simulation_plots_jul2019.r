##  Plot output of function sim().
##
##
##
##  09 jul 2019
################################################################

diam.cols <- function (diam.v) {
   ###  Palette of hues:
 cr <- colorRamp(colors = c("orange", "gray", "blue"), space = "Lab")
 cp <- colorRampPalette(colors = c("orange", "gray", "blue"), space = "Lab")
 diam.sc <- diam.v/max(diam.v)
 diam.cols <-  rgb(red = cr(diam.sc)[,1], green = cr(diam.sc)[,2], blue = cr(diam.sc)[,3], maxColorValue = 255)
 
 mget(c("cp", "diam.cols")) 
}



##  "s" is an object produced by burrSim(), and of class "burrSim".
map.burrSim <- function (s) {

 stopifnot(class(s) == "burrSim")
 
 col.lst <- diam.cols(s$all.burr$diam)

 ###  Map:
# x11(width = 9, height = 7)
 png("/home/nuse/Projects/GT_stuff/distance/manuscript/latex/images/simul_map.png", width = 8, height = 6, units = "in", res = 400)
 layout(matrix(c(1, 2), ncol = 2), widths = c(1, 0.2))
 par(mar = c(5, 4, 4, 2) + 0.1)
 plot(s$all.burr$loc.x, s$all.burr$loc.y, type = "n", xlim = c(0,s$survey$edge), ylim = c(0,s$survey$edge), xlab = "x (m)", ylab = "y (m)", xaxs = "i", yaxs = "i")

###  Plot the transect lines:
 for (i in 1:s$survey$n.lines) { 
  curve(trans.y(x, theta = s$survey$theta, spacing = s$survey$intrvl, y.int =  s$survey$y.int.min, trans.num = i), from = 0, to = s$survey$edge, add = T)
 }
 
 ###  add points:
 points(s$all.burr$loc.x, s$all.burr$loc.y, col = col.lst$diam.cols, pch = 16, cex = 1.5)

 ###  Color Legend!
 levels <- pretty(s$all.burr$diam/10, 13)  ## change to centimeters
 par(mar = c(5,0,4,4)+0.1)
 plot.new()
 plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i",
             yaxs = "i")
 rect(0, levels[-length(levels)], 1, levels[-1], col = col.lst$cp(length(levels)))
 tck <- axTicks(side = 4)
 axis(4, las = 1, at = tck, labels = paste(tck, c("cm", rep("", length(tck)-1))))

 dev.off()
}    ### END map.burrSim()



########################################################
########################################################
########################################################

 ##  "s" is an object produced by burrSim(), and of class "burrSim".
plot.p.burrSim <- function (s) {

 ###  Make a picture of p varying with distance and diameter:
# x11(width = 9, height = 7)
 png("/home/nuse/Projects/GT_stuff/distance/manuscript/latex/images/simul_plot.png", width = 8, height = 6, units = "in", res = 400)

 col.lst <- diam.cols(s$all.burr$diam)
 plot(p ~ dist, col = col.lst$diam.cols, data = s$all.burr, pch = 16, cex = 1.5, xlim = c(0,0.5*s$survey$intrvl), ylim = c(0,1), xlab = "Distance from transect (m)", ylab = "Pr(detection)")
 this.x <- 0.7*par("usr")[2]
 text(y = 0.9, x = this.x, labels = bquote(alpha[xi] == .(s$par$xi.50)), cex = 1.5, pos = 2, family = "mono")    ## trick from examples on the plotmath help page....
 text(y = 0.8, x = this.x, labels = bquote(italic(W)[xi] == .(s$par$xi.brk/10)), cex = 1.5, pos = 2, family = "mono")    ## convert to cm
 text(y = 0.7, x = this.x, labels = bquote(sigma == .(s$par$sigma)), cex = 1.5, pos = 2, family = "mono")    


 ###  Show detections on the p-dist plot:
 with(s$all.burr[s$all.burr$found == 1,], {points(x = dist, y = p, col = "black", pch = 3, cex = 3, lwd = 1.5); points(x = dist, y = p, col = "lightgray", pch = 1, cex = 1.5)})

 dev.off()
}   ###  END plot.p.burrSim() 

#plot.p.burrSim(s)

