



library(scales)
#color stuff
require(RColorBrewer)#load package
col.pal <- brewer.pal(8, "Dark2") #create a pallette which you loop over for corresponding values
seqoverall <- seq




graphics.off()
png("DemostandLarge.png", res = 900, height = 20, width = 24, units = "cm")



par(mfrow = c(3,4), 
    mar = c(3,1,3,2), 
    oma = c(1,4,2,0.1))
labels1 <- matrix("", length(Age_range), 2)
labels1[seq(5,75,10),1] <- seq(5,75,10)
labels2 <- matrix("", length(Age_range), 2)

# Population 1
par(mar=pyramid.plot(D1[1,1, ]*100,D1[1,2, ]*100,top.labels=c("", "Population I",""),ppmar=c(2,1,3,1), xlim = c(5,5),labelcex=1.2, unit = "",show.values=F, labels = labels1, lxcol = col.pal[2], rxcol = col.pal[3],space = 0,gap = 0))
labels <- matrix("", length(Age_range), 2)
legend("topright", c("Male", "Female"), col = c(col.pal[2], col.pal[3]),cex = 1, lty = 1,lwd = 5, bty = "n" )

mtext("Population", side = 2, outer = F, line = 3.5, cex = 1.3)
mtext("Age class", side = 2, outer = F, line = 0.2, cex = 1)
par(mar=pyramid.plot(D1[2,1, ]*100,D1[2,2, ]*100,top.labels=c("", "Population II",""),ppmar=c(2,1,3,1), xlim = c(3,3),labelcex=1.2, unit = "",show.values=F, labels = labels2,  lxcol = col.pal[2], rxcol = col.pal[3], space = 0,gap = 0))
mtext("Share of population per age class and gender [%]", side = 1,line = 4,at = -3, outer = F, cex = 0.9)
mtext("Population Differences", side = 3,line = 3,at = -3, outer = F, cex = 2)

# Population 2
par(mar=pyramid.plot(D2[1,1, ]*100,D2[1,2, ]*100,top.labels=c("", "Population I",""),ppmar=c(2,1,3,1), xlim = c(5,5),labelcex=1.2, unit = "",show.values=F, labels = labels2, lxcol = col.pal[2], rxcol = col.pal[3],space = 0,gap = 0))
labels <- matrix("", length(Age_range), 2)
par(mar=pyramid.plot(D2[2,1, ]*100,D2[2,2, ]*100,top.labels=c("", "Population II",""),ppmar=c(2,1,3,1), xlim = c(5,5),labelcex=1.2, unit = "",show.values=F, labels = labels2,  lxcol = col.pal[2], rxcol = col.pal[3], space = 0,gap = 0))

mtext("Share of population per age class and gender [%]", side = 1,line = 4,at = -5, outer = F, cex = 0.9)
mtext("Sampling Differences", side = 3,line = 3,at = -5, outer = F, cex = 2)



#Samples 1

par(mar=pyramid.plot(a11[1,],a11[2,],top.labels=c("", "",""),ppmar=c(2,1,3,1), xlim = c(max(c(a11,a12)+3),max(c(a11,a12))+3),labelcex=1.2, unit = "",show.values=F, labels = labels1,  lxcol = col.pal[2], rxcol = col.pal[3],space = 0,gap = 0))
mtext("Sample", side = 2, outer = F, line = 3.5, cex = 1.3)
mtext("Age class", side = 2, outer = F, line = 0.2, cex = 1)
par(mar=pyramid.plot(a12[1,],a12[2,],top.labels=c("", "",""),ppmar=c(2,1,3,1), xlim = c(max(c(a11,a12)+3),max(c(a11,a12))+3),labelcex=1.2, unit = "",show.values=F, labels = labels2,  lxcol = col.pal[2], rxcol = col.pal[3],space = 0,gap = 0))

mtext("Number of individuals per age class and gender", side = 1,line = 4,at = -22, outer = F, cex = 0.9)

#Samples 2

par(mar=pyramid.plot(a21[1,],a21[2,],top.labels=c("", "",""),ppmar=c(2,1,3,1), xlim = c(max(c(a21,a22)+3),max(c(a21,a22))+3),labelcex=1.2, unit = "",show.values=F, labels = labels2,  lxcol = col.pal[2], rxcol = col.pal[3],space = 0,gap = 0))
par(mar=pyramid.plot(a22[1,],a22[2,],top.labels=c("", "",""),ppmar=c(2,1,3,1), xlim = c(max(c(a21,a22)+3),max(c(a21,a22))+3),labelcex=1.2, unit = "",show.values=F, labels = labels2,  lxcol = col.pal[2], rxcol = col.pal[3],space = 0,gap = 0))

mtext("Number of individuals per age class and gender", side = 1,line = 4,at = -22, outer = F, cex = 0.9)






#Densities 1


dens <- density(m11samp$phi)
x1 <- min(which(dens$x >= quantile(m11samp$phi, 0)))  
x2 <- max(which(dens$x <  quantile(m11samp$phi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[4],alpha = 0.9), border = NA))
abline(v = phi1[1], lty = 2, lwd = 2)
par(new=TRUE)
dens <- density(m11samp$psi)
x1 <- min(which(dens$x >= quantile(m11samp$psi, 0)))  
x2 <- max(which(dens$x <  quantile(m11samp$psi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[4],alpha = 0.3), border = NA))
mtext("Density", side = 2,line = 3, outer = F, cex = 1)



dens <- density(m12samp$phi)
x1 <- min(which(dens$x >= quantile(m12samp$phi, 0)))  
x2 <- max(which(dens$x <  quantile(m12samp$phi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[5],alpha = 0.9), border = NA))
abline(v = phi1[2], lty = 2, lwd = 2)

par(new=TRUE)

dens <- density(m12samp$psi)
x1 <- min(which(dens$x >= quantile(m12samp$psi, 0)))
x2 <- max(which(dens$x <  quantile(m12samp$psi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[5],alpha = 0.3), border = NA))

mtext("Probability of choosing prosocial option", side = 1,line = 2.9,at=-0.1, outer = F, cex = 0.9)
legend("topleft", title = "Population", c("I","I (poststratified)", "II","II (poststratified)"), col = c(alpha(col.pal[4],alpha = 0.9),alpha(col.pal[4],alpha = 0.3),alpha(col.pal[5],alpha = 0.9),alpha(col.pal[5],alpha = 0.3)), lwd = 6, bty="n", cex = 1)







#Densities 2


dens <- density(m21samp$phi)
x1 <- min(which(dens$x >= quantile(m21samp$phi, 0)))  
x2 <- max(which(dens$x <  quantile(m21samp$phi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", yaxt = "n",ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[4],alpha = 0.9), border = NA))
abline(v = phi2[1], lty = 2, lwd = 2)
par(new=TRUE)
dens <- density(m21samp$psi)
x1 <- min(which(dens$x >= quantile(m21samp$psi, 0)))  
x2 <- max(which(dens$x <  quantile(m21samp$psi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[4],alpha = 0.3), border = NA))



dens <- density(m22samp$phi)
x1 <- min(which(dens$x >= quantile(m22samp$phi, 0)))  
x2 <- max(which(dens$x <  quantile(m22samp$phi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[5],alpha = 0.9), border = NA))
abline(v = phi2[2], lty = 2, lwd = 2)

par(new=TRUE)

dens <- density(m22samp$psi)
x1 <- min(which(dens$x >= quantile(m22samp$psi, 0)))
x2 <- max(which(dens$x <  quantile(m22samp$psi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[5],alpha = 0.3), border = NA))

mtext("Probability of choosing prosocial option", side = 1,line = 2.9,at=-0.1, outer = F, cex = 0.9)

dev.off()

