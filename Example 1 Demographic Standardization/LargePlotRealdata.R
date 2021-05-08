



library(scales)
#color stuff
require(RColorBrewer)#load package
col.pal <- brewer.pal(8, "Dark2") #create a pallette which you loop over for corresponding values
seqoverall <- seq


graphics.off()
png("AgeCurves.png", res = 900, height = 16, width = 24, units = "cm")




#Age Curves

  par(mfrow = c(2,2),
      mar = c(1,2,1.2,0.5), 
      oma = c(2.2,5.1,1,0))
  
  ###
  ##
  # Population Differences
  ##
  ###
  samp <- extract.samples(m_Berlin)
  samp_pred_p <- apply(samp$pred_p_m,2,quantile,probs=c(0.05,0.5,0.95))
  Lower <- samp_pred_p[1,]
  Median <- samp_pred_p[2,]
  Higher <- samp_pred_p[3,]
  
  Observed_m <- ifelse(c(1:20) %in% unique(d_list_Berlin$age[d_list_Berlin$gender==1]),"Observed Age","Unobserved Age")
  Observed_f <- ifelse(c(1:20) %in% unique(d_list_Berlin$age[d_list_Berlin$gender==2]),"Observed Age","Unobserved Age")
  
  plot(NA, type = "l", xlim = c(0,20), ylim = c(0,1), xlab = "", ylab = "", lwd = 2, lty = 2, main= "Population I")
  points(x= 1:20 - 0.3, Median, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.3)), pch = 16)
  segments(x0 = 1:20-0.3, y0 = Lower, x1=  1:20-0.3, y1 = Higher, lwd = 2, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.3)))
  
  legend("topleft", c("Male","Female"), col = c(alpha(col.pal[2],alpha = 1), alpha(col.pal[3],alpha = 1)), pch = c(16,17), lwd = 2, lty=1, bty="n", cex = 1.2)
  
  
  samp_pred_p <- apply(samp$pred_p_f,2,quantile,probs=c(0.05,0.5,0.95))
  Lower <- samp_pred_p[1,]
  Median <- samp_pred_p[2,]
  Higher <- samp_pred_p[3,]
  
  par(new = TRUE)
  plot(NA, type = "l", xlim = c(0,20), ylim = c(0,1), xlab = "", ylab = "", lwd = 2, lty = 2, main= "")
  points(x= 1:20 + 0.3,Median, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.3)), pch = 17)
  segments(x0 = 1:20+0.3, y0 = Lower, x1=  1:20+0.3, y1 = Higher, lwd = 2, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.3)))
  
  
  mtext("Population Differences", side = 2, outer = FALSE, line = 6, cex = 1.2)
  
  
  
  
  
  
  
  samp <- extract.samples(m_popdiff2)
  samp_pred_p <- apply(samp$pred_p_m,2,quantile,probs=c(0.05,0.5,0.95))
  Lower <- samp_pred_p[1,]
  Median <- samp_pred_p[2,]
  Higher <- samp_pred_p[3,]
  
  Observed_m <- ifelse(c(1:d2_popdiff$MA) %in% unique(d2_popdiff$age[d2_popdiff$gender==1]),"Observed Age","Unobserved Age")
  Observed_f <- ifelse(c(1:d2_popdiff$MA) %in% unique(d2_popdiff$age[d2_popdiff$gender==2]),"Observed Age","Unobserved Age")
  
  plot(NA, type = "l", xlim = c(0,max(Age_range)), ylim = c(0,1), xlab = "", ylab = "", lwd = 2, lty = 2, main= "Population II")
  points(Median, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.3)), pch = 16)
  segments(x0 = 1:d2_popdiff$MA, y0 = Lower, x1=  1:d2_popdiff$MA, y1 = Higher, lwd = 2, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.3)))
  

  samp_pred_p <- apply(samp$pred_p_f,2,quantile,probs=c(0.05,0.5,0.95))
  Lower <- samp_pred_p[1,]
  Median <- samp_pred_p[2,]
  Higher <- samp_pred_p[3,]
  
  par(new = TRUE)
  plot(NA, type = "l", xlim = c(0,max(Age_range)), ylim = c(0,1), xlab = "", ylab = "", lwd = 2, lty = 2, main= "")
  points(Median, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.3)), pch = 17)
  segments(x0 = 1:d2_popdiff$MA, y0 = Lower, x1=  1:d2_popdiff$MA, y1 = Higher, lwd = 2, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.3)))
  
  
  
  
  
  ###
  ##
  # Sampling Differences
  ##
  ###
  
  samp <- extract.samples(m_samplediff1)
  samp_pred_p <- apply(samp$pred_p_m,2,quantile,probs=c(0.05,0.5,0.95))
  Lower <- samp_pred_p[1,]
  Median <- samp_pred_p[2,]
  Higher <- samp_pred_p[3,]
  
  Observed_m <- ifelse(c(1:d1_samplediff$MA) %in% unique(d1_samplediff$age[d1_samplediff$gender==1]),"Observed Age","Unobserved Age")
  Observed_f <- ifelse(c(1:d1_samplediff$MA) %in% unique(d1_samplediff$age[d1_samplediff$gender==2]),"Observed Age","Unobserved Age")
  
  plot(NA, type = "l", xlim = c(0,max(Age_range)), ylim = c(0,1), xlab = "", ylab = "", lwd = 2, lty = 2, main= "")
  points(Median, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.3)), pch = 16)
  segments(x0 = 1:d1_samplediff$MA, y0 = Lower, x1=  1:d1_samplediff$MA, y1 = Higher, lwd = 2, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.3)))
  
  
  mtext("Sampling Differences", side = 2, outer = FALSE,line = 6, cex = 1.2)
  

  
  
  samp_pred_p <- apply(samp$pred_p_f,2,quantile,probs=c(0.05,0.5,0.95))
  Lower <- samp_pred_p[1,]
  Median <- samp_pred_p[2,]
  Higher <- samp_pred_p[3,]
  
  par(new = TRUE)
  plot(NA, type = "l", xlim = c(0,max(Age_range)), ylim = c(0,1), xlab = "", ylab = "", lwd = 2, lty = 2, main= "")
  points(Median, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.3)), pch = 17)
  segments(x0 = 1:d1_samplediff$MA, y0 = Lower, x1=  1:d1_samplediff$MA, y1 = Higher, lwd = 2, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.3)))
  
  
  
  
  
  
  
  
  
  samp <- extract.samples(m_samplediff2)
  samp_pred_p <- apply(samp$pred_p_m,2,quantile,probs=c(0.05,0.5,0.95))
  Lower <- samp_pred_p[1,]
  Median <- samp_pred_p[2,]
  Higher <- samp_pred_p[3,]
  
  Observed_m <- ifelse(c(1:d2_samplediff$MA) %in% unique(d2_samplediff$age[d2_samplediff$gender==1]),"Observed Age","Unobserved Age")
  Observed_f <- ifelse(c(1:d2_samplediff$MA) %in% unique(d2_samplediff$age[d2_samplediff$gender==2]),"Observed Age","Unobserved Age")
  
  plot(NA, type = "l", xlim = c(0,max(Age_range)), ylim = c(0,1), xlab = "", ylab = "", lwd = 2, lty = 2, main= "")
  points(Median, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.3)), pch = 16)
  segments(x0 = 1:d2_samplediff$MA, y0 = Lower, x1=  1:d2_samplediff$MA, y1 = Higher, lwd = 2, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.3)))
  
  
  
  
  
  samp_pred_p <- apply(samp$pred_p_f,2,quantile,probs=c(0.05,0.5,0.95))
  Lower <- samp_pred_p[1,]
  Median <- samp_pred_p[2,]
  Higher <- samp_pred_p[3,]
  
  par(new = TRUE)
  plot(NA, type = "l", xlim = c(0,max(Age_range)), ylim = c(0,1), xlab = "", ylab = "", lwd = 2, lty = 2, main= "")
  points(Median, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.3)), pch = 17)
  segments(x0 = 1:d2_samplediff$MA, y0 = Lower, x1=  1:d2_samplediff$MA, y1 = Higher, lwd = 2, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.3)))
  
  
  
  
  
  
  
  
  mtext("Age", side = 1, outer = TRUE, line = 1, cex = 1.2)
  mtext("Probability of choosing prosocial option", side = 2, line = 1.2,outer = TRUE)
  
  
  
  dev.off()
  
  
  
  
{  

graphics.off()
png("DemostandLarge2.png", res = 900, height = 25, width = 25, units = "cm")



par(mfrow = c(3,4), 
    mar = c(3,1,3,2), 
    oma = c(1,4,2,0.1))
labels1 <- matrix("", length(Age_range), 2)
labels1[seq(5,75,10),1] <- seq(5,75,10)
labels2 <- matrix("", length(Age_range), 2)

# Population 1
par(mar=pyramid.plot(D_popdiff[1,1, ]*100,D_popdiff[1,2, ]*100,top.labels=c("", "Population I",""),ppmar=c(2,1,3,1), xlim = c(5,5),labelcex=1.2, unit = "",show.values=F, labels = labels1, lxcol = col.pal[2], rxcol = col.pal[3],space = 0,gap = 0))
labels <- matrix("", length(Age_range), 2)
legend("topright", c("Male", "Female"), col = c(col.pal[2], col.pal[3]),cex = 1, lty = 1,lwd = 5, bty = "n" )

mtext("Population", side = 2, outer = F, line = 3.5, cex = 1.3)
mtext("Age class", side = 2, outer = F, line = 0.2, cex = 1)
par(mar=pyramid.plot(D_popdiff[2,1, ]*100,D_popdiff[2,2, ]*100,top.labels=c("", "Population II",""),ppmar=c(2,1,3,1), xlim = c(3,3),labelcex=1.2, unit = "",show.values=F, labels = labels2,  lxcol = col.pal[2], rxcol = col.pal[3], space = 0,gap = 0))
mtext("Share of population per age class and gender [%]", side = 1,line = 4,at = -3, outer = F, cex = 0.9)
mtext("Population Differences", side = 3,line = 3,at = -3, outer = F, cex = 2)

# Population 2
par(mar=pyramid.plot(D_samplediff[1,1, ]*100,D_samplediff[1,2, ]*100,top.labels=c("", "Population I",""),ppmar=c(2,1,3,1), xlim = c(5,5),labelcex=1.2, unit = "",show.values=F, labels = labels2, lxcol = col.pal[2], rxcol = col.pal[3],space = 0,gap = 0))
labels <- matrix("", length(Age_range), 2)
par(mar=pyramid.plot(D_samplediff[2,1, ]*100,D_samplediff[2,2, ]*100,top.labels=c("", "Population II",""),ppmar=c(2,1,3,1), xlim = c(5,5),labelcex=1.2, unit = "",show.values=F, labels = labels2,  lxcol = col.pal[2], rxcol = col.pal[3], space = 0,gap = 0))

mtext("Share of population per age class and gender [%]", side = 1,line = 4,at = -5, outer = F, cex = 0.9)
mtext("Sampling Differences", side = 3,line = 3,at = -5, outer = F, cex = 2)



#Samples 1

par(mar=pyramid.plot(SampleD_popdiff[1,1,],SampleD_popdiff[1,2,],top.labels=c("", "",""),ppmar=c(2,1,3,1), xlim = c(max(SampleD_popdiff[1,,]+3),max(SampleD_popdiff[1,,]+3)),labelcex=1.2, unit = "",show.values=F, labels = labels1,  lxcol = col.pal[2], rxcol = col.pal[3],space = 0,gap = 0))
mtext("Sample", side = 2, outer = F, line = 3.5, cex = 1.3)
mtext("Age class", side = 2, outer = F, line = 0.2, cex = 1)
par(mar=pyramid.plot(SampleD_popdiff[2,1,],SampleD_popdiff[2,2,],top.labels=c("", "",""),ppmar=c(2,1,3,1), xlim = c(max(SampleD_popdiff[1,,]+3),max(SampleD_popdiff[1,,]+3)),labelcex=1.2, unit = "",show.values=F, labels = labels2,  lxcol = col.pal[2], rxcol = col.pal[3],space = 0,gap = 0))

mtext("Number of individuals per age class and gender", side = 1,line = 4,at = -22, outer = F, cex = 0.9)

#Samples 2

par(mar=pyramid.plot(SampleD_samplediff[1,1,],SampleD_samplediff[1,2,],top.labels=c("", "",""),ppmar=c(2,1,3,1), xlim = c(max(SampleD_samplediff[2,,]+3),max(SampleD_samplediff[2,,]+3)),labelcex=1.2, unit = "",show.values=F, labels = labels2,  lxcol = col.pal[2], rxcol = col.pal[3],space = 0,gap = 0))
par(mar=pyramid.plot(SampleD_samplediff[2,1,],SampleD_samplediff[2,2,],top.labels=c("", "",""),ppmar=c(2,1,3,1), xlim = c(max(SampleD_samplediff[2,,]+3),max(SampleD_samplediff[2,,]+3)),labelcex=1.2, unit = "",show.values=F, labels = labels2,  lxcol = col.pal[2], rxcol = col.pal[3],space = 0,gap = 0))

mtext("Number of individuals per age class and gender", side = 1,line = 4,at = -22, outer = F, cex = 0.9)


#Densities 1

dens <- density(s_popdiff1$phi)
x1 <- min(which(dens$x >= quantile(s_popdiff1$phi, 0)))  
x2 <- max(which(dens$x <  quantile(s_popdiff1$phi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[4],alpha = 0.9), border = NA))
abline(v = phi_popdiff[1], lty = 2, lwd = 2)
par(new=TRUE)
dens <- density(s_popdiff1$psi)
x1 <- min(which(dens$x >= quantile(s_popdiff1$psi, 0)))  
x2 <- max(which(dens$x <  quantile(s_popdiff1$psi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[4],alpha = 0.3), border = NA))
mtext("Density", side = 2,line = 3, outer = F, cex = 1)



dens <- density(s_popdiff2$phi)
x1 <- min(which(dens$x >= quantile(s_popdiff2$phi, 0)))  
x2 <- max(which(dens$x <  quantile(s_popdiff2$phi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[5],alpha = 0.9), border = NA))
abline(v = phi_popdiff[2], lty = 2, lwd = 2)

par(new=TRUE)

dens <- density(s_popdiff2$psi)
x1 <- min(which(dens$x >= quantile(s_popdiff2$psi, 0)))
x2 <- max(which(dens$x <  quantile(s_popdiff2$psi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[5],alpha = 0.3), border = NA))

mtext("Probability of choosing prosocial option", side = 1,line = 2.9,at=-0.1, outer = F, cex = 0.9)
legend("topleft", title = "Population", c("I","I (poststratified)", "II","II (poststratified)"), col = c(alpha(col.pal[4],alpha = 0.9),alpha(col.pal[4],alpha = 0.3),alpha(col.pal[5],alpha = 0.9),alpha(col.pal[5],alpha = 0.3)), lwd = 6, bty="n", cex = 1)







#Densities 2


dens <- density(s_samplediff1$phi)
x1 <- min(which(dens$x >= quantile(s_samplediff1$phi, 0)))  
x2 <- max(which(dens$x <  quantile(s_samplediff1$phi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", yaxt = "n",ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[4],alpha = 0.9), border = NA))
abline(v = phi_samplediff[1], lty = 2, lwd = 2)
par(new=TRUE)
dens <- density(s_samplediff1$psi)
x1 <- min(which(dens$x >= quantile(s_samplediff1$psi, 0)))  
x2 <- max(which(dens$x <  quantile(s_samplediff1$psi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[4],alpha = 0.3), border = NA))



dens <- density(s_samplediff2$phi)
x1 <- min(which(dens$x >= quantile(s_samplediff2$phi, 0)))  
x2 <- max(which(dens$x <  quantile(s_samplediff2$phi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[5],alpha = 0.9), border = NA))
abline(v = phi_samplediff[2], lty = 2, lwd = 2)

par(new=TRUE)

dens <- density(s_samplediff2$psi)
x1 <- min(which(dens$x >= quantile(s_samplediff2$psi, 0)))
x2 <- max(which(dens$x <  quantile(s_samplediff2$psi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[5],alpha = 0.3), border = NA))

mtext("Probability of choosing prosocial option", side = 1,line = 2.9,at=-0.1, outer = F, cex = 0.9)

dev.off()

}
