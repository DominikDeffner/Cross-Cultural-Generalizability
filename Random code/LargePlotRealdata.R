



library(scales)
#color stuff
require(RColorBrewer)#load package
col.pal <- brewer.pal(8, "Dark2") #create a pallette which you loop over for corresponding values
seqoverall <- seq


graphics.off()
png("AgeCurvesReal.png", res = 900, height = 10, width = 22, units = "cm")



#Age Curves

  par(mfrow = c(1,2),
      mar = c(1,2,1.2,0.5), 
      oma = c(2.5,2,1,0))
  
  ###
  ##
  # Berlin
  ##
  ###
  samp <- extract.samples(m_Berlin)
  samp_pred_p <- apply(samp$pred_p_m,2,quantile,probs=c(0.055,0.5,0.945))
  Lower <- samp_pred_p[1,]
  Median <- samp_pred_p[2,]
  Higher <- samp_pred_p[3,]
  
  Observed_m <- ifelse(c(1:20) %in% unique(d_list_Berlin$age[d_list_Berlin$gender==1]),"Observed Age","Unobserved Age")
  Observed_f <- ifelse(c(1:20) %in% unique(d_list_Berlin$age[d_list_Berlin$gender==2]),"Observed Age","Unobserved Age")
  
  plot(NA, type = "l", xlim = c(0,20), ylim = c(0,1), xlab = "", ylab = "", lwd = 1, lty = 2, main= "Berlin (Germany)")
  points(x= 1:20 - 0.25, Median, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.2)), pch = 16)
  segments(x0 = 1:20-0.25, y0 = Lower, x1=  1:20-0.25, y1 = Higher, lwd = 1, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.2)))
  
  
  
  samp_pred_p <- apply(samp$pred_p_f,2,quantile,probs=c(0.055,0.5,0.945))
  Lower <- samp_pred_p[1,]
  Median <- samp_pred_p[2,]
  Higher <- samp_pred_p[3,]
  
  par(new = TRUE)
  plot(NA, type = "l", xlim = c(0,20), ylim = c(0,1), xlab = "", ylab = "", lwd = 1, lty = 2, main= "")
  points(x= 1:20 + 0.25,Median, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.2)), pch = 17)
  segments(x0 = 1:20+0.25, y0 = Lower, x1=  1:20+0.25, y1 = Higher, lwd = 1, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.2)))
  
  

  
  
  
  
  
  ###
  ##
  # Vanuatu
  ##
  ###
  samp <- extract.samples(m_Vanuatu)
  samp_pred_p <- apply(samp$pred_p_m,2,quantile,probs=c(0.055,0.5,0.945))
  Lower <- samp_pred_p[1,]
  Median <- samp_pred_p[2,]
  Higher <- samp_pred_p[3,]
  
  Observed_m <- ifelse(c(1:20) %in% unique(d_list_Vanuatu$age[d_list_Vanuatu$gender==1]),"Observed Age","Unobserved Age")
  Observed_f <- ifelse(c(1:20) %in% unique(d_list_Vanuatu$age[d_list_Vanuatu$gender==2]),"Observed Age","Unobserved Age")
  
  plot(NA, type = "l", xlim = c(0,20), ylim = c(0,1), xlab = "", ylab = "", lwd = 1, lty = 2, main= "Tanna (Vanuatu)")
  points(x= 1:20 - 0.25, Median, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.2)), pch = 16)
  segments(x0 = 1:20-0.25, y0 = Lower, x1=  1:20-0.25, y1 = Higher, lwd = 1, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.2)))
  
  legend("topleft", c("Male","Female"), col = c(alpha(col.pal[2],alpha = 1), alpha(col.pal[3],alpha = 1)), pch = c(16,17), lwd = 1, lty=1, bty="n", cex = 1.2)
  
  samp_pred_p <- apply(samp$pred_p_f,2,quantile,probs=c(0.055,0.5,0.945))
  Lower <- samp_pred_p[1,]
  Median <- samp_pred_p[2,]
  Higher <- samp_pred_p[3,]
  
  par(new = TRUE)
  plot(NA, type = "l", xlim = c(0,20), ylim = c(0,1), xlab = "", ylab = "", lwd = 1, lty = 2, main= "")
  points(x= 1:20 + 0.25,Median, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.2)), pch = 17)
  segments(x0 = 1:20+0.25, y0 = Lower, x1=  1:20+0.25, y1 = Higher, lwd = 1, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.2)))
  
  
  
  
  
  
  
  
  mtext("Age category", side = 1, outer = TRUE, line = 1.5, cex = 1.2)
  mtext("Probability of choosing prosocial option", side = 2, line = 1,outer = TRUE)
  
  
  
  dev.off()
  
  
  
  


graphics.off()
png("DemostandBERVAT.png", res = 900, height = 19, width = 15, units = "cm")



par(mfrow = c(3,2), 
    mar = c(4,1,2,2), 
    oma = c(0.1,4,0,0.1))




par(mar=pyramid.plot(Pop_Vanuatu[,1],Pop_Vanuatu[,2],top.labels=c("", "Tanna (Vanuatu)",""), ppmar=c(2,1,3,1), xlim = c(10,10),labelcex=1, unit = "",show.values=F, labels = labels1, lxcol = col.pal[2], rxcol = col.pal[3],space = 0.2,gap = 0))
mtext("Population", side = 2, outer = F, line = 3.5, cex = 1.3)
mtext("Age class", side = 2, outer = F, line = 0.2, cex = 1)


par(mar=pyramid.plot(Pop_Berlin[,1],Pop_Berlin[,2],top.labels=c("", "Berlin (Germany)",""), ppmar=c(2,1,3,1), xlim = c(10,10),labelcex=1, unit = "",show.values=F, labels = labels2, lxcol = col.pal[2], rxcol = col.pal[3],space = 0.2,gap = 0))
mtext("Share of population per age class and gender [%]", side = 1,line = 4.5,at = -11, outer = F, cex = 0.9)

legend("topleft", c("Male", "Female"), col = c(col.pal[2], col.pal[3]),cex = 1, lty = 1,lwd = 5, bty = "n" )


par(mar=pyramid.plot(Sample_Vanuatu[,1],Sample_Vanuatu[,2],top.labels=c("", "",""), ppmar=c(2,1,3,1), xlim = c(30,30),labelcex=1, unit = "",show.values=F, labels = labels1, lxcol = col.pal[2], rxcol = col.pal[3],space = 0.2,gap = 0))
mtext("Sample", side = 2, outer = F, line = 3.5, cex = 1.3)
mtext("Age class", side = 2, outer = F, line = 0.2, cex = 1)

par(mar=pyramid.plot(Sample_Berlin[,1],Sample_Berlin[,2],top.labels=c("", "",""), ppmar=c(2,1,3,1), xlim = c(30,30),labelcex=1, unit = "",show.values=F, labels = labels2, lxcol = col.pal[2], rxcol = col.pal[3],space = 0.2,gap = 0))



mtext("Number of individuals per age class and gender", side = 1,line = 4.5,at = -33, outer = F, cex = 0.9)



#Densities 1
s_Berlin <- extract.samples(m_Berlin)
s_Vanuatu <- extract.samples(m_Vanuatu)

dens <- density(s_Vanuatu$p_source)
x1 <- min(which(dens$x >= quantile(s_Vanuatu$p_source, 0)))  
x2 <- max(which(dens$x <  quantile(s_Vanuatu$p_source, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,15), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[1],alpha = 0.5), border = NA))
par(new=TRUE)
dens <- density(s_Vanuatu$p_pop)
x1 <- min(which(dens$x >= quantile(s_Vanuatu$p_pop, 0)))  
x2 <- max(which(dens$x <  quantile(s_Vanuatu$p_pop, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,15), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[6],alpha = 0.5), border = NA))
mtext("Density", side = 2,line = 3, outer = F, cex = 1)

par(new=TRUE)
dens <- density(s_Vanuatu$p_other)
x1 <- min(which(dens$x >= quantile(s_Vanuatu$p_other, 0)))  
x2 <- max(which(dens$x <  quantile(s_Vanuatu$p_other, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,15), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[8],alpha = 0.5), border = NA))

legend("topleft", c("Empirical estimate","Poststratified to sample population", "Poststratified to other population"), col = c(alpha(col.pal[1],alpha = 0.5),alpha(col.pal[6],alpha = 0.5),alpha(col.pal[8],alpha = 0.5)), lwd = 6, bty="n", cex = 1)




dens <- density(s_Berlin$p_source)
x1 <- min(which(dens$x >= quantile(s_Berlin$p_source, 0)))  
x2 <- max(which(dens$x <  quantile(s_Berlin$p_source, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,15), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[1],alpha = 0.5), border = NA))
par(new=TRUE)
dens <- density(s_Berlin$p_pop)
x1 <- min(which(dens$x >= quantile(s_Berlin$p_pop, 0)))  
x2 <- max(which(dens$x <  quantile(s_Berlin$p_pop, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,15), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[6],alpha = 0.5), border = NA))


par(new=TRUE)
dens <- density(s_Berlin$p_other)
x1 <- min(which(dens$x >= quantile(s_Berlin$p_other, 0)))  
x2 <- max(which(dens$x <  quantile(s_Berlin$p_other, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,15), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[8],alpha = 0.5), border = NA))



mtext("Probability of choosing prosocial option", side = 1,line = 2.9,at=-0.1, outer = F, cex = 0.9)



dev.off()

