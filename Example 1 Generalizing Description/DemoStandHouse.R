
#Generalizing description: Cross-cultural comparisons and demographic standardization
#Real data example based on House et al., 2020
#The script prepares the data, runs the multilevel regression with poststratification analysis and creates plot in Fig.3 in the manuscript

#Load some packages and set working directory to main folder
library(readr)
library(rethinking)
library(plotrix)
library(scales)
require(RColorBrewer)
setwd("~/GitHub/Cross-Cultural-Generalizability")

#Population distributions (age and gender) of Vanuatu and Berlin
Pop_Berlin <- read.csv("data/Berlin-2020.csv")
Pop_Vanuatu <- read.csv("data/Vanuatu-2019.csv")
Pop_Vanuatu <- Pop_Vanuatu[-nrow(Pop_Vanuatu),] #Delete empty cells (100+)
Pop_Berlin <- as.matrix(Pop_Berlin[,c(2,3)])
Pop_Vanuatu <- as.matrix(Pop_Vanuatu[,c(2,3)])

#Change to proportion of different cells in populations
Pop_Berlin_raw <- Pop_Berlin
Pop_Berlin <- (Pop_Berlin/sum(Pop_Berlin))*100

Pop_Vanuatu_raw <- Pop_Vanuatu
Pop_Vanuatu <- (Pop_Vanuatu/sum(Pop_Vanuatu))*100


#Load data from adults and children from House et al.,2020
data_adult <- read.csv("data/Model_1a_1b_1c_data.csv")
data_adult <- data_adult[,c("SUBJECT_ID","GENDER_1female","fieldid", "AGE_in_years","T1_ad_choice_1yes")]
data_adult$choice <- data_adult$T1_ad_choice_1yes
data_adult$T1_ad_choice_1yes <- NULL

data_children <- read.csv("data/Model_4a_4b_4c_4d_data.csv")
data_children <- data_children[,c("SUBJECT_ID","GENDER_1female","fieldid", "AGE_in_years","T1_choice_1yes")]
data_children$choice <- data_children$T1_choice_1yes
data_children$T1_choice_1yes <- NULL


#Combine both datasets
data_comb <- rbind(data_children, data_adult)
data_comb <- data_comb[with(data_comb, order(data_comb$fieldid, data_comb$AGE_in_years)),]
data_comb$gender <- data_comb$GENDER_1female +1
MA <- max(round(data_comb$AGE_in_years))

#Extract data from Berlin and Tanna, Vanuatu
d_Berlin <- data_comb[which(data_comb$fieldid == 1),]
d_Vanuatu <- data_comb[which(data_comb$fieldid == 7),]

#Categorize data into 20 age categories spanning 5 years each
d_Berlin$AGE_binned <- sapply(1:nrow(d_Berlin), function(i) which.max( 1/(seq(5, 100, 5) - d_Berlin$AGE_in_years[i]) ) )
d_Vanuatu$AGE_binned <- sapply(1:nrow(d_Vanuatu), function(i) which.max( 1/(seq(5, 100, 5) - d_Vanuatu$AGE_in_years[i]) ) )


#Sample distributions (age and gender) of Vanuatu and Berlin
Sample_Berlin <- matrix(0, 20, 2)
Sample_Vanuatu <- matrix(0, 20, 2)

for (i in 1:20) {
  for (g in 1:2) {
    Sample_Berlin[i,g] <- length(which(d_Berlin$AGE_binned==i & d_Berlin$gender == g))
    Sample_Vanuatu[i,g] <- length(which(d_Vanuatu$AGE_binned==i & d_Vanuatu$gender == g))
  }
}



#Prepare data lists as input for Stan 

d_list_Berlin <- list(N = nrow(d_Berlin),              #Sample size
                      MA = 20,                         #Number of age categories
                      age = d_Berlin$AGE_binned,       #Binned ages
                      outcome = d_Berlin$choice,       #Choices in dictator game
                      gender = d_Berlin$gender,        #Gender
                      P_Pop = t(Pop_Berlin_raw),       #Population demography of Berlin
                      P_other = t(Pop_Vanuatu_raw))    #Population demography of Vanuatu


d_list_Vanuatu <- list(N = nrow(d_Vanuatu),            #Sample size
                       MA = 20,                         #Number of age categories  
                       age = d_Vanuatu$AGE_binned,      #Binned ages
                       outcome = d_Vanuatu$choice,      #Choices in dictator game
                       gender = d_Vanuatu$gender,       #Gender
                       P_Pop = t(Pop_Vanuatu_raw),      #Population demography of Vanuatu
                       P_other = t(Pop_Berlin_raw))     #Population demography of Berlin



# The next section invokes the stan model code in the example 1 folder

# The first model ("model_empirical") is a simple Bernoulli model that provides unbiased or "empirical" estimate in each population. 

# The second model ("model_MRpoststratification") uses Gaussian processes to compute age and gender specific estimates and poststratifies to population 
# from which sample was taken and to the other population



#Pass code to Rstan and run models
m_Berlin  <- stan( file  = "Example 1 Generalizing Description/model_MRpoststratification.stan" , data=d_list_Berlin  ,iter = 5000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.99, max_treedepth = 13))  
m_Vanuatu <- stan( file  = "Example 1 Generalizing Description/model_MRpoststratification.stan" , data=d_list_Vanuatu ,iter = 5000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.99, max_treedepth = 13))  

m_Berlin_basic  <- stan( file  = "Example 1 Generalizing Description/model_empirical.stan" , data=d_list_Berlin  ,iter = 5000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.99, max_treedepth = 13))  
m_Vanuatu_basic <- stan( file  = "Example 1 Generalizing Description/model_empirical.stan" , data=d_list_Vanuatu ,iter = 5000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.99, max_treedepth = 13))  


######
####
###
##
# Plotting Code
##
###
####
####


###
##
#Gender and age-specific estimates (Fig.S4 in the appendix)
##
###

#graphics.off()
#png("AgeCurvesReal.png", res = 900, height = 10, width = 22, units = "cm")

col.pal <- brewer.pal(8, "Dark2") 
seqoverall <- seq

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

#dev.off()



###
##
# Demographic standardization comparing Tanna, Vanuatu, and Berlin, Germany (Fig.3)
##
###


#graphics.off()
#png("DemostandBERVAT.png", res = 900, height = 19, width = 15, units = "cm")

par(mfrow = c(3,2), 
    mar = c(4,1,2,2), 
    oma = c(0.1,4,0,0.1))

# Create pyramid plot labels
labels1 <- matrix("", 20, 2)
labels1[,1] <- c("1-5", "","11-15","","21-25","","31-35","","41-45","",
                 "51-55","","61-65","","71-75","","81-85","","91-95","")
labels2 <- matrix("", 20, 2)

#Population Pyramids
par(mar=pyramid.plot(Pop_Vanuatu[,1],Pop_Vanuatu[,2],top.labels=c("", "Tanna (Vanuatu)",""), ppmar=c(2,1,3,1),raxlab = seq(0,10,2),laxlab = seq(0,10,2), xlim = c(10,10),labelcex=1, unit = "",show.values=F, labels = labels1, lxcol = col.pal[2], rxcol = col.pal[3],space = 0.2,gap = 0))
mtext("Population", side = 2, outer = F, line = 3.5, cex = 1.3)
mtext("Age class", side = 2, outer = F, line = 0.2, cex = 1)
par(mar=pyramid.plot(Pop_Berlin[,1],Pop_Berlin[,2],top.labels=c("", "Berlin (Germany)",""), ppmar=c(2,1,3,1),raxlab = seq(0,10,2),laxlab = seq(0,10,2), xlim = c(10,10),labelcex=1, unit = "",show.values=F, labels = labels2, lxcol = col.pal[2], rxcol = col.pal[3],space = 0.2,gap = 0))
mtext("Share of population per age class and gender [%]", side = 1,line = 4.5,at = -11, outer = F, cex = 0.9)
legend("topleft", c("Male", "Female"), col = c(col.pal[2], col.pal[3]),cex = 1.1, lty = 1,lwd = 8, bty = "n" )

#Sample Pyramids
par(mar=pyramid.plot(Sample_Vanuatu[,1],Sample_Vanuatu[,2],top.labels=c("", "",""), ppmar=c(2,1,3,1),raxlab = seq(0,30,10),laxlab = seq(0,30,10), xlim = c(30,30),labelcex=1, unit = "",show.values=F, labels = labels1, lxcol = col.pal[2], rxcol = col.pal[3],space = 0.2,gap = 0))
mtext("Sample", side = 2, outer = F, line = 3.5, cex = 1.3)
mtext("Age class", side = 2, outer = F, line = 0.2, cex = 1)
par(mar=pyramid.plot(Sample_Berlin[,1],Sample_Berlin[,2],top.labels=c("", "",""), ppmar=c(2,1,3,1),raxlab = seq(0,30,10),laxlab = seq(0,30,10), xlim = c(30,30),labelcex=1, unit = "",show.values=F, labels = labels2, lxcol = col.pal[2], rxcol = col.pal[3],space = 0.2,gap = 0))
mtext("Number of individuals per age class and gender", side = 1,line = 4.5,at = -33, outer = F, cex = 0.9)


#Densities
s_Berlin <- extract.samples(m_Berlin)
s_Vanuatu <- extract.samples(m_Vanuatu)
s_Berlin_basic <- extract.samples(m_Berlin_basic)
s_Vanuatu_basic <- extract.samples(m_Vanuatu_basic)

#Densities Tanna

dens <- density(inv_logit(s_Vanuatu_basic$p))
x1 <- min(which(dens$x >= quantile(inv_logit(s_Vanuatu_basic$p), 0)))  
x2 <- max(which(dens$x <  quantile(inv_logit(s_Vanuatu_basic$p), 1)))
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

#Densities Berlin
dens <- density(inv_logit(s_Berlin_basic$p))
x1 <- min(which(dens$x >= quantile(inv_logit(s_Berlin_basic$p), 0)))  
x2 <- max(which(dens$x <  quantile(inv_logit(s_Berlin_basic$p), 1)))
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

#dev.off()