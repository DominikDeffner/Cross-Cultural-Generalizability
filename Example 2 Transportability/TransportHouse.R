
#Generalizing experimental results: Transportability of Causal Effects
#Real data example based on House et al., 2020
#The script prepares the data, runs the transport analysis and creates plot in Fig.5 in the manuscript

#Load some packages and set working directory to main folder
library(readr)
library(rethinking)
library(plotrix)
library(scales)
library(RColorBrewer)
setwd("~/GitHub/Cross-Cultural-Generalizability")

#Load the data from House et al., 2020
data <- read.csv("data/Model_6a_6b_6c_6d_data.csv")

#Create new variable for experimental condition, round ages and recode gender such that 1 = male and 2 = female
data$condition <- sapply(1:nrow(data), function(x) which(c(data$CONDITION_1_1yes[x],data$CONDITION_2_1yes[x],data$CONDITION_3_1yes[x]) == 1 ))
data <- data[,c("SUBJECT_ID","GENDER_1female","fieldid","AGE_in_years","condition","T1_choice_1yes")]
data$choice <- data$T1_choice_1yes
data$T1_choice_1yes <- NULL
data$age <- round(data$AGE_in_years)
data$gender <- data$GENDER_1female + 1

#Exclude data from "both ok", we focus on comparison between "Generous" and "Selfish" conditions
data <- subset(data, data$condition != 3)

#Create matrix with demography of study samples
Demo <- matrix(0, length(unique(data$fieldid)), length(unique(data$age)))
for (j in sort(unique(data$fieldid))) {
  for (i in min(data$age):max(data$age)) {
      Demo[j, which(min(data$age):max(data$age) == i)] <- length(which(data$fieldid==j & data$age==i))
  }
}          


#Prepare data list as input for Stan 
d_list <- list(N = nrow(data),                         #Number of unique choices/participants    
               N_pop = length(unique(data$fieldid)),   #Number of field sites
               MA = 12,                                #Maximum age (recoded, see below)
               pop_id = data$fieldid,                  #Numerical field id (1=Berlin, 2=La Plata, 3=Phoenix, 4=Pune, 5=Shuar, 6=Wiichi)
               age = data$age - 3,                     #We recode age, which ranges from 4-15 in the data, to 1-12 to make indexing easier
               condition = data$condition-1,           #Dummy code conditions
               outcome = data$choice,                  #Dictator game choice
               gender = data$gender,                   #Gender: 1 = male and 2 = female
               Demo = Demo,                            #Demography of 6 samples
               Ref = 6)                                #This is the target population for the transport (we choose the Wiichi as an example)


# The next section invokes the stan model code in the example 2 folder
# The first model ("model_basic") is a simple fixed effects model that provides unbiased or "empirical" estimates of the 
# causal effect in each population. As the causal effect is defined on the outcome scale (difference in choice probabilities)
# and does not directly correspond to any model parameter,it is computed from the model parameters in the generated quantities section

#The second model ("model_transport") uses Gaussian processes to compute age-specific causal effects for each population.
#It then uses these strata-specific effects to transport effects to any arbitrary target population (here Ref = 6, which corresponds to the Wiichi)


#Pass code to Rstan and run models
m_empirical <- stan(file  = "Example 2 Transportability/model_basic.stan" , data= d_list ,iter = 5000, cores = 4, chains=4, control = list(adapt_delta=0.99, max_treedepth = 15))  
m_transport <- stan(file  = "Example 2 Transportability/model_transport.stan" , data= d_list ,iter = 5000, cores = 4,chains=4, control = list(adapt_delta=0.99, max_treedepth = 15))  


#Extract samples from the posterior
s_basic <- extract.samples(m_empirical)
s <- extract.samples(m_transport)


#Code for Fig.5 in the manuscript

col.pal <- brewer.pal(8, "Dark2") #create a palette which you loop over for corresponding values
seqoverall <- seq
Society <- c("Berlin (GER)","La Plata (ARG)","Phoenix (USA)", "Pune (IND)", "Shuar (ECU)", "Wichí (ARG)")

#Uncomment to save plot as png file in working directory
#graphics.off()
#png("Transport.png", res = 900, height = 12, width = 16, units = "cm")

par(mfrow= c(3,2),
    oma=c(3,0,3,0),
    mar=c(2,0,1,1))

for (i in 1:6) {
  
  dens <- density(s_basic$empirical_p[,i])
  x1 <- min(which(dens$x >= quantile(s_basic$empirical_p[,i], 0)))
  x2 <- max(which(dens$x <  quantile(s_basic$empirical_p[,i], 1)))
  plot(dens, xlim = c(-0.3,0.9), ylim = c(0,max(dens$y)), type="n", ann = FALSE, bty = "n", yaxt = "n")
  with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[i],alpha = 0.7), border = NA))
  abline(v = 0, lty = 2)
  mtext(Society[i], side = 2, line = -3, cex = 0.7)
  if (i ==1) legend("topright", c("Empirical estimates", "Transported to the Wichí"), col = c(alpha("black",alpha = 0.7),alpha("black",alpha = 0.2)), cex = 0.7, lwd = 8, lty = 1, bty = "n")
  par(new = TRUE)
  
  dens <- density(s$transport_p[,i])
  x1 <- min(which(dens$x >= quantile(s$transport_p[,i], 0)))
  x2 <- max(which(dens$x <  quantile(s$transport_p[,i], 1)))
  plot(dens, xlim = c(-0.3,0.9), ylim = c(0,max(dens$y)), type="n", ann = FALSE, bty = "n", yaxt = "n")
  with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[i],alpha = 0.2), border = NA))
  abline(v = 0, lty = 2)

}

mtext("Effect of norm prime on prosocial choices", side = 1,line = 1.5,outer = TRUE, cex = 0.8)
mtext("'Transport' of causal effects across populations", side = 3,line = 1.5,outer = TRUE, cex = 1)

#dev.off()

