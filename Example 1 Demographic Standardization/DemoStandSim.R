
# Generalizing description: Cross-cultural comparisons and demographic standardization
# Simulation Example
# The aim here is to simulate and analyze exemplary data sets generated from different processes
# a) Disparities arise only from real demographic differences among populations (e.g. one is older)
# b) Disparities arise only from differences in sampling procedure (e.g. gender of researcher influences gender of volunteer participants)

#Load packages
library(rethinking)
library(scales)
library(sn)
library(plotrix)
library(RColorBrewer)

N <- 500   # Sample size
Age_range <- c(1:90) #Age range of participants


####                     
###
##
# a) Real differences
##
###
####

#Create population array [pop, sex, age]
D_popdiff <- array(NA,c(2,2,length(Age_range)))

#Exponential distribution (similar to many growing populations)
D_popdiff[1,1, ] <- dexp(seq(0, 100, 1), 0.04)[Age_range]
D_popdiff[1,2, ] <- dexp(seq(0, 100, 1), 0.04)[Age_range]

#Skewed normal distribution (similar to many shrinking populations)
D_popdiff[2,1, ] <- dsn(1:100, xi = 30, omega = 50, alpha = 1)[Age_range]
D_popdiff[2,2, ] <- dsn(1:100, xi = 30, omega = 50, alpha = 1)[Age_range]


#Generate data
#Construct dataframe
d <- data.frame(id = 1:(2*N), soc_id = c(rep(1,N), rep(2,N)), age = NA, gender = NA, outcome = NA)

#Simulate ages from above distributions
for (pop_id in 1:2) {
  d$age[d$soc_id==pop_id] <- sample(Age_range, N, replace = TRUE, prob = D_popdiff[pop_id,1,])
}

#Simulate Genders
p_male <- 0.5
d$gender[d$soc_id==1] <- sample(c(1,2),N, replace = TRUE, prob = c(p_male, 1-p_male)) 
d$gender[d$soc_id==2] <- sample(c(1,2),N, replace = TRUE, prob = c(p_male, 1-p_male)) 

#Create demographic statistics of samples
SampleD_popdiff <- array(NA,c(2,2,length(Age_range)))
for (pop_id in 1:2) {
  for (gender in 1:2) {
    for (i in Age_range) {
      SampleD_popdiff[pop_id, gender, i] <- length(which(d$age[d$soc_id == pop_id]==i & d$gender[d$soc_id == pop_id] == gender))
    }  
  }
}

# Create "Real" cultural differences: (Demography-independent) Probabilities in both populations on logit scale
p_logit_culture <-c(-3, -1.5)

#Effects of age and gender on prosociality
b_age <- 0.04
b_gender <- 2

#Generate observations
for(i in 1:(2*N)) d$outcome[i] <- rbinom(1, 1, inv_logit(p_logit_culture[d$soc_id[i]] + b_age*d$age[i] + b_gender*(d$gender[i]-1)) )


#Compute true expected prosociality in both populations
True_Value_popdiff <- c()
for (pop in 1:2) {
  expect_pos = 0
  total = 0
  for (a in 1:2){
    for (b in 1:max(Age_range)){
      total = total + D_popdiff[pop,a,b];
      expect_pos = expect_pos + D_popdiff[pop,a,b] * inv_logit(p_logit_culture[pop] + b_gender*(a-1) + b_age*b);
    }
  }  
  True_Value_popdiff[pop] = expect_pos / total
}


#Create data lists for both populations 
d1_popdiff <- list(N = N,                    #Sample size
                   MA = max(Age_range),              #Maximum age in data set
                   gender = d$gender[d$soc_id==1],   #Gender
                   age = d$age[d$soc_id==1],         #Age
                   outcome = d$outcome[d$soc_id==1]  #Choices in dictator game
)

#Population demography of both sites: Stan model, as we coded it, requires number (as integer) of individuals in each cell (age/gender combination) for poststratification,
#so we multiply by large number. This does not make a difference because we only care about the relative sizes of cells in the target population

d1_popdiff$Pop <-  D_popdiff[2,, ]* 1e9
mode(d1_popdiff$Pop) <- 'integer'


d2_popdiff <- list(N = N,
                   MA = max(Age_range),
                   gender = d$gender[d$soc_id==2], 
                   age = d$age[d$soc_id==2],
                   outcome = d$outcome[d$soc_id==2]
)

d2_popdiff$Pop <-  D_popdiff[1,, ]* 1e9
mode(d2_popdiff$Pop) <- 'integer'





####                     
###
##
# b) Differences in sampling of ages/genders
##
###
####


#Create population array [pop, sex, age]
D_samplediff <- array(NA,c(2,2,length(Age_range)))

#Exponential distribution (similar to many growing populations)
D_samplediff[1,1, ] <- dexp(seq(0, 100, 1), 0.04)[Age_range]
D_samplediff[1,2, ] <- dexp(seq(0, 100, 1), 0.04)[Age_range]

D_samplediff[2,1, ] <- dexp(seq(0, 100, 1), 0.04)[Age_range]
D_samplediff[2,2, ] <- dexp(seq(0, 100, 1), 0.04)[Age_range]


#Generate data
#Construct dataframe
d <- data.frame(id = 1:(2*N), soc_id = c(rep(1,N), rep(2,N)), age = NA, gender = NA, outcome = NA)

#Simulate ages from above distributions
for (pop_id in 1:2) {
  d$age[d$soc_id==pop_id] <- sample(Age_range, N, replace = TRUE, prob = D_samplediff[pop_id,1,])
}

#Simulate Genders, we assume that genders have different probabilities to be sampled in both populations

d$gender[d$soc_id==1] <- sample(c(1,2),N, replace = TRUE, prob = c(0.3, 0.7)) 
d$gender[d$soc_id==2] <- sample(c(1,2),N, replace = TRUE, prob = c(0.8, 0.2)) 


#Create demographic statistics of samples
SampleD_samplediff <- array(NA,c(2,2,length(Age_range)))
for (pop_id in 1:2) {
  for (gender in 1:2) {
    for (i in Age_range) {
      SampleD_samplediff[pop_id, gender, i] <- length(which(d$age[d$soc_id == pop_id]==i & d$gender[d$soc_id == pop_id] == gender))
    }  
  }
}



#Generate observations
for(i in 1:(2*N)) d$outcome[i] <- rbinom(1, 1, inv_logit(p_logit_culture[d$soc_id[i]] + b_age*d$age[i] + b_gender*(d$gender[i]-1)) )

#Compute true expected prosociality in both populations
True_Value_samplediff <- c()
for (pop in 1:2) {
  expect_pos = 0
  total = 0
  for (a in 1:2){
    for (b in 1:max(Age_range)){
      total = total + D_samplediff[pop,a,b];
      expect_pos = expect_pos + D_samplediff[pop,a,b] * inv_logit(p_logit_culture[pop] + b_gender*(a-1) + b_age*b);
    }
  }  
  True_Value_samplediff[pop] = expect_pos / total
}






# Prepare lists for stan
d1_samplediff <- list(N = N,
                      MA = max(Age_range),
                      gender = d$gender[d$soc_id==1], 
                      age = d$age[d$soc_id==1],
                      outcome = d$outcome[d$soc_id==1]
)

d2_samplediff <- list(N = N,
                      MA = max(Age_range),
                      gender = d$gender[d$soc_id==2], 
                      age = d$age[d$soc_id==2],
                      outcome = d$outcome[d$soc_id==2]
)


#Population demography of both sites: Stan model, as we coded it, requires number (as integer) of individuals 
#in each cell (age/gender combination) for poststratification, so we multiply by large number. 
#This does not make a difference because we only care about the relative sizes of cells in the target population

d1_samplediff$Pop <-  D_samplediff[1,,] * 1e9
mode(d1_samplediff$Pop) <- 'integer'


d2_samplediff$Pop <-  D_samplediff[2,,] * 1e9
mode(d2_samplediff$Pop) <- 'integer'




#Here, we code the stan model for a simple Bernoulli model. This is used for "empirical" estimates


Basic_sim <- "

data {
  int N;
  int outcome[N];
}

parameters {
  real p;            
}

model {
 outcome ~ bernoulli_logit(p);
}

"



#Here, we code the stan model for our Gaussian Process Multilevel Regression with Poststratification

MRP_sim <- "

//Function for Gaussian Process kernel

functions{
  matrix GPL(int K, real C, real D, real S){
   matrix[K,K] Rho;                       
   real KR;                               
   KR = K;                             

   for(i in 1:(K-1)){
   for(j in (i+1):K){  
    Rho[i,j] = C * exp(-D * ( (j-i)^2 / KR^2) );           
    Rho[j,i] = Rho[i,j];                       
    }}

   for (i in 1:K){
    Rho[i,i] = 1;                               
    }

   return S*cholesky_decompose(Rho);
   }
}

data {
  int N;
  int MA;
  int age[N];
  int outcome[N];
  int gender[N];
  int<lower = 0> Pop[2,MA];   
  }

parameters {
  vector[2] alpha;            //Gender-specific intercepts
  matrix[2,MA] age_effect;    //Matrix for Gaussian process age effects
  
  //Here we define the Control parameters for the Gaussian processes; they determine how covariance changes with increasing distance in age
  real<lower=0> eta[2];
  real<lower=0> sigma[2];
  real<lower=0, upper=1> rho[2];
}

model {
  vector[N] p;

  //Priors
  alpha ~ normal(0, 2);
  eta ~ exponential(2);
  sigma ~ exponential(1);
  rho ~ beta(10, 1);

  //We compute age-specific offsets for each sex
  for ( i in 1:2){
   age_effect[i,] ~ multi_normal_cholesky( rep_vector(0, MA) , GPL(MA, rho[i], eta[i], sigma[i]) );
  }  

  //This is the linear model: Choice probabilities are composed of (gender-specific) intercept and gender-specific offset for each age category
  for ( i in 1:N ) {
   p[i] = alpha[gender[i]] + age_effect[gender[i],age[i]];
  }
  
 //Binomial likelihood for observed choices
   outcome ~ binomial_logit(1, p);
}

 //We use the generated quantities section for the poststratification
 
generated quantities{
   real<lower = 0, upper = 1> p_target;  // This is value for p in the target population 
   real expect_pos = 0;
   int total = 0;
   
  vector[MA] pred_p_m;
   vector[MA] pred_p_f;
   
  //Here we compute predictions for each age and gender class
  pred_p_m = inv_logit(alpha[1] + age_effect[1,]');
  pred_p_f = inv_logit(alpha[2] + age_effect[2,]');


   //Here we do the actual poststratification
      for (a in 1:2)
         for (b in 1:MA){
          total += Pop[a,b];
          expect_pos += Pop[a,b] * inv_logit(alpha[a] + age_effect[a,b]);
      }
    p_target = expect_pos / total;
}

"


#Pass code to Rstan, run models and extract the samples from the posterior
#Empirical Estimates
m_popdiff1_basic <- stan( model_code  = Basic_sim , data=d1_popdiff ,iter = 5000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.99, max_treedepth = 13))  
m_popdiff2_basic <- stan( model_code  = Basic_sim , data=d2_popdiff ,iter = 5000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.99, max_treedepth = 13))  
s_popdiff1_basic <- extract.samples(m_popdiff1_basic)
s_popdiff2_basic <- extract.samples(m_popdiff2_basic)


m_samplediff1_basic <- stan( model_code  = Basic_sim , data=d1_samplediff ,iter = 5000, cores = 4, seed=4, chains=4, control = list(adapt_delta=0.99, max_treedepth = 13))  
m_samplediff2_basic <- stan( model_code  = Basic_sim , data=d2_samplediff ,iter = 5000, cores = 4, seed=4, chains=4, control = list(adapt_delta=0.99, max_treedepth = 13))  
s_samplediff1_basic <- extract.samples(m_samplediff1_basic)
s_samplediff2_basic <- extract.samples(m_samplediff2_basic)

#Mister P
m_popdiff1 <- stan( model_code  = MRP_sim , data=d1_popdiff ,iter = 5000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.99, max_treedepth = 13))  
m_popdiff2 <- stan( model_code  = MRP_sim , data=d2_popdiff ,iter = 5000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.99, max_treedepth = 13))  
s_popdiff1 <- extract.samples(m_popdiff1)
s_popdiff2 <- extract.samples(m_popdiff2)


m_samplediff1 <- stan( model_code  = MRP_sim , data=d1_samplediff ,iter = 5000, cores = 4, seed=4, chains=4, control = list(adapt_delta=0.99, max_treedepth = 13))  
m_samplediff2 <- stan( model_code  = MRP_sim , data=d2_samplediff ,iter = 5000, cores = 4, seed=4, chains=4, control = list(adapt_delta=0.99, max_treedepth = 13))  
s_samplediff1 <- extract.samples(m_samplediff1)
s_samplediff2 <- extract.samples(m_samplediff2)






######
####
###
##
# Plotting Code
##
###
####
####





col.pal <- brewer.pal(8, "Dark2")
seqoverall <- seq


###
##
# Code for Figure S2 in the appendix
##
###


#graphics.off()
#png("AgeCurves.png", res = 900, height = 16, width = 24, units = "cm")


par(mfrow = c(2,2),
    mar = c(1,2,1.2,0.5), 
    oma = c(2.2,5.1,1,0))

###
##
# Population Differences
##
###
samp <- extract.samples(m_popdiff1)
samp_pred_p <- apply(samp$pred_p_m,2,quantile,probs=c(0.05,0.5,0.95))
Lower <- samp_pred_p[1,]
Median <- samp_pred_p[2,]
Higher <- samp_pred_p[3,]

Observed_m <- ifelse(c(1:d1_popdiff$MA) %in% unique(d1_popdiff$age[d1_popdiff$gender==1]),"Observed Age","Unobserved Age")
Observed_f <- ifelse(c(1:d1_popdiff$MA) %in% unique(d1_popdiff$age[d1_popdiff$gender==2]),"Observed Age","Unobserved Age")

plot(NA, type = "l", xlim = c(0,max(Age_range)), ylim = c(0,1), xlab = "", ylab = "", lwd = 2, lty = 2, main= "Population I")
points(Median, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.3)), pch = 16)
segments(x0 = 1:d1_popdiff$MA, y0 = Lower, x1=  1:d1_popdiff$MA, y1 = Higher, lwd = 2, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.3)))

legend("topleft", c("Male","Female"), col = c(alpha(col.pal[2],alpha = 1), alpha(col.pal[3],alpha = 1)), pch = c(16,17), lwd = 2, lty=1, bty="n", cex = 1.2)
samp_pred_p <- apply(samp$pred_p_f,2,quantile,probs=c(0.05,0.5,0.95))
Lower <- samp_pred_p[1,]
Median <- samp_pred_p[2,]
Higher <- samp_pred_p[3,]

par(new = TRUE)
plot(NA, type = "l", xlim = c(0,max(Age_range)), ylim = c(0,1), xlab = "", ylab = "", lwd = 2, lty = 2, main= "")
points(Median, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.3)), pch = 17)
segments(x0 = 1:d1_popdiff$MA, y0 = Lower, x1=  1:d1_popdiff$MA, y1 = Higher, lwd = 2, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.3)))

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

#dev.off()



###
##
# Code for Figure S1 in the appendix
##
###

#graphics.off()
#png("DemostandLarge.png", res = 900, height = 22, width = 25, units = "cm")

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
dens <- density(inv_logit(s_popdiff1_basic$p))
x1 <- min(which(dens$x >= quantile(inv_logit(s_popdiff1_basic$p), 0)))  
x2 <- max(which(dens$x <  quantile(inv_logit(s_popdiff1_basic$p), 1)))
plot(dens, xlim = c(0,1), ylim = c(0,25), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[1],alpha = 0.9), border = NA))
abline(v = True_Value_popdiff[1], lty = 2, lwd = 2)
par(new=TRUE)
dens <- density(s_popdiff1$p_target)
x1 <- min(which(dens$x >= quantile(s_popdiff1$p_target, 0)))  
x2 <- max(which(dens$x <  quantile(s_popdiff1$p_target, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,25), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[1],alpha = 0.3), border = NA))
mtext("Density", side = 2,line = 3, outer = F, cex = 1)
legend("topright", c("Empirical","Poststr. to Pop II"), col = c(alpha(col.pal[1],alpha = 0.9),alpha(col.pal[1],alpha = 0.3)), lwd = 6, bty="n", cex = 0.9)

dens <- density(inv_logit(s_popdiff2_basic$p))
x1 <- min(which(dens$x >= quantile(inv_logit(s_popdiff2_basic$p), 0)))  
x2 <- max(which(dens$x <  quantile(inv_logit(s_popdiff2_basic$p), 1)))
plot(dens, xlim = c(0,1), ylim = c(0,25), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[4],alpha = 0.9), border = NA))
abline(v = True_Value_popdiff[2], lty = 2, lwd = 2)
par(new=TRUE)
dens <- density(s_popdiff2$p_target)
x1 <- min(which(dens$x >= quantile(s_popdiff2$p_target, 0)))
x2 <- max(which(dens$x <  quantile(s_popdiff2$p_target, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,25), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[4],alpha = 0.3), border = NA))
mtext("Probability of choosing prosocial option", side = 1,line = 2.9,at=-0.1, outer = F, cex = 0.9)
legend("topleft", c("Empirical","Poststr. to Pop I"), col = c(alpha(col.pal[4],alpha = 0.9),alpha(col.pal[4],alpha = 0.3)), lwd = 6, bty="n", cex =  0.9)


#Densities 2
dens <- density(inv_logit(s_samplediff1_basic$p))
x1 <- min(which(dens$x >= quantile(inv_logit(s_samplediff1_basic$p), 0)))  
x2 <- max(which(dens$x <  quantile(inv_logit(s_samplediff1_basic$p), 1)))
plot(dens, xlim = c(0,1), ylim = c(0,25), type="n", yaxt = "n",ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[5],alpha = 0.9), border = NA))
abline(v = True_Value_samplediff[1], lty = 2, lwd = 2)
par(new=TRUE)
dens <- density(s_samplediff1$p_target)
x1 <- min(which(dens$x >= quantile(s_samplediff1$p_target, 0)))  
x2 <- max(which(dens$x <  quantile(s_samplediff1$p_target, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,25), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[5],alpha = 0.3), border = NA))
legend("topright", c("Empirical","Poststr. to Pop I"), col = c(alpha(col.pal[5],alpha = 0.9),alpha(col.pal[5],alpha = 0.3)), lwd = 6, bty="n", cex = 0.9)

dens <- density(inv_logit(s_samplediff2_basic$p))
x1 <- min(which(dens$x >= quantile(inv_logit(s_samplediff2_basic$p), 0)))  
x2 <- max(which(dens$x <  quantile(inv_logit(s_samplediff2_basic$p), 1)))
plot(dens, xlim = c(0,1), ylim = c(0,25), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[6],alpha = 0.9), border = NA))
abline(v = True_Value_samplediff[2], lty = 2, lwd = 2)
par(new=TRUE)
dens <- density(s_samplediff2$p_target)
x1 <- min(which(dens$x >= quantile(s_samplediff2$p_target, 0)))
x2 <- max(which(dens$x <  quantile(s_samplediff2$p_target, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,25), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[6],alpha = 0.3), border = NA))

legend("topleft", c("Empirical","Poststr. to Pop II"), col = c(alpha(col.pal[6],alpha = 0.9),alpha(col.pal[6],alpha = 0.3)), lwd = 6, bty="n", cex = 0.9)
mtext("Probability of choosing prosocial option", side = 1,line = 2.9,at=-0.1, outer = F, cex = 0.9)

#dev.off()




