
# Example 1: Demographic Standardization
# The aim here is to simulate and analyze exemplary data sets generated from different processes
# a) Disparities arise only from real demographic differences among populations (e.g. one is older/more men)
# b) Disparities arise only from differences in sampling procedure 

library(rethinking)
library(scales)
library(sn)
library(plotrix)

N <- 500
Age_range <- c(1:90)


####                     
###
##
# a) Real differences
##
###
####


#Create population array [pop, sex, age]
D_popdiff <- array(NA,c(2,2,length(Age_range)))

X <- 1:100

age_seq <- seq(0, 100, 1)

#Exponential distribution (similar to many growing populations)
D_popdiff[1,1, ] <- dexp(age_seq, 0.04)[Age_range]
D_popdiff[1,2, ] <- dexp(age_seq, 0.04)[Age_range]

#Skewed normal distribution (similar to many shrinking populations)
D_popdiff[2,1, ] <- dsn(X, xi = 30, omega = 50, alpha = 1)[Age_range]
D_popdiff[2,2, ] <- dsn(X, xi = 30, omega = 50, alpha = 1)[Age_range]


#Generate data
#Construct dataframe
d <- data.frame(id = 1:(2*N), soc_id = c(rep(1,N), rep(2,N)), age = NA, gender = NA, outcome = NA)

#Simulate ages from above distributions
for (pop_id in 1:2) {
  d$age[d$soc_id==pop_id] <- sample(Age_range, N, replace = TRUE, prob = D_popdiff[pop_id,1,])
}

#Simulate Genders

d$gender[d$soc_id==1] <- sample(c(1,2),N, replace = TRUE, prob = c(0.5, 0.5)) 
d$gender[d$soc_id==2] <- sample(c(1,2),N, replace = TRUE, prob = c(0.5, 0.5)) 

SampleD_popdiff <- array(NA,c(2,2,length(Age_range)))

for (pop_id in 1:2) {
  for (gender in 1:2) {
    for (i in Age_range) {
    SampleD_popdiff[pop_id, gender, i] <- length(which(d$age[d$soc_id == pop_id]==i & d$gender[d$soc_id == pop_id] == gender))
  }  
 }
}


p_logit_culture <-c(-3, -1.5)

#Effects of age and gender
b_age <- 0.04
b_gender <- 2

#Generate observations
stand_age<- c()
for (i in 1:(2*N)) stand_age[i] <- (d$age[i] - mean(d$age[d$soc_id==1]))/sd(d$age[d$soc_id==1])

for(i in 1:(2*N)) d$outcome[i] <- rbinom(1, 1, inv_logit(p_logit_culture[d$soc_id[i]] + b_age*d$age[i] + b_gender*(d$gender[i]-1)) )


phi_popdiff <- c()
for (pop in 1:2) {
  expect_pos = 0
  total = 0
  for (a in 1:2){
    for (b in 1:max(Age_range)){
      total = total + D_popdiff[pop,a,b];
      expect_pos = expect_pos + D_popdiff[pop,a,b] * inv_logit(p_logit_culture[pop] + b_gender*(a-1) + b_age*b);
    }
  }  
  phi_popdiff[pop] = expect_pos / total
}







d1_popdiff <- list(N = N,
           MA = max(Age_range),
           gender = d$gender[d$soc_id==1], 
           age = d$age[d$soc_id==1],
           outcome = d$outcome[d$soc_id==1]
)

d2_popdiff <- list(N = N,
           MA = max(Age_range),
           gender = d$gender[d$soc_id==2], 
           age = d$age[d$soc_id==2],
           outcome = d$outcome[d$soc_id==2]
)




d1_popdiff$P_same  <- SampleD_popdiff[1,,]
d1_popdiff$P_other <-  SampleD_popdiff[2,,]


d2_popdiff$P_same  <- SampleD_popdiff[2,,]
d2_popdiff$P_other <-  SampleD_popdiff[1,,]




{
  
  m2a_MRP_GP_gender_same <- "
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
  int<lower = 0> P_same[2,MA];   // Here we enter data matrix with demographic constitution of target population
  int<lower = 0> P_other[2,MA];   // Here we enter data matrix with demographic constitution of target population

  }

parameters {
  vector[2] alpha;
  matrix[2,MA] age_effect;    //Matrix for GP age effects
  
  real<lower=0> eta[2];
  real<lower=0> sigma[2];
  real<lower=0, upper=1> rho[2];
}

model {
  vector[N] p;

  alpha ~ normal(0, 2);

  eta ~ exponential(2);
  sigma ~ exponential(1);
  rho ~ beta(10, 1);


  for ( i in 1:2){
   age_effect[i,] ~ multi_normal_cholesky( rep_vector(0, MA) , GPL(MA, rho[i], eta[i], sigma[i]) );
  }  

  for ( i in 1:N ) {
   p[i] = alpha[gender[i]] + age_effect[gender[i],age[i]];
  }

 outcome ~ binomial_logit(1, p);
}
 
generated quantities{
   real expect_pos = 0;
   real<lower = 0, upper = 1> phi;  // This is value for p in the source population 
   real<lower = 0, upper = 1> psi;  // This is value for p in the target population 

   int total = 0;
   vector[MA] pred_p_m;
   vector[MA] pred_p_f;

   for (a in 1:2)
        for (b in 1:MA){
         total += P_same[a,b];
         expect_pos += P_same[a,b] * inv_logit(alpha[a] + age_effect[a,b]);
     }
     
   phi = expect_pos / total;
   
    total = 0;
    expect_pos = 0;
    
       for (a in 1:2)
         for (b in 1:MA){
          total += P_other[a,b];
          expect_pos += P_other[a,b] * inv_logit(alpha[a] + age_effect[a,b]);
      }
      
    psi = expect_pos / total;
   
  pred_p_m = inv_logit(alpha[1] + age_effect[1,]');
  pred_p_f = inv_logit(alpha[2] + age_effect[2,]');
}

"

}



m_popdiff1 <- stan( model_code  = m2a_MRP_GP_gender_same , data=d1_popdiff ,iter = 4000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.95, max_treedepth = 13))  

m_popdiff2 <- stan( model_code  = m2a_MRP_GP_gender_same , data=d2_popdiff ,iter = 4000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.95, max_treedepth = 13))  

s_popdiff1 <- extract.samples(m_popdiff1)
s_popdiff2 <- extract.samples(m_popdiff2)





####                     
###
##
# b) Differences in sampling of ages/genders
##
###
####


#Create population array [pop, sex, age]
D_samplediff <- array(NA,c(2,2,length(Age_range)))

X <- 1:100

age_seq <- seq(0, 100, 1)

#Exponential distribution (similar to many growing populations)
D_samplediff[1,1, ] <- dexp(age_seq, 0.04)[Age_range]
D_samplediff[1,2, ] <- dexp(age_seq, 0.04)[Age_range]

#Skewed normal distribution (similar to many shrinking populations)
D_samplediff[2,1, ] <- dexp(age_seq, 0.04)[Age_range]
D_samplediff[2,2, ] <-dexp(age_seq, 0.04)[Age_range]



#Generate data
#Construct dataframe
d <- data.frame(id = 1:(2*N), soc_id = c(rep(1,N), rep(2,N)), age = NA, gender = NA, outcome = NA)

#Simulate ages from above distributions
for (pop_id in 1:2) {
  d$age[d$soc_id==pop_id] <- sample(Age_range, N, replace = TRUE, prob = D_samplediff[pop_id,1,])
}

#Simulate Genders

d$gender[d$soc_id==1] <- sample(c(1,2),N, replace = TRUE, prob = c(0.3, 0.7)) 
d$gender[d$soc_id==2] <- sample(c(1,2),N, replace = TRUE, prob = c(0.8, 0.2)) 


SampleD_samplediff <- array(NA,c(2,2,length(Age_range)))

for (pop_id in 1:2) {
  for (gender in 1:2) {
    for (i in Age_range) {
      SampleD_samplediff[pop_id, gender, i] <- length(which(d$age[d$soc_id == pop_id]==i & d$gender[d$soc_id == pop_id] == gender))
    }  
  }
}




#Generate observations

#Effect of "culture" that's independent of demography
p_logit_culture <-c(-3, -1.5)

#Effects of age and gender
b_age <- 0.04
b_gender <- 2


for(i in 1:(2*N)) d$outcome[i] <- rbinom(1, 1, inv_logit(p_logit_culture[d$soc_id[i]] + b_age*d$age[i] + b_gender*(d$gender[i]-1)) )

#Compute expected prosociality in both populations
phi_samplediff <- c()
for (pop in 1:2) {
  expect_pos = 0
  total = 0
  for (a in 1:2){
    for (b in 1:max(Age_range)){
      total = total + D_samplediff[pop,a,b];
      expect_pos = expect_pos + D_samplediff[pop,a,b] * inv_logit(p_logit_culture[pop] + b_gender*(a-1) + b_age*b);
    }
  }  
  phi_samplediff[pop] = expect_pos / total
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


#
d1_samplediff$P_same <- SampleD_samplediff[1,,]
d1_samplediff$P_other <-  D_samplediff[1,,] * 1e9
mode(d1_samplediff$P_other) <- 'integer'


d2_samplediff$P_same <- SampleD_samplediff[2,,]
d2_samplediff$P_other <-  D_samplediff[2,,] * 1e9
mode(d2_samplediff$P_other) <- 'integer'







m_samplediff1 <- stan( model_code  = m2a_MRP_GP_gender_same , data=d1_samplediff ,iter = 4000, cores = 4, seed=4, chains=1, control = list(adapt_delta=0.95, max_treedepth = 13))  

m_samplediff2 <- stan( model_code  = m2a_MRP_GP_gender_same , data=d2_samplediff ,iter = 4000, cores = 4, seed=4, chains=1, control = list(adapt_delta=0.95, max_treedepth = 13))  

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



dev.off()






graphics.off()
png("DemostandLarge.png", res = 900, height = 22, width = 25, units = "cm")



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
plot(dens, xlim = c(0,1), ylim = c(0,25), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[1],alpha = 0.9), border = NA))
abline(v = phi_popdiff[1], lty = 2, lwd = 2)
par(new=TRUE)
dens <- density(s_popdiff1$psi)
x1 <- min(which(dens$x >= quantile(s_popdiff1$psi, 0)))  
x2 <- max(which(dens$x <  quantile(s_popdiff1$psi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,25), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[1],alpha = 0.3), border = NA))
mtext("Density", side = 2,line = 3, outer = F, cex = 1)

legend("topright", c("Empirical","Poststr. to Pop II"), col = c(alpha(col.pal[1],alpha = 0.9),alpha(col.pal[1],alpha = 0.3)), lwd = 6, bty="n", cex = 0.9)


dens <- density(s_popdiff2$phi)
x1 <- min(which(dens$x >= quantile(s_popdiff2$phi, 0)))  
x2 <- max(which(dens$x <  quantile(s_popdiff2$phi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,25), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[4],alpha = 0.9), border = NA))
abline(v = phi_popdiff[2], lty = 2, lwd = 2)

par(new=TRUE)

dens <- density(s_popdiff2$psi)
x1 <- min(which(dens$x >= quantile(s_popdiff2$psi, 0)))
x2 <- max(which(dens$x <  quantile(s_popdiff2$psi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,25), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[4],alpha = 0.3), border = NA))

mtext("Probability of choosing prosocial option", side = 1,line = 2.9,at=-0.1, outer = F, cex = 0.9)

legend("topleft", c("Empirical","Poststr. to Pop I"), col = c(alpha(col.pal[4],alpha = 0.9),alpha(col.pal[4],alpha = 0.3)), lwd = 6, bty="n", cex =  0.9)







#Densities 2


dens <- density(s_samplediff1$phi)
x1 <- min(which(dens$x >= quantile(s_samplediff1$phi, 0)))  
x2 <- max(which(dens$x <  quantile(s_samplediff1$phi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,25), type="n", yaxt = "n",ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[5],alpha = 0.9), border = NA))
abline(v = phi_samplediff[1], lty = 2, lwd = 2)
par(new=TRUE)
dens <- density(s_samplediff1$psi)
x1 <- min(which(dens$x >= quantile(s_samplediff1$psi, 0)))  
x2 <- max(which(dens$x <  quantile(s_samplediff1$psi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,25), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[5],alpha = 0.3), border = NA))

legend("topright", c("Empirical","Poststr. to Pop I"), col = c(alpha(col.pal[5],alpha = 0.9),alpha(col.pal[5],alpha = 0.3)), lwd = 6, bty="n", cex = 0.9)



dens <- density(s_samplediff2$phi)
x1 <- min(which(dens$x >= quantile(s_samplediff2$phi, 0)))  
x2 <- max(which(dens$x <  quantile(s_samplediff2$phi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,25), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[6],alpha = 0.9), border = NA))
abline(v = phi_samplediff[2], lty = 2, lwd = 2)

par(new=TRUE)

dens <- density(s_samplediff2$psi)
x1 <- min(which(dens$x >= quantile(s_samplediff2$psi, 0)))
x2 <- max(which(dens$x <  quantile(s_samplediff2$psi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,25), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[6],alpha = 0.3), border = NA))

legend("topleft", c("Empirical","Poststr. to Pop II"), col = c(alpha(col.pal[6],alpha = 0.9),alpha(col.pal[6],alpha = 0.3)), lwd = 6, bty="n", cex = 0.9)

mtext("Probability of choosing prosocial option", side = 1,line = 2.9,at=-0.1, outer = F, cex = 0.9)

dev.off()




