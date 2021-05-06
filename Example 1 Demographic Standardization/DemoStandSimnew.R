
# Example 1: Demographic Standardization
# The aim here is to simulate and analyze exemplary data sets generated from different processes
# a) Disparities arise only from real demographic differences among populations (e.g. one is older/more men)
# b) Disparities arise only from differences in sampling procedure 
# c) Disparities arise from both processes

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


p_logit_culture <-c(-3, -1)

#Effects of age and gender
b_age <- 0.03
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


d1_popdiff$P_same  <- SampleD_popdiff[2,,]
d1_popdiff$P_other <-  SampleD_popdiff[1,,]




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
  
  real<lower=0> eta;
  real<lower=0> sigma;
  real<lower=0, upper=1> rho;
}

model {
  vector[N] p;

  alpha ~ normal(0, 1);

  eta ~ exponential(2);
  sigma ~ exponential(1);
  rho ~ beta(10, 1);


  for ( i in 1:2){
   age_effect[i,] ~ multi_normal_cholesky( rep_vector(0, MA) , GPL(MA, rho, eta, sigma) );
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



m_popdiff1 <- stan( model_code  = m2a_MRP_GP_gender_same , data=d1_popdiff ,iter = 4000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.95, max_treedepth = 13))  

m_popdiff2 <- stan( model_code  = m2a_MRP_GP_gender_same , data=d1_popdiff ,iter = 4000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.95, max_treedepth = 13))  

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

d$gender[d$soc_id==1] <- sample(c(1,2),N, replace = TRUE, prob = c(0.35, 0.65)) 
d$gender[d$soc_id==2] <- sample(c(1,2),N, replace = TRUE, prob = c(0.75, 0.25)) 


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
p_logit_culture <-c(-3, -1)

#Effects of age and gender
b_age <- 0.03
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
d1_samplediff$P_other <-  D_samplediff[1,,] * 1000000000
mode(d1_samplediff$P_other) <- 'integer'


d2_samplediff$P_same <- SampleD_samplediff[2,,]
d2_samplediff$P_other <-  D_samplediff[2,,] * 1000000000
mode(d2_samplediff$P_other) <- 'integer'







m_samplediff1 <- stan( model_code  = m2a_MRP_GP_gender_same , data=d1_samplediff ,iter = 4000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.95, max_treedepth = 13))  

m_samplediff2 <- stan( model_code  = m2a_MRP_GP_gender_same , data=d2_samplediff ,iter = 4000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.95, max_treedepth = 13))  

s_samplediff1 <- extract.samples(m_samplediff1)
s_samplediff2 <- extract.samples(m_samplediff2)







