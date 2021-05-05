
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
D1 <- array(NA,c(2,2,length(Age_range)))

X <- 1:100

age_seq <- seq(0, 100, 1)

#Exponential distribution (similar to many growing populations)
D1[1,1, ] <- dexp(age_seq, 0.04)[Age_range]
D1[1,2, ] <- dexp(age_seq, 0.04)[Age_range]

#Skewed normal distribution (similar to many shrinking populations)
D1[2,1, ] <- dsn(X, xi = 30, omega = 50, alpha = 1)[Age_range]
D1[2,2, ] <- dsn(X, xi = 30, omega = 50, alpha = 1)[Age_range]


#Generate data
#Construct dataframe
d <- data.frame(id = 1:(2*N), soc_id = c(rep(1,N), rep(2,N)), age = NA, gender = NA, outcome = NA)

#Simulate ages from above distributions
for (pop_id in 1:2) {
  d$age[d$soc_id==pop_id] <- sample(Age_range, N, replace = TRUE, prob = D1[pop_id,1,])
}

#Simulate Genders

d$gender[d$soc_id==1] <- sample(c(1,2),N, replace = TRUE, prob = c(0.5, 0.5)) 
d$gender[d$soc_id==2] <- sample(c(1,2),N, replace = TRUE, prob = c(0.5, 0.5)) 



a11 <- matrix(0, 2, length(Age_range))
for (i in Age_range) {
     a11[1,which(Age_range == i)] <- length(which(d$age[d$soc_id == 1]==i & d$gender[d$soc_id == 1] == 1))
     a11[2,which(Age_range == i)] <- length(which(d$age[d$soc_id == 1]==i & d$gender[d$soc_id == 1] == 2))
}
a12 <- matrix(0, 2, length(Age_range))
for (i in Age_range) {
  a12[1,which(Age_range == i)] <- length(which(d$age[d$soc_id == 2]==i & d$gender[d$soc_id == 2] == 1))
  a12[2,which(Age_range == i)] <- length(which(d$age[d$soc_id == 2]==i & d$gender[d$soc_id == 2] == 2))
}



p_logit_culture <-c(-1, 0)

b_age <- 0.04
b_gender <- 0

#Generate observations
stand_age<- c()
for (i in 1:(2*N)) stand_age[i] <- (d$age[i] - mean(d$age[d$soc_id==1]))/sd(d$age[d$soc_id==1])

for(i in 1:(2*N)) d$outcome[i] <- rbinom(1, 1, inv_logit(p_logit_culture[d$soc_id[i]] + b_age*d$age[i] + b_gender*(d$gender[i]-1)) )


phi1 <- c()
for (pop in 1:2) {
  expect_pos = 0
  total = 0
  for (a in 1:2){
    for (b in 1:max(Age_range)){
      total = total + D1[pop,a,b];
      expect_pos = expect_pos + D1[pop,a,b] * inv_logit(p_logit_culture[pop] + b_gender*(a-1) + b_age*b);
    }
  }  
  phi1[pop] = expect_pos / total
}







d1 <- list(N = N,
           MA = max(Age_range),
           gender = d$gender[d$soc_id==1], 
           age = d$age[d$soc_id==1],
           outcome = d$outcome[d$soc_id==1]
)

d2 <- list(N = N,
           MA = max(Age_range),
           gender = d$gender[d$soc_id==2], 
           age = d$age[d$soc_id==2],
           outcome = d$outcome[d$soc_id==2]
)


d1$P_same <-a11
d1$P_other <-  a12


d2$P_same <- a12

d2$P_other <-  a11




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



m11 <- stan( model_code  = m2a_MRP_GP_gender_same , data=d1 ,iter = 2000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.95, max_treedepth = 13))  

m12 <- stan( model_code  = m2a_MRP_GP_gender_same , data=d2 ,iter = 2000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.95, max_treedepth = 13))  

m11samp <- extract.samples(m11)
m12samp <- extract.samples(m12)





####                     
###
##
# b) Differences in sampling of ages/genders
##
###
####


#Create population array [pop, sex, age]
D2 <- array(NA,c(2,2,length(Age_range)))

X <- 1:100

age_seq <- seq(0, 100, 1)

#Exponential distribution (similar to many growing populations)
D2[1,1, ] <- dexp(age_seq, 0.04)[Age_range]
D2[1,2, ] <- dexp(age_seq, 0.04)[Age_range]

#Skewed normal distribution (similar to many shrinking populations)
D2[2,1, ] <- dexp(age_seq, 0.04)[Age_range]
D2[2,2, ] <-dexp(age_seq, 0.04)[Age_range]



#Generate data
#Construct dataframe
d <- data.frame(id = 1:(2*N), soc_id = c(rep(1,N), rep(2,N)), age = NA, gender = NA, outcome = NA)

#Simulate ages from above distributions
for (pop_id in 1:2) {
  d$age[d$soc_id==pop_id] <- sample(Age_range, N, replace = TRUE, prob = D2[pop_id,1,])
}

#Simulate Genders

d$gender[d$soc_id==1] <- sample(c(1,2),N, replace = TRUE, prob = c(0.4, 0.6)) 
d$gender[d$soc_id==2] <- sample(c(1,2),N, replace = TRUE, prob = c(0.7, 0.3)) 



a21 <- matrix(0, 2, length(Age_range))
for (i in Age_range) {
  a21[1,which(Age_range == i)] <- length(which(d$age[d$soc_id == 1]==i & d$gender[d$soc_id == 1] == 1))
  a21[2,which(Age_range == i)] <- length(which(d$age[d$soc_id == 1]==i & d$gender[d$soc_id == 1] == 2))
}
a22 <- matrix(0, 2, length(Age_range))
for (i in Age_range) {
  a22[1,which(Age_range == i)] <- length(which(d$age[d$soc_id == 2]==i & d$gender[d$soc_id == 2] == 1))
  a22[2,which(Age_range == i)] <- length(which(d$age[d$soc_id == 2]==i & d$gender[d$soc_id == 2] == 2))
}





#Generate observations

#Effect of "culture" that's independent of demography
p_logit_culture <-c(-3, -1)

#Effects of age and gender
b_age <- 0.02
b_gender <- 2.5


for(i in 1:(2*N)) d$outcome[i] <- rbinom(1, 1, inv_logit(p_logit_culture[d$soc_id[i]] + b_age*d$age[i] + b_gender*(d$gender[i]-1)) )

#Compute expected prosociality in both populations
phi2 <- c()
for (pop in 1:2) {
  expect_pos = 0
  total = 0
  for (a in 1:2){
    for (b in 1:max(Age_range)){
      total = total + D2[pop,a,b];
      expect_pos = expect_pos + D2[pop,a,b] * inv_logit(p_logit_culture[pop] + b_gender*(a-1) + b_age*b);
    }
  }  
  phi2[pop] = expect_pos / total
}






# Prepare lists for stan
d1 <- list(N = N,
           MA = max(Age_range),
           gender = d$gender[d$soc_id==1], 
           age = d$age[d$soc_id==1],
           outcome = d$outcome[d$soc_id==1]
)

d2 <- list(N = N,
           MA = max(Age_range),
           gender = d$gender[d$soc_id==2], 
           age = d$age[d$soc_id==2],
           outcome = d$outcome[d$soc_id==2]
)


#
d1$P_same <-a21
d1$P_other <-  D2[1,,] * 1000000000
mode(d1$P_other) <- 'integer'

d2$P_same <- a22

d2$P_other <-  D2[2,,] * 1000000000

mode(d2$P_other) <- 'integer'







m21 <- stan( model_code  = m2a_MRP_GP_gender_same , data=d1 ,iter = 2000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.95, max_treedepth = 13))  

m22 <- stan( model_code  = m2a_MRP_GP_gender_same , data=d2 ,iter = 2000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.95, max_treedepth = 13))  

m21samp <- extract.samples(m21)
m22samp <- extract.samples(m22)







