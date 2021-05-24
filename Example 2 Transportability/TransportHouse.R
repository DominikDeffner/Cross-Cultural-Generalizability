
# House data for transport example


library(readr)
library(rethinking)
library(plotrix)

setwd("~/GitHub/Cross_Cultural_Generalizability")

data <- read.csv("House_data/Model_6a_6b_6c_6d_data.csv")

data$condition <- sapply(1:nrow(data), function(x) which(c(data$CONDITION_1_1yes[x],data$CONDITION_2_1yes[x],data$CONDITION_3_1yes[x]) == 1 ))

data <- data[,c("SUBJECT_ID","GENDER_1female","fieldid","AGE_in_years","condition","T1_choice_1yes")]
data$choice <- data$T1_choice_1yes
data$T1_choice_1yes <- NULL
data$age <- round(data$AGE_in_years)
data$gender <- data$GENDER_1female + 1
# Exclude data from "both ok"
data <- subset(data, data$condition != 3)

data_Phoenix <- data[which(data$fieldid == 3),]
data_Pune    <- data[which(data$fieldid == 4),]


Demo <- array(0, dim = c(6 , 12, 2))
for (j in sort(unique(data$fieldid))) {
  for (i in 4:15) {
    for (g in 1:2) {
      Demo[j, which(4:14 == i),g] <- length(which(data$fieldid==j & data$age==i & data$gender== g))
    }
  }
}          



Demo <- matrix(0, 6, 12)
for (j in sort(unique(data$fieldid))) {
  for (i in 4:15) {
      Demo[j, which(4:15 == i)] <- length(which(data$fieldid==j & data$age==i))
    
  }
}          





Demo_Phoenix <- matrix(0, 12, 2)
Demo_Pune    <- matrix(0, 12, 2)

for (i in 4:15) {
  for (g in 1:2) {
    Demo_Phoenix[which(4:15 == i),g] <- length(which(data_Phoenix$age==i & data_Phoenix$gender== g))
    Demo_Pune[which(4:15 == i),g]     <- length(which(data_Pune$age==i & data_Pune$gender == g))
    
  }
}



#Prepare for stan

d_list_Phoenix <- list(N = nrow(data_Phoenix), 
                       MA = 11,
                       age = data_Phoenix$age - 3,
                       condition = data_Phoenix$condition -1,
                       outcome = data_Phoenix$choice, 
                       gender = data_Phoenix$gender, 
                       P_empirical = t(Demo_Phoenix),
                       P_other = t(Demo_Pune))


d_list_Pune <- list(N = nrow(data_Pune), 
                    MA = 11,
                    age = data_Pune$age - 3,
                    condition = data_Pune$condition -1,
                    outcome = data_Pune$choice, 
                    gender = data_Pune$gender, 
                    P_empirical = t(Demo_Pune),
                    P_other = t(Demo_Phoenix))




d_list <- list(N = nrow(data), 
               N_pop = length(unique(data$fieldid)),
               MA = 12,
               pop_id = data$fieldid,
               age = data$age - 3,
               condition = data$condition -1,
               outcome = data$choice, 
               gender = data$gender,
               Demo = Demo)









  
 model <- "
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
  int N_pop;
  int MA;
  int age[N];
  int pop_id[N];
  int condition[N];
  int outcome[N];
  int gender[N];
  int Demo[6, MA];
  }

parameters {
  real alpha[N_pop];
  real b_prime[N_pop];
  
  matrix[N_pop, MA] age_effect;    //Vector for GP age effects
  
  real<lower=0> eta[N_pop];
  real<lower=0> sigma[N_pop];
  real<lower=0, upper=1> rho[N_pop];
}

model {
  vector[N] p;

  alpha ~ normal(0, 1);
  b_prime ~ normal(0, 1);
  eta ~ exponential(2);
  sigma ~ exponential(1);
  rho ~ beta(10, 1);

 for (h in 1:N_pop){
    age_effect[h,] ~ multi_normal_cholesky( rep_vector(0, MA) , GPL(MA, rho[h], eta[h], sigma[h]) );
 }
    

  for ( i in 1:N ) {
   p[i] = alpha[pop_id[i]] + (b_prime[pop_id[i]] + age_effect[pop_id[i],age[i]]) * condition[i];
  }

 outcome ~ binomial_logit(1, p);
}


generated quantities{

real empirical_p[N_pop];
real transport_p[N_pop];
  
 for (h in 1:N_pop){
  real expect_pos = 0;
  int total = 0;
    for (b in 1:MA){
      total += Demo[h,b];
      expect_pos += Demo[h,b] * inv_logit(b_prime[h] + age_effect[h,b]);
    }
    
  empirical_p[h] = expect_pos / total;
 }

 for (h in 1:N_pop){
  real expect_pos = 0;
  int total = 0;
    for (b in 1:MA){
      total += Demo[1,b];
      expect_pos += Demo[1,b] * inv_logit(b_prime[h] + age_effect[h,b]);
    }
    
  transport_p[h] = expect_pos / total;
}

}


"






library(rstan)
m4 <- stan( model_code  = model , data= d_list ,iter = 2000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.9, max_treedepth = 13))  





m_Phoenix <- stan( model_code  = m2 , data= d_list_Phoenix ,iter = 5000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.9, max_treedepth = 13))  




m_Pune <- stan( model_code  = m2 , data= d_list_Pune ,iter = 5000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.9, max_treedepth = 13))  


