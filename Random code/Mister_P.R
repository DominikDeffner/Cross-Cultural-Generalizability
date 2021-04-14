# Try to do simple MRP
# If we want to compare distribution of trait, we need to account for demographic makeup

setwd("~/GitHub/Cross_Cultural_Generalizability")

data <- read.csv("House_data/Model_1a_1b_1c_data.csv")

data <- subset(data, data$fieldid == 1)

P <- matrix(10, nrow = 4, ncol = 2)
P[,1] <- 10
P[,2] <- 100

P[2,2] <- 0

# plot

d <- list(y = data$T1_ad_choice_1yes,
          ind_id = sapply(1:nrow(data), function (i) which(unique(data$SUBJECT_ID) == data$SUBJECT_ID[i])),
          soc_id = data$fieldid,
          age = data$AGE_in_years,
          gender = as.integer(data$GENDER_1female +1), # 2 female now
          N = nrow(data),
          P = P
)


d$age <- sapply(1:length(d$age), function (i){
  if (d$age[i] <= 30) return(1)
  if (d$age[i] > 30 & d$age[i] <=40 ) return(2)
  if (d$age[i] > 40 & d$age[i] <=50 ) return(3)
  if (d$age[i] > 50) return(4)})



library(rethinking)


### Mister P stan model with gender and discrete age categories
{
MRP_stan <- "

data {
  int N;
  int age[N];
  int gender[N];
  int soc_id[N];
  int y[N];
  int<lower = 0> P[4,2];   // Here we enter data matrix with demographic constitution of target population
}
parameters {
  real alpha;
  
  real<lower = 0> sigma_beta;  
  vector<multiplier = sigma_beta>[4] beta;
  
  real<lower = 0> sigma_gamma;
  vector<multiplier = sigma_gamma>[2] gamma;
}

model {
  y ~ bernoulli_logit(alpha + beta[age] + gamma[gender]);
  alpha ~ normal(0, 2);
  beta ~ normal(0, sigma_beta);
  gamma ~ normal(0, sigma_gamma);
  sigma_beta ~ exponential(1); 
  sigma_gamma ~ exponential(1);
}
generated quantities {
   real expect_pos = 0;
   real<lower = 0, upper = 1> phi;
   int total = 0;
   for (b in 1:4)
     for (c in 1:2){
         total += P[b, c];
         expect_pos += P[b, c] * inv_logit(alpha + beta[b] + gamma[c]);
     }
     
   phi = expect_pos / total;
}

"
}





d$age <- data$AGE_in_years
d$AgeMat <- matrix(nrow = d$N, ncol = d$N)
for (i in 1:d$N) {
  for (j in 1:d$N) {
    d$AgeMat[i,j] <- abs(d$age[i]-d$age[j])
  }
}


d$x_pred <- c(10,50,80)

d$PredMat <- matrix(nrow = length(d$x_pred), ncol = length(d$x_pred))
for (i in 1:length(d$x_pred)) {
  for (j in 1:length(d$x_pred)) {
    d$PredMat[i,j] <- abs(d$x_pred[i]-d$x_pred[j])
  }
}

d$N_pred <- length(d$x_pred)


### Mister P model with Gaussian process on age
{
  MRP_GP <- "
functions{
  
  matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
    int N = dims(x)[1];
    matrix[N, N] K;
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha + delta;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha + delta;
    return K;
  }
}

data {
  int N;
  real age[N];
  int gender[N];
  int soc_id[N];
  int y[N];
  int<lower = 0> P[4,2];   // Here we enter data matrix with demographic constitution of target population
  matrix[N,N] AgeMat;
  int N_pred;
  real x_pred[N_pred];
  matrix[N_pred,N_pred] PredMat;

}
parameters {
  real alpha;
  
  vector[N] k;
  real<lower=0> etasq;
  real<lower=0> rhosq;
}

model {
  vector[N] p;
  matrix[N,N] SIGMA;

  alpha ~ normal(0, 2);

  rhosq ~ exponential( 0.5 );
  etasq ~ exponential( 2 );

  SIGMA = cov_GPL2(AgeMat, etasq, rhosq, 0.01);
  k ~ multi_normal( rep_vector(0,N) , SIGMA );

  for ( i in 1:N ) {
   p[i] = inv_logit( alpha + k[i]);
  }

y ~ binomial(1, p);

}


 
 generated quantities {
   matrix[N_pred,N_pred] SIGMA_pred;
   vector[N_pred] k_pred;
   vector[N_pred] p_pred;

   SIGMA_pred = cov_GPL2(PredMat, etasq, rhosq, 0.01);
   k_pred = multi_normal_rng( rep_vector(0,N_pred) , SIGMA_pred );

  for ( i in 1:N_pred ) {
    p_pred[i] = inv_logit( alpha + k_pred[i]);
  }
}
 
 "
 }





#m <- stan( model_code  = MRP_stan , data=d ,iter = 2000, cores = 1, chains=1, control = list(adapt_delta=0.8, max_treedepth = 10))  

m_GP <- stan( model_code  = MRP_GP , data=d ,iter = 2000, cores = 1, chains=1, control = list(adapt_delta=0.8, max_treedepth = 10))  


