# Try to do simple MRP
# If we want to compare distribution of trait, we need to account for demographic makeup

setwd("~/GitHub/Cross_Cultural_Generalizability")

data <- read.csv("House_data/Model_1a_1b_1c_data.csv")

data <- subset(data, data$fieldid == 1)

P <- matrix(10, nrow = 4, ncol = 2)
P[,1] <- 10
P[,2] <- 100
# plot

d <- list(y = data$T1_ad_choice_1yes,
          ind_id = sapply(1:nrow(data), function (i) which(unique(data$SUBJECT_ID) == data$SUBJECT_ID[i])),
          soc_id = data$fieldid,
          age = data$AGE_in_years,
          gender = data$GENDER_1female +1, # 2 female now
          N = nrow(data),
          P = P
)


d$age <- sapply(1:length(d$age), function (i){
  if (d$age[i] <= 30) return(1)
  if (d$age[i] > 30 & d$age[i] <=40 ) return(2)
  if (d$age[i] > 40 & d$age[i] <=50 ) return(3)
  if (d$age[i] > 50) return(4)})



library(rethinking)

MRP_stan <- "

data {
  int N;
  int age[N];
  int gender[N];
  int soc_id[N];
  int y[N];
  int<lower = 0> P[4,2];
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


m <- stan( model_code  = MRP_stan , data=d ,iter = 2000, cores = 1, chains=1, control = list(adapt_delta=0.8, max_treedepth = 10))  


