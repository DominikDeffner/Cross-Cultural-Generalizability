# Try to do simple MRP with real data
# If we want to compare distribution of trait, we need to account for demographic makeup

library(rethinking)

setwd("~/GitHub/Cross_Cultural_Generalizability")

data_adult <- read.csv("House_data/Model_1a_1b_1c_data.csv")
data_adult <- data_adult[,c("SUBJECT_ID","GENDER_1female","fieldid", "AGE_in_years","T1_ad_choice_1yes")]
data_adult$choice <- data_adult$T1_ad_choice_1yes
data_adult$T1_ad_choice_1yes <- NULL


data_children <- read.csv("House_data/Model_4a_4b_4c_4d_data.csv")
data_children <- data_children[,c("SUBJECT_ID","GENDER_1female","fieldid", "AGE_in_years","T1_choice_1yes")]
data_children$choice <- data_children$T1_choice_1yes
data_children$T1_choice_1yes <- NULL


data_comb <- rbind(data_children, data_adult)
data_comb <- data_comb[with(data_comb, order(data_comb$fieldid, data_comb$AGE_in_years)),]


data_comb <- data_comb[which(data_comb$fieldid == 1),]
d  <- list(y = data_comb$choice,
           ind_id = sapply(1:nrow(data_comb), function (i) which(unique(data_comb$SUBJECT_ID) == data_comb$SUBJECT_ID[i])),
           soc_id = data_comb$fieldid,
           age = round(data_comb$AGE_in_years),
           gender = as.integer(data_comb$GENDER_1female +1), # 2 female now
           N = nrow(data_comb),
           N_age = 9,
           MA = max(round(data_comb$AGE_in_years))
           )

# Categorize ppl in x age classes


d$age <- sapply(1:length(d$age), function (i){
  if (d$age[i] <= 5) return(1)
  if (d$age[i] > 5 & d$age[i] <=10 ) return(2)
  if (d$age[i] > 10 & d$age[i] <=15 ) return(3)
  if (d$age[i] > 15 & d$age[i] <=25 ) return(4)
  if (d$age[i] > 25 & d$age[i] <=35 ) return(5)
  if (d$age[i] > 35 & d$age[i] <=45 ) return(6)
  if (d$age[i] > 45 & d$age[i] <=55 ) return(7)
  if (d$age[i] > 55 & d$age[i] <=65 ) return(8)
  if (d$age[i] > 65) return(9)
})




# Here we construct a data matrix for number of individuals in demographic classes in target population

P <- matrix(9, nrow = d$N_age, ncol = 2)
P[,1] <- 10
P[,2] <- 100

P[2,2] <- 0

d$P <- P


### Mister P stan model with gender and discrete age categories
{
MRP_stan <- "

data {
  int N;
  int N_age;
  int age[N];
  int<lower = 1, upper = 2> gender[N];
  int soc_id[N];
  int y[N];
  int<lower = 0> P[N_age,2];   // Here we enter data matrix with demographic constitution of target population
}
parameters {
  vector[2] alpha;
  
  real<lower = 0> sigma_age;  
  vector<multiplier = sigma_age>[N_age] age_effect;
  
}

model {
  y ~ bernoulli_logit(alpha[gender] + age_effect[age]);
  alpha ~ normal(0, 2);
  age_effect ~ normal(0, sigma_age);
  sigma_age ~ exponential(1); 

}
generated quantities {
   real expect_pos = 0;
   real<lower = 0, upper = 1> phi;
   int total = 0;
   for (a in 1:N_age){
     for (b in 1:2){
         total += P[a, b];
         expect_pos += P[a, b] * inv_logit(alpha[b] + age_effect[a]);
     }
    }
   phi = expect_pos / total;
}

"
}


m <- stan( model_code  = MRP_stan , data=d ,iter = 2000, cores = 1, chains=1, control = list(adapt_delta=0.8, max_treedepth = 10))  

precis(m,2)
