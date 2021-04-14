# Gaussian process post-stratification over different 
# age distributions (with/without gender)


library(rethinking)

###
##
# Model 1: Gaussian process on age
##
###


#Simulate data
N <- 400
MA <- 60    #Max Age
age = round(runif(N, 1,MA))

pred_p <- logistic(1.95*sin(seq(0,6, length.out=MA)))   #Prob of choosing 1 for different ages
plot(pred_p, type = "l", xlab = "Age")
outcome <- rep(NA,N)

#Generate outcome
for(i in 1:N) outcome[i] <- rbinom(1, 1, pred_p[age[i]] )

d<-data.frame(age=age, outcome=outcome)

#Drop data from certain age groups
ages_2_drop <- c(5,6, 13,14,15, 27,28, 37,39, 46:55)
d2 <- d[-which(d$age %in% ages_2_drop),]

#Construct data list
N2 <- nrow(d2)
d <- list(N = N2,
          MA = MA,
          age = d2$age,
          outcome = d2$outcome
)


#Here we specify age distribution in target population that we use for weighting
d$P <- rep(0, MA)
d$P[40:60] <- 100


# Gaussian Process on age
{
m1_MRP_GP <- "
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
  int<lower = 0> P[MA];   // Here we enter data matrix with demographic constitution of target population
}

parameters {
  real alpha;
  
  vector[MA] age_effect;    //Vector for GP age effects
  
  real<lower=0> eta;
  real<lower=0> sigma;
  real<lower=0, upper=1> rho;
}

model {
  vector[N] p;

  alpha ~ normal(0, 3);

  eta ~ exponential(1);
  sigma ~ exponential(1);
  rho ~ beta(10, 1);

  age_effect ~ multi_normal_cholesky( rep_vector(0, MA) , GPL(MA, rho, eta, sigma) );

  for ( i in 1:N ) {
   p[i] = alpha + age_effect[age[i]];
  }

 outcome ~ binomial_logit(1, p);
}
 
generated quantities{
   real expect_pos = 0;
   real<lower = 0, upper = 1> phi;  // This is value for p in the target population 
   int total = 0;
   vector[MA] pred_p;

        for (b in 1:MA){

         total += P[b];
         expect_pos += P[b] * inv_logit(alpha + age_effect[b]);
     }
     
   phi = expect_pos / total;
   
  pred_p = inv_logit(alpha + age_effect);
}

"
}



m1 <- stan( model_code  = m1_MRP_GP , data=d ,iter = 2000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.9, max_treedepth = 13))  

#Extract samples and plot predictions
samp <- extract.samples(m1)
samp_pred_p <- apply(samp$pred_p,2,quantile,probs=c(0.05,0.5,0.95))
Lower <- samp_pred_p[1,]
Median <- samp_pred_p[2,]
Higher <- samp_pred_p[3,]
Observed <- ifelse(c(1:MA) %in% ages_2_drop,"Unobserved Age","Observed Age")

plot(pred_p, type = "l", xlim = c(0,MA), ylim = c(0,1), xlab = "Age", ylab = "p", lwd = 2, lty = 2)
points(Median, col = ifelse(Observed == "Observed Age", "black", "blue"), pch = 16)
segments(x0 = 1:MA, y0 = Lower, x1=  1:MA, y1 = Higher, lwd = 2, col = ifelse(Observed == "Observed Age", "black", "blue"))
legend("topright", c("Observed", "Unobserved"), col = c("black", "blue"), lwd = 4, lty=1, bty="n")





###
##
# Model 2a: Gaussian processes on age, separately for 2 genders (same control parameters)
# Model 2b: Gaussian processes on age, separately for 2 genders (separate control parameters)
##
###


N <- 100
MA <- 60

age = round(runif(N, 1,MA))
gender <- sample(c(1,2), size = N, replace = TRUE)

# Age distribution in target population

pred_p_m <- logistic(1.95*sin(seq(0,6, length.out=MA)) -2)
pred_p_f <- logistic(1.95*sin(seq(0,6, length.out=MA)) + 2)

plot(pred_p_m, type = "l", xlab = "Age", ylim = c(0,1))
lines(pred_p_f, col = "blue", lty = 2)
legend("topright", c("m", "f"), col = c("black", "blue"), lwd = 1, lty=c(1,2), bty="n")

outcome <- rep(NA,N)

for(i in 1:N){
  if (gender[i] == 1){
  outcome[i] <- rbinom(1, 1,  pred_p_m[age[i]] )
  } else {
  outcome[i] <- rbinom(1, 1,  pred_p_f[age[i]] )
  }
}


d<-data.frame(age=age, gender = gender, outcome=outcome)
ages_2_drop <- c(5,6, 13,14,15, 27,28, 37,39, 46:MA)
d2 <- d[-which(d$age %in% ages_2_drop),]

N2 <- nrow(d2)
d <- list(N = N2,
          MA = MA,
          gender = d2$gender, 
          age = d2$age,
          outcome = d2$outcome
)

d$P <- matrix(0, nrow = 2, ncol = MA)
d$P[1, 20:40] <- 100  # Male
d$P[2, 41:60] <- 100  # Female


#2a: Same control parameters

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
  int<lower = 0> P[2,MA];   // Here we enter data matrix with demographic constitution of target population
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
   real<lower = 0, upper = 1> phi;  // This is value for p in the target population 
   int total = 0;
   vector[MA] pred_p_m;
   vector[MA] pred_p_f;

   for (a in 1:2)
        for (b in 1:MA){
         total += P[a,b];
         expect_pos += P[a,b] * inv_logit(alpha[a] + age_effect[a,b]);
     }
     
   phi = expect_pos / total;
   
  pred_p_m = inv_logit(alpha[1] + age_effect[1,]');
  pred_p_f = inv_logit(alpha[2] + age_effect[2,]');
}

"

}



#2b: Separate control parameters
{
  
m2a_MRP_GP_gender_sep <- "
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
  int<lower = 0> P[2,MA];   // Here we enter data matrix with demographic constitution of target population
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

  alpha ~ normal(0, 1);

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
   real<lower = 0, upper = 1> phi;  // This is value for p in the target population 
   int total = 0;
   vector[MA] pred_p_m;
   vector[MA] pred_p_f;

   for (a in 1:2)
        for (b in 1:MA){
         total += P[a,b];
         expect_pos += P[a,b] * inv_logit(alpha[a] + age_effect[a,b]);
     }
     
   phi = expect_pos / total;
   
  pred_p_m = inv_logit(alpha[1] + age_effect[1,]');
  pred_p_f = inv_logit(alpha[2] + age_effect[2,]');
}

"

}

m2a <- stan( model_code  = m2a_MRP_GP_gender_same , data=d ,iter = 2000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.95, max_treedepth = 13))  
m2b <- stan( model_code  = m2a_MRP_GP_gender_sep , data=d ,iter = 2000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.95, max_treedepth = 13))  


par(mfrow = c(1,2))
samp <- extract.samples(m2a)
samp_pred_p <- apply(samp$pred_p_m,2,quantile,probs=c(0.05,0.5,0.95))
Lower <- samp_pred_p[1,]
Median <- samp_pred_p[2,]
Higher <- samp_pred_p[3,]
Observed <- ifelse(c(1:MA) %in% ages_2_drop,"Unobserved Age","Observed Age")

plot(pred_p_m, type = "l", xlim = c(0,MA), ylim = c(0,1), xlab = "Age", ylab = "p", lwd = 2, lty = 2, main= "Male")
points(Median, col = ifelse(Observed == "Observed Age", "black", "blue"), pch = 16)
segments(x0 = 1:MA, y0 = Lower, x1=  1:MA, y1 = Higher, lwd = 2, col = ifelse(Observed == "Observed Age", "black", "blue"))
legend("topright", c("Observed", "Unobserved"), col = c("black", "blue"), lwd = 4, lty=1, bty="n")


samp_pred_p <- apply(samp$pred_p_f,2,quantile,probs=c(0.05,0.5,0.95))
Lower <- samp_pred_p[1,]
Median <- samp_pred_p[2,]
Higher <- samp_pred_p[3,]
Observed <- ifelse(c(1:MA) %in% ages_2_drop,"Unobserved Age","Observed Age")

plot(pred_p_f, type = "l", xlim = c(0,MA), ylim = c(0,1), xlab = "Age", ylab = "p", lwd = 2, lty = 2, main= "Female")
points(Median, col = ifelse(Observed == "Observed Age", "black", "blue"), pch = 16)
segments(x0 = 1:MA, y0 = Lower, x1=  1:MA, y1 = Higher, lwd = 2, col = ifelse(Observed == "Observed Age", "black", "blue"))










####

###




par(mfrow = c(1,2))
samp <- extract.samples(m2b)
samp_pred_p <- apply(samp$pred_p_m,2,quantile,probs=c(0.05,0.5,0.95))
Lower <- samp_pred_p[1,]
Median <- samp_pred_p[2,]
Higher <- samp_pred_p[3,]
Observed <- ifelse(c(1:MA) %in% ages_2_drop,"Unobserved Age","Observed Age")

plot(pred_p_m, type = "l", xlim = c(0,MA), ylim = c(0,1), xlab = "Age", ylab = "p", lwd = 2, lty = 2, main= "Male")
points(Median, col = ifelse(Observed == "Observed Age", "black", "blue"), pch = 16)
segments(x0 = 1:MA, y0 = Lower, x1=  1:MA, y1 = Higher, lwd = 2, col = ifelse(Observed == "Observed Age", "black", "blue"))
legend("topright", c("Observed", "Unobserved"), col = c("black", "blue"), lwd = 4, lty=1, bty="n")


samp_pred_p <- apply(samp$pred_p_f,2,quantile,probs=c(0.05,0.5,0.95))
Lower <- samp_pred_p[1,]
Median <- samp_pred_p[2,]
Higher <- samp_pred_p[3,]
Observed <- ifelse(c(1:MA) %in% ages_2_drop,"Unobserved Age","Observed Age")

plot(pred_p_f, type = "l", xlim = c(0,MA), ylim = c(0,1), xlab = "Age", ylab = "p", lwd = 2, lty = 2, main= "Female")
points(Median, col = ifelse(Observed == "Observed Age", "black", "blue"), pch = 16)
segments(x0 = 1:MA, y0 = Lower, x1=  1:MA, y1 = Higher, lwd = 2, col = ifelse(Observed == "Observed Age", "black", "blue"))

 



