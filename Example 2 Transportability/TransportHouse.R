
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


Demo <- matrix(0, 6, 12)
for (j in sort(unique(data$fieldid))) {
  for (i in 4:15) {
      Demo[j, which(4:15 == i)] <- length(which(data$fieldid==j & data$age==i))
  }
}          



#Prepare for stan

d_list <- list(N = nrow(data), 
               N_pop = length(unique(data$fieldid)),
               MA = 12,
               pop_id = data$fieldid,
               age = data$age - 3,
               condition = data$condition-1,
               outcome = data$choice, 
               gender = data$gender,
               Demo = Demo,
               Ref = 6)



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
  int Demo[N_pop, MA];
  int Ref;
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
  b_prime ~ normal(0, 3);
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
   real empirical_total = 0;
   real transport_total = 0;

    for ( i in 1:MA){
      empirical_total+= Demo[h,i] * (inv_logit(alpha[h]) - inv_logit(alpha[h] + (b_prime[h] + age_effect[h,i])) );
      transport_total+= Demo[Ref,i] * (inv_logit(alpha[h]) - inv_logit(alpha[h] + (b_prime[h] + age_effect[h,i])) );
    }
    empirical_p[h] = empirical_total / sum(Demo[h,]);

    transport_p[h] = transport_total / sum(Demo[Ref,]);
 }



}


"



model_basic <- "

data {
  int N;
  int N_pop;
  int MA;
  int age[N];
  int pop_id[N];
  int condition[N];
  int outcome[N];
  int gender[N];
  int Demo[N_pop, MA];
  int Ref;
  }

parameters {
  real alpha[N_pop];
  real b_prime[N_pop];
}

model {
  vector[N] p;

  alpha ~ normal(0, 1);
  b_prime ~ normal(0, 3);

  for ( i in 1:N ) {
   p[i] = alpha[pop_id[i]] + b_prime[pop_id[i]] * condition[i];
  }

 outcome ~ binomial_logit(1, p);
}

generated quantities{

real empirical_p[N_pop];
 for (h in 1:N_pop){
      empirical_p[h] =  inv_logit(alpha[h]) - inv_logit(alpha[h] + b_prime[h]);
 }
 
}
 
"



model_test <- "

data {
  int N;
  int N_pop;
  int MA;
  int age[N];
  int pop_id[N];
  int condition[N];
  int outcome[N];
  int gender[N];
  int Demo[N_pop, MA];
  int Ref;
  }

parameters {
  real alpha[N_pop];
  real b_prime[N_pop];
  real b_age[N_pop];
  real b_prime_age[N_pop];

}

model {
  vector[N] p;

  alpha ~ normal(0, 3);
  b_prime ~ normal(0, 3);
  b_age ~ normal(0, 3);
  b_prime_age ~ normal(0, 3);

  for ( i in 1:N ) {
   p[i] = alpha[pop_id[i]] + b_prime[pop_id[i]] * condition[i] + b_age[pop_id[i]]*age[i] + b_prime_age[pop_id[i]]* condition[i]*age[i] ;
  }

 outcome ~ binomial_logit(1, p);
}

generated quantities{

real empirical_p[N_pop];
 for (h in 1:N_pop){
      empirical_p[h] =  inv_logit(alpha[h]) - inv_logit(alpha[h] + b_prime[h]);
 }
 
}
 
"


m_test<- stan( model_code  = model_test , data= d_list ,iter = 5000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.9, max_treedepth = 15))  


s <- extract.samples(m_test)

curve(mean(s$alpha[,i]) + mean(s$b_prime[,i] ) )


library(rstan)

m <- list()

# for (i in 1:6) {
#   d_list$Ref <- i
#   m <- stan( model_code  = model , data= d_list ,iter = 5000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.95, max_treedepth = 15))  
# }

d_list$Ref <- 6

m_empirical <- stan( model_code  = model_basic , data= d_list ,iter = 1000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.9, max_treedepth = 15))  
m_transport <- stan( model_code  = model , data= d_list ,iter = 1000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.9, max_treedepth = 15))  


s_basic <- extract.samples(m_empirical)
s <- extract.samples(m_transport)


library(scales)
#color stuff
require(RColorBrewer)#load package
col.pal <- brewer.pal(8, "Dark2") #create a pallette which you loop over for corresponding values
seqoverall <- seq

Society <- c("Berlin (GER)","La Plata (ARG)","Phoenix (USA)", "Pune (IND)", "Shuar (ECU)", "Wichí (ARG)")




graphics.off()
png("Transport.png", res = 900, height = 12, width = 16, units = "cm")

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

dev.off()



