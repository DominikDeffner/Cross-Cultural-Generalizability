# Try to do simple Gaussian process MRP with simulated data


library(rethinking)
set.seed(1)
#Simulate data
N <- 600
MA <- 70
age = round(runif(N, 1,MA))
gender = sample(c(1,2), N, replace = TRUE) 

pred_p <- logistic(3*sin(seq(3,6, length.out=MA)))
outcome <- rep(NA,N)

for(i in 1:N)
outcome[i] <- rbinom(1, 1, pred_p[age[i]])

d<-data.frame(age=age, gender= gender, outcome=outcome)
ages_2_drop <- 10:50
d2 <- d[-which(d$age %in% ages_2_drop),]

N2 <- nrow(d2)
d <- list(N = N2,
          MA = MA,
          age = d2$age,
          gender = d2$gender,
          outcome = d2$outcome
)

d$P <- matrix(0, nrow = 2, ncol = MA)
d$P[1 , 30:50] <- 100


MRP_GP <- "
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
  int gender[N];
  int outcome[N];
  int<lower = 0> P[2,MA];   // Here we enter data matrix with demographic constitution of target population
}

parameters {
  real alpha;
  
    
  vector[2] gender_effect;
  
  vector[MA] age_effect;
  
  real<lower=0> eta;
  real<lower=0> sigma;
  real<lower=0, upper=1> rho;
}

model {
  vector[N] p;

  alpha ~ normal(0, 3);
  to_vector(gender_effect) ~ normal(0,1);

  eta ~ exponential(1);
  sigma ~ exponential(1);
  rho ~ beta(5, 1);

  age_effect ~ multi_normal_cholesky( rep_vector(0, MA) , GPL(MA, rho, eta, sigma) );

  for ( i in 1:N ) {
   p[i] = alpha + age_effect[age[i]] + gender_effect[gender[i]];
  }

 outcome ~ binomial_logit(1, p);
}
 
generated quantities{
   real expect_pos = 0;
   real<lower = 0, upper = 1> phi;
   int total = 0;
   vector[MA] pred_p;

   for (a in 1:2)
        for (b in 1:MA){

         total += P[a,b];
         expect_pos += P[a,b] * inv_logit(alpha + gender_effect[a] + age_effect[b]);
     }
     
   phi = expect_pos / total;
   
  pred_p = inv_logit(alpha + age_effect);
}

"


m_GP <- stan( model_code  = MRP_GP , data=d ,iter = 2000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.9, max_treedepth = 13))  


m2d<-rstan::extract(m_GP,pars="pred_p")
sample_eff<-apply(m2d$pred_p,2,quantile,probs=c(0.05,0.5,0.95))
df_dist<-data.frame(Age=c(1:MA),
                      LI=sample_eff[1,],
                      Median=sample_eff[2,],
                      HI=sample_eff[3,],
                      Observed=ifelse(c(1:MA) %in% ages_2_drop,"Unobserved Age","Observed Age"))
                      
df_add <- data.frame(Age=c(1:MA), Median=pred_p)

ggplot() +
 geom_point(data=df_dist,aes(x=Age,y=Median, color=Observed)) +
 geom_linerange(data=df_dist,aes(x=Age,ymin=LI,ymax=HI, color=Observed)) +
 geom_line(data = df_add, aes(x=Age,y=Median),color="royalblue") + 
 scale_color_manual(values = c("Unobserved Age" = 'darkred', "Observed Age"='black')) +
 labs(y="Prob(Outcome=1)") + theme(strip.text.x = element_text(size=14,face="bold"),axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

