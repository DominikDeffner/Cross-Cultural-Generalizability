
// TRANSPORT OF CAUSAL EFFECTS ACROSS POPULATIONS

// This model uses Gaussian processes to compute age-specific causal effects for each population.
// It then uses these strata-specific effects to transport effects to any arbitrary target population


functions{
// Define function for Gaussian process; this computes how the covariance between different ages is expected to change as the distance increases
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

// Define the observed variables we feed into the model as data
data {
  int N;                     // Number of unique choices/participants
  int N_pop;                 // Number of field sites
  int MA;                    // Maximum age
  int age[N];                // Age
  int pop_id[N];             // Numerical field id
  int condition[N];          // Dummy-coded conditions
  int outcome[N];            // Choice (1 means prosocial)
  int gender[N];             // Gender (1 means male, 2 means female)
  int Demo[N_pop, MA];       // Demography of samples
  int Ref;                   // Reference population for transport
}

// Define the unobserved variables (parameters) that we estimate from the data
parameters {
  real alpha[N_pop];            //Population-specific intercept
  real b_prime[N_pop];          //Population-specific influence of norm primes on the logit scale
  matrix[N_pop, MA] age_effect; //Matrix to hold Gaussian process age effects for each population

  // Here we define the Control parameters (separately for populations) for the Gaussian processes;
  // they determine how covariance changes with increasing distance in age
  real<lower=0> eta[N_pop];
  real<lower=0> sigma[N_pop];
  real<lower=0, upper=1> rho[N_pop];
}

model {
  vector[N] p;

  // Define priors for parameters
  alpha ~ normal(0, 2);
  b_prime ~ normal(0, 2);
  eta ~ exponential(2);
  sigma ~ exponential(1);
  rho ~ beta(10, 1);

  //We compute age-specific offsets for each population
  for (h in 1:N_pop){
    age_effect[h,] ~ multi_normal_cholesky( rep_vector(0, MA) , GPL(MA, rho[h], eta[h], sigma[h]) );
  }

  //This is the linear model: Choice probabilities are composed of (site-specific) intercept and age-specific effect of norm prime condition
  for ( i in 1:N ) {
   p[i] = alpha[pop_id[i]] + (b_prime[pop_id[i]] + age_effect[pop_id[i],age[i]]) * condition[i];
  }

  //Finally, we need a likelihood function for the observed outcomes
  outcome ~ binomial_logit(1, p);
}

// We use the generated quantities section to compute causal effect and adjust it for demography of target population ("transport").
// Causal effects are defined as differences on the outcome scale (here: choice probabilities), so we first convert to outcome scale and then compute causal contrast

generated quantities{
real transport_p[N_pop];
 for (h in 1:N_pop){
   real total = 0;
    for ( i in 1:MA){
      total+= Demo[Ref,i] * (inv_logit(alpha[h]) - inv_logit(alpha[h] + (b_prime[h] + age_effect[h,i])) );
    }
   transport_p[h] = total / sum(Demo[Ref,]);
 }
}
