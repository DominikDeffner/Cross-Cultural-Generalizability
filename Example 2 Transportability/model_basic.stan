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
  alpha ~ normal(0, 2);
  b_prime ~ normal(0, 2);
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
