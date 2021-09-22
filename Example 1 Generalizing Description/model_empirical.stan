data {
  int N;
  int outcome[N];
}

parameters {
  real p;
}

model {
 outcome ~ bernoulli_logit(p);
}
