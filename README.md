# Cross-Cultural Generalizability

This repository contains the scripts to reproduce all analyses and figures in 

****Deffner, D., Rohrer, J. & McElreath, R. (submitted) A Causal Framework for Cross-Cultural Generalizability****

The preprint can be found at https://psyarxiv.com/fqukp

"Example 1 Generalizing description" corresponds to section 3.1. where we use causal thinking to accurately describe and compare the distribution of a trait across societies.

- "DemoStandSim.r" shows simulated data example for demographic standardization (Appendix B in the supplementary material)
- "DemoStandHouse" shows real data example for demographic standardization (section 3.1.2. in the main text)
- "model_empirical.stan" provides code for a simple Bernoulli model that provides unbiased or "empirical" estimate in each population. 
- "model_MRpoststratification.stan" uses Gaussian processes to compute age and gender specific estimates and poststratifies to population from which sample was taken and to the other population. Remember that arbitrary target populations are possible and populations can vary in any number of background factors.


"Example 2 Generalizing experimental results" corresponds to section 3.2. where we use the causal framework to generalize and compare causal effects across societies.

- "TransportHouse.r" shows real data example for the transport of causal effects across populations (section 3.2.2. in the main text)
- "model_basic.stan" is a simple fixed effects model that provides unbiased or "empirical" estimates of the causal effect in each population.
- "model_transport.stan" uses Gaussian processes to compute age-specific causal effects for each population. It then uses these strata-specific effects to transport effects to any arbitrary target population.

R files contain code to simulate/prepare data, run the multilevel regression with poststratification/transport models and produce the plots in the manuscript

"data" contains relevant experimental data from House et al., 2020 (https://www.nature.com/articles/s41562-019-0734-z), 
as well as demographic data from Vanuatu und Berlin used for the data example in section 3.1.

***Software requirements***
The code was written in R 4.0.3. Statistical models are fit using the Stan MCMC engine via the rstan package (2.21.2), which requires a C++ compiler. Installation        instructions are available at https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started. See also the Stan user guide at https://mc-stan.org/users/documentation. The rethinking package (2.12) is required to process fitted model outputs (installation instructions at http://xcelab.net/rm/software/).
