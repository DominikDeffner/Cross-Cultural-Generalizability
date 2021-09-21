# Cross-Cultural Generalizability

This repository contains the scripts to reproduce all results and plots in 

***Deffner, D., Rohrer, J. & McElreath, R. (submitted) A Causal Framework for Cross-Cultural Generalizability***

"Example 1 Generalizing description" corresponds to section 3.1. where we use causal thinking to accurately describe and compare the distribution of a trait across societies.

- "DemoStandSim.r" shows simulated data example for demographic standardization (Appendix B in the supplementary material)
- "DemoStandHouse" shows real data example for demographic standardization (section 3.1.2 in the main text)

"Example 2 Generalizing experimental results" corresponds to section 3.3. where we use causal thinking to generaliza and compare causal effects across societies.

- "TransportHouse" shows real data example for the transport of causal effects across populations standardization (section 3.2.2 in the main text)

Each file contains code to simulate/prepare data, run the multilevel regression with poststratification/transport models and produce the plots in the manuscript

"data" contains relevant experimental data from House et al., 2020 (https://www.nature.com/articles/s41562-019-0734-z), 
as well as demographic data from Vanuatu und Berlin used for the data example in section 3.1.