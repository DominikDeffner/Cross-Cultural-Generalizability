
# Example 1: Demographic Standardization
# The aim here is to simulate and analyze exemplary data sets generated from different processes
# a) Disparities arise only from real demographic differences among populations (e.g. one is older/more men)
# b) Disparities arise only from differences in sampling procedure 
# c) Disparities arise from both processes

library(rethinking)
library(scales)
library(sn)
library(plotrix)

N <- 500
Age_range <- c(1:90)


####                     
###
##
# a) Real differences
##
###
####


#Create population array [pop, sex, age]
D <- array(NA,c(2,2,length(Age_range)))

X <- 1:100

age_seq <- seq(0, 100, 1)

#Exponential distribution (similar to many growing populations)
D[1,1, ] <- dexp(age_seq, 0.04)[Age_range]
D[1,2, ] <- dexp(age_seq, 0.04)[Age_range]

#Skewed normal distribution (similar to many shrinking populations)
D[2,1, ] <- dsn(X, xi = 30, omega = 50, alpha = 1)[Age_range]
D[2,2, ] <- dsn(X, xi = 30, omega = 50, alpha = 1)[Age_range]

 
 par(mfrow = c(1,2))
 plot(X, dexp(X,0.04), type = "l", xlim = c(0,100), ylab= "Density", xlab = "Age")
 plot(X, dsn(X, xi = 30, omega = 70, alpha = 1), type = "l", xlim = c(0,100), ylab= "Density", xlab = "Age", ylim = c(0,0.01))


#Generate data
#Construct dataframe
d <- data.frame(id = 1:(2*N), soc_id = c(rep(1,N), rep(2,N)), age = NA, gender = NA, outcome = NA)

#Simulate ages from above distributions
for (pop_id in 1:2) {
  d$age[d$soc_id==pop_id] <- sample(Age_range, N, replace = TRUE, prob = D[pop_id,1,])
}

#Simulate Genders

d$gender[d$soc_id==1] <- sample(c(1,2),N, replace = TRUE, prob = c(0.5, 0.5)) 
d$gender[d$soc_id==2] <- sample(c(1,2),N, replace = TRUE, prob = c(0.5, 0.5)) 



a1 <- matrix(0, 2, length(Age_range))
for (i in Age_range) {
     a1[1,which(Age_range == i)] <- length(which(d$age[d$soc_id == 1]==i & d$gender[d$soc_id == 1] == 1))
     a1[2,which(Age_range == i)] <- length(which(d$age[d$soc_id == 1]==i & d$gender[d$soc_id == 1] == 2))
}
a2 <- matrix(0, 2, length(Age_range))
for (i in Age_range) {
  a2[1,which(Age_range == i)] <- length(which(d$age[d$soc_id == 2]==i & d$gender[d$soc_id == 2] == 1))
  a2[2,which(Age_range == i)] <- length(which(d$age[d$soc_id == 2]==i & d$gender[d$soc_id == 2] == 2))
}


par(mfrow = c(2,2), 
    mar = c(1,1,1,0), 
    oma = c(3.3,4,0,1))
labels1 <- matrix("", length(Age_range), 2)
labels1[seq(5,75,10),1] <- seq(5,75,10)
labels2 <- matrix("", length(Age_range), 2)


par(mar=pyramid.plot(D[1,1, ]*100,D[1,2, ]*100,top.labels=c("", "Society 1",""),ppmar=c(2,1,3,1), xlim = c(3,3),labelcex=1.2, unit = "",show.values=F, labels = labels1, lxcol = "indianred", rxcol = "darkgreen",space = 0,gap = 0))
labels <- matrix("", length(Age_range), 2)
legend("topright", c("Male", "Female"), col = c("indianred", "darkgreen"),cex = 1, lty = 1,lwd = 5, bty = "n" )

mtext("Population", side = 2, outer = F, line = 2.5, cex = 1.3)
mtext("Age class", side = 2, outer = F, line = 0.2, cex = 1)
par(mar=pyramid.plot(D[2,1, ]*100,D[2,2, ]*100,top.labels=c("", "Society 2",""),ppmar=c(2,1,3,1), xlim = c(3,3),labelcex=1.2, unit = "",show.values=F, labels = labels2, lxcol = "indianred", rxcol = "darkgreen", space = 0,gap = 0))
mtext("Share of population per age class and gender [%]", side = 1,line = 1.5,at = -3.5, outer = F, cex = 1)


par(mar=pyramid.plot(a1[1,],a1[2,],top.labels=c("", "",""),ppmar=c(2,1,3,1), xlim = c(max(c(a1,a2)),max(c(a1,a2))),labelcex=1.2, unit = "",show.values=F, labels = labels1, lxcol = "indianred", rxcol = "darkgreen",space = 0,gap = 0))
mtext("Sample", side = 2, outer = F, line = 2.5, cex = 1.3)
mtext("Age class", side = 2, outer = F, line = 0.2, cex = 1)
par(mar=pyramid.plot(a2[1,],a2[2,],top.labels=c("", "",""),ppmar=c(2,1,3,1), xlim = c(max(c(a1,a2)),max(c(a1,a2))),labelcex=1.2, unit = "",show.values=F, labels = labels2, lxcol = "indianred", rxcol = "darkgreen",space = 0,gap = 0))

mtext("Number of individuals per age class and gender", side = 1,line = 1.5,at = -15, outer = F, cex = 1)














p_logit_culture <-c(-1, -1)

b_age <- 0.04
b_gender <- 0

#Generate observations
stand_age<- c()
for (i in 1:(2*N)) stand_age[i] <- (d$age[i] - mean(d$age[d$soc_id==1]))/sd(d$age[d$soc_id==1])

for(i in 1:(2*N)) d$outcome[i] <- rbinom(1, 1, inv_logit(p_logit_culture[d$soc_id[i]] + b_age*d$age[i] + b_gender*(d$gender[i]-1)) )


phi <- c()
for (pop in 1:2) {
  expect_pos = 0
  total = 0
  for (a in 1:2){
    for (b in 1:max(Age_range)){
      total = total + D[pop,a,b];
      expect_pos = expect_pos + D[pop,a,b] * inv_logit(p_logit_culture[pop] + b_gender*(a-1) + b_age*b);
    }
  }  
  phi[pop] = expect_pos / total
}







d1 <- list(N = N,
           MA = max(Age_range),
           gender = d$gender[d$soc_id==1], 
           age = d$age[d$soc_id==1],
           outcome = d$outcome[d$soc_id==1]
)

d2 <- list(N = N,
           MA = max(Age_range),
           gender = d$gender[d$soc_id==2], 
           age = d$age[d$soc_id==2],
           outcome = d$outcome[d$soc_id==2]
)


d1$P_same <-a1
d1$P_other <-  a2
d1$P_other[1,] <- 100
d1$P_other[2,] <- 100

d2$P_same <- a2

d2$P_other <-  a2

d2$P_other[1,] <- 100
d2$P_other[2,] <- 100
















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
  int<lower = 0> P_same[2,MA];   // Here we enter data matrix with demographic constitution of target population
  int<lower = 0> P_other[2,MA];   // Here we enter data matrix with demographic constitution of target population

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
   real<lower = 0, upper = 1> phi;  // This is value for p in the source population 
   real<lower = 0, upper = 1> psi;  // This is value for p in the target population 

   int total = 0;
   vector[MA] pred_p_m;
   vector[MA] pred_p_f;

   for (a in 1:2)
        for (b in 1:MA){
         total += P_same[a,b];
         expect_pos += P_same[a,b] * inv_logit(alpha[a] + age_effect[a,b]);
     }
     
   phi = expect_pos / total;
   
    total = 0;
    expect_pos = 0;
    
       for (a in 1:2)
         for (b in 1:MA){
          total += P_other[a,b];
          expect_pos += P_other[a,b] * inv_logit(alpha[a] + age_effect[a,b]);
      }
      
    psi = expect_pos / total;
   
  pred_p_m = inv_logit(alpha[1] + age_effect[1,]');
  pred_p_f = inv_logit(alpha[2] + age_effect[2,]');
}

"

}



m11 <- stan( model_code  = m2a_MRP_GP_gender_same , data=d1 ,iter = 5000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.95, max_treedepth = 13))  

m12 <- stan( model_code  = m2a_MRP_GP_gender_same , data=d2 ,iter = 5000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.95, max_treedepth = 13))  

m11samp <- extract.samples(m11)
m12samp <- extract.samples(m12)


#Age Curves
{
par(mfrow = c(2,2))
samp <- extract.samples(m11)
samp_pred_p <- apply(samp$pred_p_m,2,quantile,probs=c(0.05,0.5,0.95))
Lower <- samp_pred_p[1,]
Median <- samp_pred_p[2,]
Higher <- samp_pred_p[3,]
Observed <- ifelse(c(1:d2$MA) %in% unique(d1$age),"Observed Age","Unobserved Age")
plot(NA, type = "l", xlim = c(0,d2$MA), ylim = c(0,1), xlab = "Age", ylab = "p", lwd = 2, lty = 2, main= "Male")
points(Median, col = ifelse(Observed == "Observed Age", "black", "blue"), pch = 16)
segments(x0 = 1:d2$MA, y0 = Lower, x1=  1:d2$MA, y1 = Higher, lwd = 2, col = ifelse(Observed == "Observed Age", "black", "blue"))
legend("topright", c("Observed", "Unobserved"), col = c("black", "blue"), lwd = 4, lty=1, bty="n")
samp_pred_p <- apply(samp$pred_p_f,2,quantile,probs=c(0.05,0.5,0.95))
Lower <- samp_pred_p[1,]
Median <- samp_pred_p[2,]
Higher <- samp_pred_p[3,]
Observed <- ifelse(c(1:d1$MA) %in% unique(d1$age),"Observed Age","Unobserved Age")

plot(NA, type = "l", xlim = c(0,d2$MA), ylim = c(0,1), xlab = "Age", ylab = "p", lwd = 2, lty = 2, main= "Female")
points(Median, col = ifelse(Observed == "Observed Age", "black", "blue"), pch = 16)
segments(x0 = 1:d2$MA, y0 = Lower, x1=  1:d2$MA, y1 = Higher, lwd = 2, col = ifelse(Observed == "Observed Age", "black", "blue"))

samp <- extract.samples(m12)
samp_pred_p <- apply(samp$pred_p_m,2,quantile,probs=c(0.05,0.5,0.95))
Lower <- samp_pred_p[1,]
Median <- samp_pred_p[2,]
Higher <- samp_pred_p[3,]
Observed <- ifelse(c(1:d2$MA) %in% unique(d1$age),"Observed Age","Unobserved Age")

plot(NA, type = "l", xlim = c(0,d2$MA), ylim = c(0,1), xlab = "Age", ylab = "p", lwd = 2, lty = 2, main= "Male")
points(Median, col = ifelse(Observed == "Observed Age", "black", "blue"), pch = 16)
segments(x0 = 1:d2$MA, y0 = Lower, x1=  1:d2$MA, y1 = Higher, lwd = 2, col = ifelse(Observed == "Observed Age", "black", "blue"))
legend("topright", c("Observed", "Unobserved"), col = c("black", "blue"), lwd = 4, lty=1, bty="n")


samp_pred_p <- apply(samp$pred_p_f,2,quantile,probs=c(0.05,0.5,0.95))
Lower <- samp_pred_p[1,]
Median <- samp_pred_p[2,]
Higher <- samp_pred_p[3,]
Observed <- ifelse(c(1:d1$MA) %in% unique(d1$age),"Observed Age","Unobserved Age")

plot(NA, type = "l", xlim = c(0,d2$MA), ylim = c(0,1), xlab = "Age", ylab = "p", lwd = 2, lty = 2, main= "Female")
points(Median, col = ifelse(Observed == "Observed Age", "black", "blue"), pch = 16)
segments(x0 = 1:d2$MA, y0 = Lower, x1=  1:d2$MA, y1 = Higher, lwd = 2, col = ifelse(Observed == "Observed Age", "black", "blue"))
}



graphics.off()
png("Demostand_nodiff.png", res = 600, height = 22, width = 18, units = "cm")


par(mfrow = c(3,2), 
    mar = c(3,1,3,2), 
    oma = c(1,4,0,0.1))
labels1 <- matrix("", length(Age_range), 2)
labels1[seq(5,75,10),1] <- seq(5,75,10)
labels2 <- matrix("", length(Age_range), 2)


par(mar=pyramid.plot(D[1,1, ]*100,D[1,2, ]*100,top.labels=c("", "Population I",""),ppmar=c(2,1,3,1), xlim = c(5,5),labelcex=1.2, unit = "",show.values=F, labels = labels1, lxcol = "indianred", rxcol = "darkgreen",space = 0,gap = 0))
labels <- matrix("", length(Age_range), 2)
legend("topright", c("Male", "Female"), col = c("indianred", "darkgreen"),cex = 1, lty = 1,lwd = 5, bty = "n" )

mtext("Population", side = 2, outer = F, line = 3.5, cex = 1.3)
mtext("Age class", side = 2, outer = F, line = 0.2, cex = 1)
par(mar=pyramid.plot(D[2,1, ]*100,D[2,2, ]*100,top.labels=c("", "Population II",""),ppmar=c(2,1,3,1), xlim = c(3,3),labelcex=1.2, unit = "",show.values=F, labels = labels2, lxcol = "indianred", rxcol = "darkgreen", space = 0,gap = 0))
mtext("Share of population per age class and gender [%]", side = 1,line = 4,at = -5.5, outer = F, cex = 1)


par(mar=pyramid.plot(a1[1,],a1[2,],top.labels=c("", "",""),ppmar=c(2,1,3,1), xlim = c(max(c(a1,a2)+3),max(c(a1,a2))+3),labelcex=1.2, unit = "",show.values=F, labels = labels1, lxcol = "indianred", rxcol = "darkgreen",space = 0,gap = 0))
mtext("Sample", side = 2, outer = F, line = 3.5, cex = 1.3)
mtext("Age class", side = 2, outer = F, line = 0.2, cex = 1)
par(mar=pyramid.plot(a2[1,],a2[2,],top.labels=c("", "",""),ppmar=c(2,1,3,1), xlim = c(max(c(a1,a2)+3),max(c(a1,a2))+3),labelcex=1.2, unit = "",show.values=F, labels = labels2, lxcol = "indianred", rxcol = "darkgreen",space = 0,gap = 0))

mtext("Number of individuals per age class and gender", side = 1,line = 4,at = -18, outer = F, cex = 1)



dens <- density(m11samp$phi)
x1 <- min(which(dens$x >= quantile(m11samp$phi, 0)))  
x2 <- max(which(dens$x <  quantile(m11samp$phi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[4],alpha = 0.9), border = NA))
abline(v = phi[1], lty = 2, lwd = 2)
par(new=TRUE)
dens <- density(m11samp$psi)
x1 <- min(which(dens$x >= quantile(m11samp$psi, 0)))  
x2 <- max(which(dens$x <  quantile(m11samp$psi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[4],alpha = 0.3), border = NA))
mtext("Density", side = 2,line = 3, outer = F, cex = 1)



dens <- density(m12samp$phi)
x1 <- min(which(dens$x >= quantile(m12samp$phi, 0)))  
x2 <- max(which(dens$x <  quantile(m12samp$phi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[5],alpha = 0.9), border = NA))
abline(v = phi[2], lty = 2, lwd = 2)

par(new=TRUE)

dens <- density(m12samp$psi)
x1 <- min(which(dens$x >= quantile(m12samp$psi, 0)))
x2 <- max(which(dens$x <  quantile(m12samp$psi, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,30), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[5],alpha = 0.3), border = NA))

mtext("Probability of choosing prosocial option", side = 1,line = 3,at=-0.1, outer = F, cex = 1)

legend("topleft", title = "Population", c("I","I (adjusted)", "II", "II (adjusted)"), col = c(alpha(col.pal[4],alpha = 0.9),alpha(col.pal[4],alpha = 0.3),alpha(col.pal[5],alpha = 0.9),alpha(col.pal[5],alpha = 0.3)), lwd = 6, bty="n")

dev.off()



















####                     
###
##
# b) Differences in samping of ages
##
###
####


#Create population array [pop, sex, age]
D <- array(NA,c(2,2,length(Age_range)))


age_seq <- seq(0, 100, 1)

#Exponential distribution for both populations

D[1,1, ] <- dexp(age_seq, 0.03)[Age_range]
D[1,2, ] <- dexp(age_seq, 0.03)[Age_range]
D[2,1, ] <- dexp(age_seq, 0.03)[Age_range]
D[2,2, ] <- dexp(age_seq, 0.03)[Age_range]


#Generate data

d <- data.frame(id = 1:(2*N), soc_id = c(rep(1,N), rep(2,N)), age = NA, gender = NA, outcome = NA)

#Simulate ages from above distributions
for (pop_id in 1:2) {
  d$age[d$soc_id==pop_id] <- sample(Age_range, N, replace = TRUE, prob = D[pop_id,1,])
}

#Simulate Genders

d$gender[d$soc_id==1] <- sample(c(1,2),N, replace = TRUE, prob = c(0.5, 0.5)) 
d$gender[d$soc_id==2] <- sample(c(1,2),N, replace = TRUE, prob = c(0.5, 0.5)) 

b_age <- 2
b_gender <- 5

#Generate observations

for(i in 1:(2*N)) d$outcome[i] <- rbinom(1, 1, inv_logit(b_age*standardize(d$age)[i] + b_gender * (d$gender[i]-1)  ) )



d1 <- list(N = N,
           MA = max(d$age[d$soc_id==1]),
           gender = d$gender[d$soc_id==1], 
           age = d$age[d$soc_id==1],
           outcome = d$outcome[d$soc_id==1]
)

d1$P <- matrix(0, nrow = 2, ncol = d1$MA)
d1$P[1, 10:30] <- 1000  # Male
d1$P[2, 41:60] <- 5  # Female















