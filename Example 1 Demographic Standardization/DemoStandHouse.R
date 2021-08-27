
#Generalizing description: Cross-cultural comparisons and demographic standardization
#Real data example based on House et al., 2020
#The script prepares the data, runs the multilevel regression with poststratification analysis and creates plot in Fig.4 in the manuscript

#Load some packages and set working directory to main folder
library(readr)
library(rethinking)
library(plotrix)
library(scales)
require(RColorBrewer)
setwd("~/GitHub/Cross_Cultural_Generalizability")

#Population distributions (age and gender) of Vanuatu and Berlin
Pop_Berlin <- read.csv("Example 1 Demographic Standardization/Berlin-2020.csv")
Pop_Vanuatu <- read.csv("Example 1 Demographic Standardization/Vanuatu-2019.csv")
Pop_Vanuatu <- Pop_Vanuatu[-nrow(Pop_Vanuatu),] #Delete empty cells (100+)
Pop_Berlin <- as.matrix(Pop_Berlin[,c(2,3)])
Pop_Vanuatu <- as.matrix(Pop_Vanuatu[,c(2,3)])

#Change to proportion of different cells in populations
Pop_Berlin_raw <- Pop_Berlin
Pop_Berlin <- (Pop_Berlin/sum(Pop_Berlin))*100

Pop_Vanuatu_raw <- Pop_Vanuatu
Pop_Vanuatu <- (Pop_Vanuatu/sum(Pop_Vanuatu))*100


#Load data from adults and children from House et al.,2020
data_adult <- read.csv("House_data/Model_1a_1b_1c_data.csv")
data_adult <- data_adult[,c("SUBJECT_ID","GENDER_1female","fieldid", "AGE_in_years","T1_ad_choice_1yes")]
data_adult$choice <- data_adult$T1_ad_choice_1yes
data_adult$T1_ad_choice_1yes <- NULL

data_children <- read.csv("House_data/Model_4a_4b_4c_4d_data.csv")
data_children <- data_children[,c("SUBJECT_ID","GENDER_1female","fieldid", "AGE_in_years","T1_choice_1yes")]
data_children$choice <- data_children$T1_choice_1yes
data_children$T1_choice_1yes <- NULL


#Combine both datasets
data_comb <- rbind(data_children, data_adult)
data_comb <- data_comb[with(data_comb, order(data_comb$fieldid, data_comb$AGE_in_years)),]
data_comb$gender <- data_comb$GENDER_1female +1
MA <- max(round(data_comb$AGE_in_years))

#Extract data from Berlin and Tanna, Vanuatu
d_Berlin <- data_comb[which(data_comb$fieldid == 1),]
d_Vanuatu <- data_comb[which(data_comb$fieldid == 7),]

#Categorize data into 20 age categories spanning 5 years each
d_Berlin$AGE_binned <- sapply(1:nrow(d_Berlin), function(i) which.max( 1/(seq(5, 100, 5) - d_Berlin$AGE_in_years[i]) ) )
d_Vanuatu$AGE_binned <- sapply(1:nrow(d_Vanuatu), function(i) which.max( 1/(seq(5, 100, 5) - d_Vanuatu$AGE_in_years[i]) ) )


#Sample distributions (age and gender) of Vanuatu and Berlin
Sample_Berlin <- matrix(0, 20, 2)
Sample_Vanuatu <- matrix(0, 20, 2)

for (i in 1:20) {
 for (g in 1:2) {
  Sample_Berlin[i,g] <- length(which(d_Berlin$AGE_binned==i & d_Berlin$gender == g))
  Sample_Vanuatu[i,g] <- length(which(d_Vanuatu$AGE_binned==i & d_Vanuatu$gender == g))
 }
}





#Prepare data lists as input for Stan 

d_list_Berlin <- list(N = nrow(d_Berlin),              #Sample size
                      MA = 20,                         #Number of age categories
                      age = d_Berlin$AGE_binned,       #Binned ages
                      outcome = d_Berlin$choice,       #Choices in dictator game
                      gender = d_Berlin$gender,        #Gender
                      P_empirical = t(Sample_Berlin),  #Sample demography of Berlin
                      P_Pop = t(Pop_Berlin_raw),       #Population demography of Berlin
                      P_other = t(Pop_Vanuatu_raw))    #Population demography of Vanuatu


d_list_Vanuatu <- list(N = nrow(d_Vanuatu),            #Sample size
                      MA = 20,                         #Number of age categories  
                      age = d_Vanuatu$AGE_binned,      #Binned ages
                      outcome = d_Vanuatu$choice,      #Choices in dictator game
                      gender = d_Vanuatu$gender,       #Gender
                      P_empirical = t(Sample_Vanuatu), #Sample demography of Tanna, Vanuatu
                      P_Pop = t(Pop_Vanuatu_raw),      #Population demography of Vanuatu
                      P_other = t(Pop_Berlin_raw))     #Population demography of Berlin


#Here, we code the stan model for our Gaussian Process Multilevel Regression with Poststratification

MRP_House <- "

//Function for Gaussian Process kernel

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
  int<lower = 0> P_empirical[2,MA];
  int<lower = 0> P_other[2,MA];  
  int<lower = 0> P_Pop[2,MA]; 

  }

parameters {
  vector[2] alpha;            //Gender-specific intercepts
  matrix[2,MA] age_effect;    //Matrix for Gaussian process age effects
  
  //Here we define the Control parameters for the Gaussian processes; they determine how covariance changes with increasing distance in age
  real<lower=0> eta[2];
  real<lower=0> sigma[2];
  real<lower=0, upper=1> rho[2];
}

model {
  vector[N] p;

  alpha ~ normal(0, 3);
  eta ~ exponential(2);
  sigma ~ exponential(1);
  rho ~ beta(10, 1);

  //We compute age-specific offsets for each sex
  for ( i in 1:2){
   age_effect[i,] ~ multi_normal_cholesky( rep_vector(0, MA) , GPL(MA, rho[i], eta[i], sigma[i]) );
  }  

  //This is the linear model: Choice probabilities are composed of (gender-specific) intercept and gender-specific offset for each age category
  for ( i in 1:N ) {
   p[i] = alpha[gender[i]] + age_effect[gender[i],age[i]];
  }

 //Binomial likelihood for observed choices
 outcome ~ binomial_logit(1, p);

}
 
 //We use the generated quantities section for the poststratification
 
generated quantities{
   real<lower = 0, upper = 1> p_source;  // This is value for p in the sample
   real<lower = 0, upper = 1> p_pop;     // This is value for p in the population from which sample was taken
   real<lower = 0, upper = 1> p_other;   // This is value for p in the other population 

   real expect_pos = 0;
   int total = 0;
   
   //Here we compute predictions for each age and gender class
   vector[MA] pred_p_m;
   vector[MA] pred_p_f;
   pred_p_m = inv_logit(alpha[1] + age_effect[1,]');
   pred_p_f = inv_logit(alpha[2] + age_effect[2,]');


   //Here we do the actual poststratification
    
   //Sample  
   for (a in 1:2)
        for (b in 1:MA){
         total += P_empirical[a,b];
         expect_pos += P_empirical[a,b] * inv_logit(alpha[a] + age_effect[a,b]);
     }
    p_source = expect_pos / total;
   
   
   
    //Poststratified to sample population  
    total = 0;
    expect_pos = 0;
    for (a in 1:2)
      for (b in 1:MA){
        total += P_Pop[a,b];
        expect_pos += P_Pop[a,b] * inv_logit(alpha[a] + age_effect[a,b]);
      }
    p_pop = expect_pos / total;

    
    
    //Poststratified to other population  
    total = 0;
    expect_pos = 0;
       for (a in 1:2)
         for (b in 1:MA){
          total += P_other[a,b];
          expect_pos += P_other[a,b] * inv_logit(alpha[a] + age_effect[a,b]);
      }
     p_other = expect_pos / total;
   
   

}

"


#Pass code to Rstan and run models
m_Berlin  <- stan( model_code  = MRP_House , data=d_list_Berlin  ,iter = 5000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.99, max_treedepth = 13))  
m_Vanuatu <- stan( model_code  = MRP_House , data=d_list_Vanuatu ,iter = 5000, cores = 4, seed=1, chains=4, control = list(adapt_delta=0.99, max_treedepth = 13))  



######
####
###
##
# Plotting Code
##
###
####
####


###
##
#Gender and age-specific estimates (Fig.S4 in the appendix)
##
###

#graphics.off()
#png("AgeCurvesReal.png", res = 900, height = 10, width = 22, units = "cm")

col.pal <- brewer.pal(8, "Dark2") 
seqoverall <- seq

par(mfrow = c(1,2),
    mar = c(1,2,1.2,0.5), 
    oma = c(2.5,2,1,0))

###
##
# Berlin
##
###
samp <- extract.samples(m_Berlin)
samp_pred_p <- apply(samp$pred_p_m,2,quantile,probs=c(0.055,0.5,0.945))
Lower <- samp_pred_p[1,]
Median <- samp_pred_p[2,]
Higher <- samp_pred_p[3,]
Observed_m <- ifelse(c(1:20) %in% unique(d_list_Berlin$age[d_list_Berlin$gender==1]),"Observed Age","Unobserved Age")
Observed_f <- ifelse(c(1:20) %in% unique(d_list_Berlin$age[d_list_Berlin$gender==2]),"Observed Age","Unobserved Age")
plot(NA, type = "l", xlim = c(0,20), ylim = c(0,1), xlab = "", ylab = "", lwd = 1, lty = 2, main= "Berlin (Germany)")
points(x= 1:20 - 0.25, Median, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.2)), pch = 16)
segments(x0 = 1:20-0.25, y0 = Lower, x1=  1:20-0.25, y1 = Higher, lwd = 1, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.2)))
samp_pred_p <- apply(samp$pred_p_f,2,quantile,probs=c(0.055,0.5,0.945))
Lower <- samp_pred_p[1,]
Median <- samp_pred_p[2,]
Higher <- samp_pred_p[3,]
par(new = TRUE)
plot(NA, type = "l", xlim = c(0,20), ylim = c(0,1), xlab = "", ylab = "", lwd = 1, lty = 2, main= "")
points(x= 1:20 + 0.25,Median, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.2)), pch = 17)
segments(x0 = 1:20+0.25, y0 = Lower, x1=  1:20+0.25, y1 = Higher, lwd = 1, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.2)))


###
##
# Vanuatu
##
###
samp <- extract.samples(m_Vanuatu)
samp_pred_p <- apply(samp$pred_p_m,2,quantile,probs=c(0.055,0.5,0.945))
Lower <- samp_pred_p[1,]
Median <- samp_pred_p[2,]
Higher <- samp_pred_p[3,]
Observed_m <- ifelse(c(1:20) %in% unique(d_list_Vanuatu$age[d_list_Vanuatu$gender==1]),"Observed Age","Unobserved Age")
Observed_f <- ifelse(c(1:20) %in% unique(d_list_Vanuatu$age[d_list_Vanuatu$gender==2]),"Observed Age","Unobserved Age")
plot(NA, type = "l", xlim = c(0,20), ylim = c(0,1), xlab = "", ylab = "", lwd = 1, lty = 2, main= "Tanna (Vanuatu)")
points(x= 1:20 - 0.25, Median, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.2)), pch = 16)
segments(x0 = 1:20-0.25, y0 = Lower, x1=  1:20-0.25, y1 = Higher, lwd = 1, col = ifelse(Observed_m == "Observed Age", alpha(col.pal[2],alpha = 1), alpha(col.pal[2],alpha = 0.2)))
legend("topleft", c("Male","Female"), col = c(alpha(col.pal[2],alpha = 1), alpha(col.pal[3],alpha = 1)), pch = c(16,17), lwd = 1, lty=1, bty="n", cex = 1.2)
samp_pred_p <- apply(samp$pred_p_f,2,quantile,probs=c(0.055,0.5,0.945))
Lower <- samp_pred_p[1,]
Median <- samp_pred_p[2,]
Higher <- samp_pred_p[3,]
par(new = TRUE)
plot(NA, type = "l", xlim = c(0,20), ylim = c(0,1), xlab = "", ylab = "", lwd = 1, lty = 2, main= "")
points(x= 1:20 + 0.25,Median, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.2)), pch = 17)
segments(x0 = 1:20+0.25, y0 = Lower, x1=  1:20+0.25, y1 = Higher, lwd = 1, col = ifelse(Observed_f == "Observed Age", alpha(col.pal[3],alpha = 1), alpha(col.pal[3],alpha = 0.2)))

mtext("Age category", side = 1, outer = TRUE, line = 1.5, cex = 1.2)
mtext("Probability of choosing prosocial option", side = 2, line = 1,outer = TRUE)

#dev.off()



###
##
# Demographic standardization comparing Tanna, Vanuatu, and Berlin, Germany (Fig.4)
##
###

#graphics.off()
#png("DemostandBERVAT.png", res = 900, height = 19, width = 15, units = "cm")

par(mfrow = c(3,2), 
    mar = c(4,1,2,2), 
    oma = c(0.1,4,0,0.1))

# Create pyramid plot labels
labels1 <- matrix("", 20, 2)
labels1[,1] <- c("1-5", "","11-15","","21-25","","31-35","","41-45","",
                 "51-55","","61-65","","71-75","","81-85","","91-95","")
labels2 <- matrix("", 20, 2)

#Population Pyramids
par(mar=pyramid.plot(Pop_Vanuatu[,1],Pop_Vanuatu[,2],top.labels=c("", "Tanna (Vanuatu)",""), ppmar=c(2,1,3,1),raxlab = seq(0,10,2),laxlab = seq(0,10,2), xlim = c(10,10),labelcex=1, unit = "",show.values=F, labels = labels1, lxcol = col.pal[2], rxcol = col.pal[3],space = 0.2,gap = 0))
mtext("Population", side = 2, outer = F, line = 3.5, cex = 1.3)
mtext("Age class", side = 2, outer = F, line = 0.2, cex = 1)
par(mar=pyramid.plot(Pop_Berlin[,1],Pop_Berlin[,2],top.labels=c("", "Berlin (Germany)",""), ppmar=c(2,1,3,1),raxlab = seq(0,10,2),laxlab = seq(0,10,2), xlim = c(10,10),labelcex=1, unit = "",show.values=F, labels = labels2, lxcol = col.pal[2], rxcol = col.pal[3],space = 0.2,gap = 0))
mtext("Share of population per age class and gender [%]", side = 1,line = 4.5,at = -11, outer = F, cex = 0.9)
legend("topleft", c("Male", "Female"), col = c(col.pal[2], col.pal[3]),cex = 1.1, lty = 1,lwd = 8, bty = "n" )

#Sample Pyramids
par(mar=pyramid.plot(Sample_Vanuatu[,1],Sample_Vanuatu[,2],top.labels=c("", "",""), ppmar=c(2,1,3,1),raxlab = seq(0,30,10),laxlab = seq(0,30,10), xlim = c(30,30),labelcex=1, unit = "",show.values=F, labels = labels1, lxcol = col.pal[2], rxcol = col.pal[3],space = 0.2,gap = 0))
mtext("Sample", side = 2, outer = F, line = 3.5, cex = 1.3)
mtext("Age class", side = 2, outer = F, line = 0.2, cex = 1)
par(mar=pyramid.plot(Sample_Berlin[,1],Sample_Berlin[,2],top.labels=c("", "",""), ppmar=c(2,1,3,1),raxlab = seq(0,30,10),laxlab = seq(0,30,10), xlim = c(30,30),labelcex=1, unit = "",show.values=F, labels = labels2, lxcol = col.pal[2], rxcol = col.pal[3],space = 0.2,gap = 0))
mtext("Number of individuals per age class and gender", side = 1,line = 4.5,at = -33, outer = F, cex = 0.9)


#Densities Tanna
s_Berlin <- extract.samples(m_Berlin)
s_Vanuatu <- extract.samples(m_Vanuatu)
dens <- density(s_Vanuatu$p_source)
x1 <- min(which(dens$x >= quantile(s_Vanuatu$p_source, 0)))  
x2 <- max(which(dens$x <  quantile(s_Vanuatu$p_source, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,15), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[1],alpha = 0.5), border = NA))
par(new=TRUE)
dens <- density(s_Vanuatu$p_pop)
x1 <- min(which(dens$x >= quantile(s_Vanuatu$p_pop, 0)))  
x2 <- max(which(dens$x <  quantile(s_Vanuatu$p_pop, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,15), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[6],alpha = 0.5), border = NA))
mtext("Density", side = 2,line = 3, outer = F, cex = 1)
par(new=TRUE)
dens <- density(s_Vanuatu$p_other)
x1 <- min(which(dens$x >= quantile(s_Vanuatu$p_other, 0)))  
x2 <- max(which(dens$x <  quantile(s_Vanuatu$p_other, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,15), type="n", ann = FALSE, bty = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[8],alpha = 0.5), border = NA))
legend("topleft", c("Empirical estimate","Poststratified to sample population", "Poststratified to other population"), col = c(alpha(col.pal[1],alpha = 0.5),alpha(col.pal[6],alpha = 0.5),alpha(col.pal[8],alpha = 0.5)), lwd = 6, bty="n", cex = 1)


#Densities Berlin
dens <- density(s_Berlin$p_source)
x1 <- min(which(dens$x >= quantile(s_Berlin$p_source, 0)))  
x2 <- max(which(dens$x <  quantile(s_Berlin$p_source, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,15), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[1],alpha = 0.5), border = NA))
par(new=TRUE)
dens <- density(s_Berlin$p_pop)
x1 <- min(which(dens$x >= quantile(s_Berlin$p_pop, 0)))  
x2 <- max(which(dens$x <  quantile(s_Berlin$p_pop, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,15), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[6],alpha = 0.5), border = NA))
par(new=TRUE)
dens <- density(s_Berlin$p_other)
x1 <- min(which(dens$x >= quantile(s_Berlin$p_other, 0)))  
x2 <- max(which(dens$x <  quantile(s_Berlin$p_other, 1)))
plot(dens, xlim = c(0,1), ylim = c(0,15), type="n", ann = FALSE, bty = "n", yaxt = "n")
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=alpha(col.pal[8],alpha = 0.5), border = NA))
mtext("Probability of choosing prosocial option", side = 1,line = 2.9,at=-0.1, outer = F, cex = 0.9)

#dev.off()







