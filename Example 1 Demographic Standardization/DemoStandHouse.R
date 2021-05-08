
library(readr)
library(plotrix)

#Population Distribution of Vanuatu and Berlin
Pop_Berlin <- read.csv("C:/Users/Dominik/Desktop/Berlin-2020.csv")
Pop_Vanuatu <- read.csv("C:/Users/Dominik/Desktop/Vanuatu-2017.csv")
Pop_Vanuatu <- Pop_Vanuatu[-nrow(Pop_Vanuatu),]

# Creat pyramid plot labels
labels1 <- matrix("", 20, 2)
labels1[,1] <- Pop_Vanuatu[,1]
labels2 <- matrix("", 20, 2)

Pop_Berlin <- as.matrix(Pop_Berlin[,c(2,3)])
Pop_Vanuatu <- as.matrix(Pop_Vanuatu[,c(2,3)])

Pop_Berlin_raw <- Pop_Berlin
Pop_Berlin <- (Pop_Berlin/sum(Pop_Berlin))*100
Pop_Vanuatu_raw <- Pop_Vanuatu

Pop_Vanuatu <- (Pop_Vanuatu/sum(Pop_Vanuatu))*100

#Sample Distribution of Vanuatu and Berlin
#Create demography age pyramids
setwd("~/GitHub/Cross_Cultural_Generalizability")

data_adult <- read.csv("House_data/Model_1a_1b_1c_data.csv")
data_adult <- data_adult[,c("SUBJECT_ID","GENDER_1female","fieldid", "AGE_in_years","T1_ad_choice_1yes")]
data_adult$choice <- data_adult$T1_ad_choice_1yes
data_adult$T1_ad_choice_1yes <- NULL


data_children <- read.csv("House_data/Model_4a_4b_4c_4d_data.csv")
data_children <- data_children[,c("SUBJECT_ID","GENDER_1female","fieldid", "AGE_in_years","T1_choice_1yes")]
data_children$choice <- data_children$T1_choice_1yes
data_children$T1_choice_1yes <- NULL


Society <- c("Berlin (GER)","La Plata (ARG)","Phoenix (USA)", "Pune (IND)", "Shuar (ECU)", "Wiichi (ARG)", "Tanna (WUT)", "Hadza (TZA)")


data_comb <- rbind(data_children, data_adult)
data_comb <- data_comb[with(data_comb, order(data_comb$fieldid, data_comb$AGE_in_years)),]
data_comb$gender <- data_comb$GENDER_1female +1
MA <- max(round(data_comb$AGE_in_years))

d_Berlin <- data_comb[which(data_comb$fieldid == 1),]
d_Vanuatu <- data_comb[which(data_comb$fieldid == 7),]

age_upper <- seq(5, 100, 5)

d_Berlin$AGE_binned <- sapply(1:nrow(d_Berlin), function(i) which.max( 1/(age_upper - d_Berlin$AGE_in_years[i]) ) )
d_Vanuatu$AGE_binned <- sapply(1:nrow(d_Vanuatu), function(i) which.max( 1/(age_upper - d_Vanuatu$AGE_in_years[i]) ) )


Sample_Berlin <- matrix(0, 20, 2)
Sample_Vanuatu <- matrix(0, 20, 2)

for (i in 1:20) {
 for (g in 1:2) {
  Sample_Berlin[i,g] <- length(which(d_Berlin$AGE_binned==i & d_Berlin$gender == g))
  Sample_Vanuatu[i,g] <- length(which(d_Vanuatu$AGE_binned==i & d_Vanuatu$gender == g))
  
 }
}










#Prepare for stan

d_list_Berlin <- list(N = nrow(d_Berlin), 
                      MA = 20,
                      age = d_Berlin$AGE_binned,
                      outcome = d_Berlin$choice, 
                      gender = d_Berlin$gender, 
                      P_empirical = t(Sample_Berlin),
                      P_Pop = t(Pop_Berlin_raw),
                      P_other = t(Sample_Vanuatu))


d_list_Vanuatu <- list(N = nrow(d_Vanuatu), 
                      MA = 20,
                      age = d_Vanuatu$AGE_binned,
                      outcome = d_Vanuatu$choice, 
                      gender = d_Vanuatu$gender, 
                      P_empirical = t(Sample_Vanuatu),
                      P_Pop = t(Pop_Vanuatu_raw),
                      P_other = t(Sample_Berlin))


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
  int<lower = 0> P_empirical[2,MA];   // Here we enter data matrix with demographic constitution of target population
  int<lower = 0> P_other[2,MA];   // Here we enter data matrix with demographic constitution of target population
  int<lower = 0> P_Pop[2,MA];   // Here we enter data matrix with demographic constitution of target population

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
   real<lower = 0, upper = 1> p_source;  // This is value for p in the source population 
   real<lower = 0, upper = 1> p_pop;  // This is value for p in the target population 
   real<lower = 0, upper = 1> p_other;  // This is value for p in the target population 

   int total = 0;
   vector[MA] pred_p_m;
   vector[MA] pred_p_f;

   for (a in 1:2)
        for (b in 1:MA){
         total += P_empirical[a,b];
         expect_pos += P_empirical[a,b] * inv_logit(alpha[a] + age_effect[a,b]);
     }
     
   p_source = expect_pos / total;
   
    total = 0;
    expect_pos = 0;
    
       for (a in 1:2)
         for (b in 1:MA){
          total += P_other[a,b];
          expect_pos += P_other[a,b] * inv_logit(alpha[a] + age_effect[a,b]);
      }
      
    p_other = expect_pos / total;
   
   
       total = 0;
    expect_pos = 0;
    
       for (a in 1:2)
         for (b in 1:MA){
          total += P_Pop[a,b];
          expect_pos += P_Pop[a,b] * inv_logit(alpha[a] + age_effect[a,b]);
      }
      
    p_pop = expect_pos / total;
   
   
  pred_p_m = inv_logit(alpha[1] + age_effect[1,]');
  pred_p_f = inv_logit(alpha[2] + age_effect[2,]');
}

"

}




library(rstan)
m_Berlin <- stan( model_code  = m2a_MRP_GP_gender_same , data=d_list_Berlin ,iter = 5000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.95, max_treedepth = 13))  
m_Vanuatu <- stan( model_code  = m2a_MRP_GP_gender_same , data=d_list_Vanuatu ,iter = 5000, cores = 1, seed=1, chains=1, control = list(adapt_delta=0.95, max_treedepth = 13))  


















library(scales)
#color stuff
require(RColorBrewer)#load package
col.pal <- brewer.pal(8, "Dark2") #create a pallette which you loop over for corresponding values


graphics.off()
png("DemostandLarge2.png", res = 900, height = 25, width = 25, units = "cm")


par(mfrow = c(3,2), 
    mar = c(3,1,3,2), 
    oma = c(1,4,2,0.1))
# Population 1
par(mar=pyramid.plot(Pop_Vanuatu[,1],Pop_Vanuatu[,2],top.labels=c("", "Vanuatu",""), ppmar=c(2,1,3,1), xlim = c(10,10),labelcex=1, unit = "",show.values=F, labels = labels1, lxcol = col.pal[2], rxcol = col.pal[3],space = 0.2,gap = 0))
mtext("Population", side = 2, outer = F, line = 3.5, cex = 1.3)
mtext("Age class", side = 2, outer = F, line = 0.2, cex = 1)

legend("topright", c("Male", "Female"), col = c(col.pal[2], col.pal[3]),cex = 1, lty = 1,lwd = 5, bty = "n" )
par(mar=pyramid.plot(Pop_Berlin[,1],Pop_Berlin[,2],top.labels=c("", "Berlin",""), ppmar=c(2,1,3,1), xlim = c(10,10),labelcex=1, unit = "",show.values=F, labels = labels2, lxcol = col.pal[2], rxcol = col.pal[3],space = 0.2,gap = 0))
mtext("Share of population per age class and gender [%]", side = 1,line = 4,at = -10, outer = F, cex = 0.9)



par(mar=pyramid.plot(Sample_Vanuatu[,1],Sample_Vanuatu[,2],top.labels=c("", "",""), ppmar=c(2,1,3,1), xlim = c(30,30),labelcex=1, unit = "",show.values=F, labels = labels1, lxcol = col.pal[2], rxcol = col.pal[3],space = 0.2,gap = 0))
par(mar=pyramid.plot(Sample_Berlin[,1],Sample_Berlin[,2],top.labels=c("", "",""), ppmar=c(2,1,3,1), xlim = c(30,30),labelcex=1, unit = "",show.values=F, labels = labels2, lxcol = col.pal[2], rxcol = col.pal[3],space = 0.2,gap = 0))





dev.off()