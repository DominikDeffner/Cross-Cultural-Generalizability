#Create demography age pyramids

library(plotrix)

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

graphics.off()
png("Pyramids.png", width = 16,height = 26, units = "cm", res = 900)
par(mfrow = c(4,2), 
    mar = c(1,1,1,0), 
    oma = c(3.3,3.4,0,1))
for (j in 1:8) {
  
data_comb <- rbind(data_children, data_adult)
data_comb <- data_comb[with(data_comb, order(data_comb$fieldid, data_comb$AGE_in_years)),]
MA <- max(round(data_comb$AGE_in_years))
data_comb <- data_comb[which(data_comb$fieldid == j),]

d  <- list(outcome = data_comb$choice,
           ind_id = sapply(1:nrow(data_comb), function (i) which(unique(data_comb$SUBJECT_ID) == data_comb$SUBJECT_ID[i])),
           soc_id = data_comb$fieldid,
           age = round(data_comb$AGE_in_years),
           gender = as.integer(data_comb$GENDER_1female +1), # 2 female now
           N = nrow(data_comb),
           MA = max(round(data_comb$AGE_in_years))
)

labels <- matrix("", MA, 2)
if (j %in% c(1,3,5,7)){
labels[5,1] <- 5
labels[15,1] <- 15
labels[25,1] <- 25
labels[35,1] <- 35
labels[45,1] <- 45
labels[55,1] <- 55
labels[65,1] <- 65
}

a <- matrix(0, 2, MA)
for (i in 1:MA) {
  a[1,i] <- length(which(d$age==i & d$gender == 1))
  a[2,i] <- length(which(d$age==i & d$gender == 2))
}

par(mar=pyramid.plot(a[1,],a[2,],top.labels=c("",Society[j],""),ppmar=c(2,1,3,1), xlim = c(9,9),labelcex=1.2, unit = "",show.values=F, labels = labels, lxcol = "indianred", rxcol = "darkgreen", space = 0,gap = 0))

if (j == 2) legend("topright", c("Male", "Female"), col = c("indianred", "darkgreen"),cex = 1.4, lty = 1,lwd = 5, bty = "n" )

}

mtext("Number of individuals per age and gender", side = 1,line = 2, outer = TRUE, cex = 1.3)
mtext("Age class", side = 2, outer = TRUE, line = 2, cex = 1.3)
dev.off()
