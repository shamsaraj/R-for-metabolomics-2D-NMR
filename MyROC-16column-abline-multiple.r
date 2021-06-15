# ROC curve(s) with ROCR
# Chupakhin Vladimir (chupvl@gmail.com)
# loading ROCR library
#loading enrichvs library
library("ROCR")# package for ROC curve 
library (enrichvs)# package for calculation of enrichment factor
setwd("D:/test4/ada/roc")##############working directory
header <- c ("name", "v" , "c" , "d", "a" , "vcda", "vcd" , "vda", "vca", "cda" , "vc", "vd" , "va", "cd" , "ca", "da")#manualy defined headers
X<- c(2,12)
for (i in (X)){#a loop for reading all columns 
title <- header[i]
activesfile<-"actives.csv"###### inputfile name for actives
decoysfile<-"decoys.csv"######	input file name for decoys
table1 <- read.csv(activesfile,  sep=",", header=T)
score1 <- table1 [,i]
GGG1 <- score1[!is.na(score1)]#omiting NAs
table2 <- read.csv(decoysfile,  sep=",", header=T)
score2 <- table2 [,i]
GGG2 <- score2[!is.na(score2)]
y <- length (GGG1)
z <- length (GGG2)
x <- c("actives", "decoys")
xx <- c(1,0)
GGG <- c (GGG1, GGG2)
labels <- (rep(x ,c(y,z)))#build a vector for ROCR labels
labels2 <- (rep(xx ,c(y,z)))#build a vector for enrichvs labels
if (i==2){
pred1 <- prediction (GGG, labels)
perf1 <- performance (pred1, "tpr", "fpr")
}
if (i==12){
pred2 <- prediction (GGG, labels)
perf2 <- performance (pred2, "tpr", "fpr")
}
}
# changing params for the ROC plot - width, etc
par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=1.3,cex.lab=1.4)
#plotting the ROC curve
pdftitle<-paste(activesfile, title, "pdf", sep=".")
pdf (pdftitle)
plot(perf1,col="blue",lty=1, lwd=2)
plot(perf2,add = TRUE, col="red",lty=3, lwd=4)
abline(a=0, b= 1, lty=3, lwd=2)
# calculating AUC
auc1 <- performance(pred1,"auc")
# now converting S4 class to vector
auc1 <- unlist(slot(auc1, "y.values"))
# adding min and max ROC AUC to the center of the plot
AUC1 <- round(auc, digits = 3)
AUC1 <- paste(c("AUC Vina = "),AUC1,sep="")
legend(0.5,0.3,c(AUC1),border="white",cex=1,box.col = "white")
# calculating AUC
auc2 <- performance(pred2,"auc")
# now converting S4 class to vector
auc2 <- unlist(slot(auc2, "y.values"))
# adding min and max ROC AUC to the center of the plot
AUC2 <- round(auc2, digits = 3)
AUC2 <- paste(c("AUC Vina + DSX = "),AUC2,sep="")
legend(0.5,0.2,c(AUC2),border="white",cex=1,box.col = "white")
dev.off()
dev.off()





