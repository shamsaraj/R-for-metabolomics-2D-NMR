# ROC curve(s) with ROCR
# Chupakhin Vladimir (chupvl@gmail.com)
# loading ROCR library
#loading enrichvs library
library("ROCR")# package for ROC curve 
library (enrichvs)# package for calculation of enrichment factor
setwd("E:/")##############working directory
header <- c ("name", "v" , "c" , "d", "a" , "vcda", "vcd" , "vda", "vca", "cda" , "vc", "vd" , "va", "cd" , "ca", "da")#manualy defined headers
X<- c(1,2,3,4)
for (i in (X)){#a loop for reading all files
title <- header[i]
if (i==1){
activesfile<-"E:/SP/3F17/actives.csv"###### inputfile name for actives
decoysfile<-"E:/SP/3F17/decoys.csv"######	input file name for decoys
}
if (i==2){
activesfile<-"E:/SP/3F19/actives.csv"###### inputfile name for actives
decoysfile<-"E:/SP/3F19/decoys.csv"######	input file name for decoys
}
if (i==3){
activesfile<-"E:/SP/3N2V/actives.csv"###### inputfile name for actives
decoysfile<-"E:/SP/3N2V/decoys.csv"######	input file name for decoys
}
if (i==4){
activesfile<-"E:/SP/3NX7/actives.csv"###### inputfile name for actives
decoysfile<-"E:/SP/3NX7/decoys.csv"######	input file name for decoys
}
table1 <- read.csv(activesfile,  sep=";", header=T)
score1 <- table1 [,2]
GGG1 <- score1[!is.na(score1)]#omiting NAs
table2 <- read.csv(decoysfile,  sep=";", header=T)
score2 <- table2 [,2]
GGG2 <- score2[!is.na(score2)]
y <- length (GGG1)
z <- length (GGG2)
x <- c("actives", "decoys")
xx <- c(1,0)
GGG <- c (GGG1, GGG2)
labels <- (rep(x ,c(y,z)))#build a vector for ROCR labels
labels2 <- (rep(xx ,c(y,z)))#build a vector for enrichvs labels
if (i==1){
pred1 <- prediction (GGG, labels)
perf1 <- performance (pred1, "tpr", "fpr")
}
if (i==2){
pred2 <- prediction (GGG, labels)
perf2 <- performance (pred2, "tpr", "fpr")
}
if (i==3){
pred3 <- prediction (GGG, labels)
perf3 <- performance (pred3, "tpr", "fpr")
}
if (i==4){
pred4 <- prediction (GGG, labels)
perf4 <- performance (pred4, "tpr", "fpr")
}
}
# changing params for the ROC plot - width, etc
par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=1.3,cex.lab=1.4)
#plotting the ROC curve
pdftitle<-"cross.pdf"#paste(activesfile, title, "pdf", sep=".")
pdf (pdftitle)
plot(perf1,col="black",lty=1, lwd=1)
plot(perf2,add = TRUE, col="black",lty=3, lwd=3)
plot(perf3,add = TRUE, col="black",lty=2, lwd=2)
plot(perf4,add = TRUE, col="black",lty=3, lwd=5)
abline(a=0, b= 1, lty=3, lwd=1)
segments(0.8, 0.36, 0.92, 0.36, col= 'black',lty=1, lwd=1)
segments(0.8, 0.28, 0.92, 0.28, col= 'black',lty=3, lwd=3)
segments(0.8, 0.20, 0.92, 0.20, col= 'black',lty=2, lwd=2)
segments(0.8, 0.12, 0.92, 0.12, col= 'black',lty=3, lwd=5)

# calculating AUC 1
auc1 <- performance(pred1,"auc")
# now converting S4 class to vector
auc1 <- unlist(slot(auc1, "y.values"))
# adding min and max ROC AUC to the center of the plot
AUC1 <- round(auc1, digits = 3)
AUC1 <- paste(c("AUC for 3F17 = "),AUC1,sep="")
legend(0.4,0.4,c(AUC1),border="white",cex=1,box.col = "white")
# calculating AUC 2
auc2 <- performance(pred2,"auc")
# now converting S4 class to vector
auc2 <- unlist(slot(auc2, "y.values"))
# adding min and max ROC AUC to the center of the plot
AUC2 <- round(auc2, digits = 3)
AUC2 <- paste(c("AUC for 3F19 = "),AUC2,sep="")
legend(0.4,0.32,c(AUC2),border="white",cex=1,box.col = "white")
# calculating AUC 3
auc3 <- performance(pred3,"auc")
# now converting S4 class to vector
auc3 <- unlist(slot(auc3, "y.values"))
# adding min and max ROC AUC to the center of the plot
AUC3 <- round(auc3, digits = 3)
AUC3 <- paste(c("AUC for 3N2V = "),AUC3,sep="")
legend(0.4,0.24,c(AUC3),border="white",cex=1,box.col = "white")
# calculating AUC 4
auc4 <- performance(pred4,"auc")
# now converting S4 class to vector
auc4 <- unlist(slot(auc4, "y.values"))
# adding min and max ROC AUC to the center of the plot
AUC4 <- round(auc4, digits = 3)
AUC4 <- paste(c("AUC for 3NX7 = "),AUC4,sep="")
legend(0.4,0.16,c(AUC4),border="white",cex=1,box.col = "white")
dev.off()
dev.off()





