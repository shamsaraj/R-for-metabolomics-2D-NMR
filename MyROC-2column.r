# ROC curve(s) with ROCR
# Chupakhin Vladimir (chupvl@gmail.com)
# loading ROCR library
#loading enrichvs library
library("ROCR")# package for ROC curve 
library (enrichvs)# package for calculation of enrichment factor
setwd("E:/1/roc")##############working directory
header <- c ("name", "energy")#manualy defined headers

for (i in (1:2)){#a loop for reading all columns 
title <- header[i]
activesfile<-"actives.csv"###### inputfile name for actives
decoysfile<-"decoys.csv"######	input file name for decoys
table1 <- read.csv(activesfile,  sep=";", header=T)
score1 <- table1 [,i]
GGG1 <- score1[!is.na(score1)]#omiting NAs
table2 <- read.csv(decoysfile,  sep=";", header=T)
score2 <- table2 [,i]
GGG2 <- score2[!is.na(score2)]
y <- length (GGG1)
z <- length (GGG2)
x <- c("actives", "decoys")
xx <- c(1,0)
GGG <- c (GGG1, GGG2)
labels <- (rep(x ,c(y,z)))#build a vector for ROCR labels
labels2 <- (rep(xx ,c(y,z)))#build a vector for enrichvs labels

#calculation of enrichment factors at different levels and auc using enrichvs
ef100<-enrichment_factor(GGG, labels2, top=1.0, decreasing=FALSE)
ef20<-enrichment_factor(GGG, labels2, top=0.2, decreasing=FALSE)
ef10<-enrichment_factor(GGG, labels2, top=0.1, decreasing=FALSE)
ef2<-enrichment_factor(GGG, labels2, top=0.02, decreasing=FALSE)
ef1<-enrichment_factor(GGG, labels2, top=0.01, decreasing=FALSE)
ef0.2<-enrichment_factor(GGG, labels2, top=0.002, decreasing=FALSE)
ef0.1<-enrichment_factor(GGG, labels2, top=0.001, decreasing=FALSE)
auc1<-auc(GGG, labels2, decreasing=FALSE, top=1.0)

#write headers for output file
if (i==1){
filelines<-paste("Scoring", "AUC","EF100%", "EF20%", "EF10%", "EF2%", "EF1%", "EF0.2%", "EF0.1%", sep=",")
write (activesfile, file= "output.txt", append = TRUE)
write (filelines, file= "output.txt", append = TRUE)
}

#omiting second column
#if (i!=2){


filelines<-paste(title, auc1, ef100, ef20, ef10, ef2, ef1, ef0.2, ef0.1, sep=",")
write (filelines, file= "output.txt", append = TRUE) 

#calculation of AUC and ploting ROC curve by ROCR
pred <- prediction (GGG, labels)
perf <- performance (pred, "tpr", "fpr")
#plot (perf, colonize=T)
# changing params for the ROC plot - width, etc
par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=1.3,cex.lab=1.4)
# plotting the ROC curve
#imagetitle<-paste(activesfile, title, "png", sep=".")
#png (imagetitle)
imagetitle<-paste(activesfile, title, "pdf", sep=".")
pdf (imagetitle)
plot(perf,col="black",lty=3, lwd=3)
# calculating AUC
auc <- performance(pred,"auc")
# now converting S4 class to vector
auc <- unlist(slot(auc, "y.values"))
# adding min and max ROC AUC to the center of the plot
AUC <- round(auc, digits = 2)
AUC <- paste(c("AUC = "),AUC,sep="")
legend(0.5,0.3,c(AUC),border="white",cex=1.7,box.col = "white")

dev.off()

}
#}

dev.off()
