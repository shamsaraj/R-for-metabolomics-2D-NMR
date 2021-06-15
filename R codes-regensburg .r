#source("http://bioconductor.org/biocLite.R")
#biocLite()
library (vsn)
library(multtest)
library(JADE)
library (RColorBrewer)
library(ggplot2)
#setwd("C:/Downloads")
#setwd ("/home/shj06387/Desktop")

x<-read.table(file.choose(), sep=",", header=T, row.names=1)
#x<-read.table("/home/shj06387/Desktop/export_txt-1d", sep=",", header=T, row.names=1)
#x<-read.table("E:/Scratch/export_txt-1d.csv", sep=",", header=T, row.names=1)
x1<-as.matrix (t(x))
kick.out<-apply(x1,1,function(z){all(z==0)})#find zeros (zero in all samples)
x2<-x1[!kick.out,] #kick zeros out

fit<-vsn2(x2)#variance stabilization 
x3=predict(fit,newdata=x2)#create an object of the same type as x2 from vsn object
qqnorm(x3) #control
meanSdPlot(x3) #control

l1<-numeric(10)
for (i in 6:10){l1[i]=1} # setting of labels
t<-mt.teststat(x3[,1:10],classlabel=l1,test="t") # only first 2 groups
pt<-2*pt(-abs(t),df=8)
pdf('hist-1d.pdf')
hist(pt,50)
dev.off()

ica1<-JADE(x2,n.comp = 2,maxiter = 12200)
#Par (mfrow = c (1,2)) 
#pdf('ICA-1d.pdf')
#plot(ica1$A[1:5,1],ica1$A[1:5,2],xlab="IC1",ylab="IC2", pch=19, cex= 1.4, col="green")
#plot the control samples in green circle
#points(ica1$A[6:10,1],ica1$A[6:10,2],pch=17,col="red")#plot the patient samples in red triangle
plot(ica1$A[,1],ica1$A[,2], xlab="IC1",ylab="IC2", pch=19, cex=1.4, col= c("green", "green", "green", "green", "green", "red", "red", "red", "red", "red"))
#coloring for two groups each have four members
#dev.off() 

#PCA
pcdat1 = prcomp(t(x2))
#screeplot(pcdat1)
#biplot(pcdat1)
pdf('PCA-1d.pdf')
#plot(pcdat1$x[1:5,1],pcdat1$x[1:5,2],xlab="PC1",ylab="PC2", pch=19, cex= 1.4,col="green")
#plotthe control samples in green circle
#points(pcdat1$x[6:10,1],pcdat1$x[6:10,2],pch=17,col="red")
#plot the patient samples in red triangle
plot(pcdat1$x[,1], pcdat1$x[,2], xlab="PC1",ylab="PC2", pch=19, cex=1.4, col= c("green", "green", "green", "green", "green", "red", "red", "red", "red", "red"))
dev.off() 
getwd()

#Loadings
ica1<-JADE(t(x2),n.comp = 2,maxiter = 12200)
#save the frequencies in a numeric format
frequenzenx2<-signif(as.numeric(sub('X........', '',dimnames(x2)[[1]][])),4)
#plot the loadings of the first dimension
pdf ('Loadings-IC1.pdf')
plot(frequenzenx2,ica1$A[,1],type="s",xlab="Frequency /ppm",ylab="Loading IC1",xlim=c(5.3,0.75))
dev.off()

#pcdat2 = princomp(x2)#This command needs more sample than variables. Advantage is it can return loadings
#loadings(pcdat2) 
#plot(pcdat2$scores[,1], pcdat2$scores[,2])

#The code to sort the loadings according to their importance
ica1<-JADE(t(x2),n.comp = 2,maxiter = 12200)#Showes the components based on the loadings
Loadings1<-sort.int(abs(ica1$A[,1]),decreasing=T,index.return=T)#Sorts the Loadings and gives the index to  them
sortedLoadings1<-matrix(data=NA,length(x2[,1]),3,dimnames=list(1:length(x2[,1]),c("Loading Factor", "Frequency","P-value")))
sortedLoadings1[,2]<-signif(as.numeric(sub('X........', '',dimnames(x2)[[1]][Loadings1$ix])),4)#adds frequencies to the second column
sortedLoadings1[,1]<-signif(ica1$A[Loadings1$ix,1],digits=5)#adds sorted weight factors to the second column
#This loop is for calculating p values. It adds them to the third column
for (i in 1:length(x2[,1])){
 j<-Loadings1$ix[i]
 sortedLoadings1[i,3]<-t.test(x2[j,1:5],x2[j,6:10])$p.value
}
#These codes are for IC2 
Loadings2<-sort.int(abs(ica1$A[,2]),decreasing=T,index.return=T)#gewichtsfaktoren sortieren
sortedLoadings2<-matrix(data=NA,length(x2[,1]),3,dimnames=list(1:length(x2[,1]),c("Loading Factor", "Frequency", "P-value")))
sortedLoadings2[,2]<-signif(as.numeric(sub('X........', '',dimnames(x2)[[1]][Loadings2$ix])),4)
sortedLoadings2[,1]<-signif(ica1$A[Loadings2$ix,2],digits=5)#gewichtsfaktoren
for (i in 1:length(x2[,1])){
 j<-Loadings2$ix[i]
 sortedLoadings2[i,3]<-t.test(x2[j,1:5],x2[j,6:10])$p.value
}

write.table(sortedLoadings1,file="IC1Loadings.csv",sep=",",quote=F)
write.table(sortedLoadings2,file="IC2Loadings.csv",sep=",",quote=F)

#Heatmap for highest 30 ranked bins of (both) IC1
bucketnumbers=30
Heatmapmatrix<-mat.or.vec(10, bucketnumbers)
Heatmapmatrix[,1:bucketnumbers]<-t(x2) [1:10,Loadings1$ix[1:bucketnumbers]]
#Heatmapmatrix[,bucketnumbers+1:bucketnumbers*2]<-t(x2) [1:10,Loadings2$ix[1:50]]
Heatmapmatrix<- t(Heatmapmatrix)
##newrownames<-read.csv(file.choose())
##attach(newrownames)
##row.names (Heatmapmatrix) <- frequencyinppm
row.names (Heatmapmatrix) [1:bucketnumbers]<- sortedLoadings1[1:bucketnumbers,2] 
#rownames (Heatmapmatrix) [bucketnumbers+1:bucketnumbers*2]<- sortedLoadings2[1:bucketnumbers,2] 
Columnvector<-c("Control1","Control2", "Control3", "Control4", "Control5", "NASH1", "NASH2","NASH3","NASH4", "NASH5")
colnames (Heatmapmatrix)<- Columnvector
pdf ("heatmap.pdf")
mypalette<-brewer.pal(11,"RdBu")
#mypalette<-brewer.pal(9,"Blues")
NASH_heatmap <- heatmap(Heatmapmatrix, Rowv=NULL, Colv=NULL,
#col = mypalette, scale="row", margins=c(5,10))
#col = heat.colors(256), scale="row", margins=c(5,10))
#col = cm.colors(256), scale="row", margins=c(5,10))
#col = rev(heat.colors(16)), scale="row", margins=c(5,10))
col = rev(mypalette), scale="row", margins=c(5,10), cexRow=1)
dev.off()

#Rowvector<-c("1.185 Isopropyl Alcohol ?", "1.175 Isopropyl Alcohol ?" ,"2.235 Acetone ?", "3.235 Betaine - Glycerophosphocholine - Glucose*", "1.285 ?", "1.145 ?", "1.155 ?", "0.865 ?", "3.425 Taurine - Glucose?", "1.305 ?", "1.295 ?", "3.895 Glucose -  Betaine - Glycerophosphocholine?*", "0.855 ?", "1.325 Threonine - Lactate ?*", "1.315 Threonine - Lactate?*", "3.465 Glucose ?", "3.485 Glucose ?", "3.725 Glucose", "1.335 Threonine - Lactate?*", "3.495 Glucose ?", "3.255 Glucose - Betaine - Glycerophosphocholine - Taurine*", "3.415 Glucose - Taurine ?", "3.545 Glucose ?", "3.845 Glucose - Glucose - Serine*", "3.225 Betaine - Glycerophosphocholine - Glucose - Arginine*", "0.875 ?", "1.275 ?", "3.455 Glucose ?", "3.735 Glucose", "3.915 Glucose - Glycerophosphocholine ?*")

#Heatmap by ggplot2
Heatdataframe=as.data.frame(Heatmapmatrix)
ggstructure(Heatdataframe) + scale_fill_continuous(low = "white", high = "steelblue") 
#ggstructure(Heatdataframe) + scale_fill_continuous(low = "green", high = "red") 

#Calculation of adjusted p values according to B and H method
p1dbinesvector<-sortedLoadings1[,3]#Original p values vector
adjustedpvlues<-mat.or.vec(751, 2)
adjustedpvlues [,1]<-sortedLoadings1[,2]
adjustedpvlues [,2]<-p.adjust(p1dbinesvector, method = "BH", n = length (p1dbinesvector))
write.table(adjustedpvlues,file="adjustedpvlues.csv",sep=",",quote=F)



