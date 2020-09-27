`+` <- function(e1, e2) {
  if (is.character(e1) | is.character(e2)) {
    paste0(e1, e2)
  } else {
    base::`+`(e1, e2)
  }
}

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(survivalROC)



beta_file='/Users/dilrajkaur/Desktop/UCEC/UCEC_HRmedian_raw.csv'
#tcga_datafile="/Users/dilrajkaur/Desktop/UCEC/UCEC.RDS"
tcga_datafile2="/Users/dilrajkaur/Desktop/UCEC/UCEC-final.csv"

beta<-read.csv(file=beta_file,header =TRUE, sep = ",", dec = ".")
data1<-read.csv(file=tcga_datafile2,header =TRUE, sep = ",", dec = ".")
#data1<-readRDS(tcga_datafile)
genes <- scan(file = '/Users/dilrajkaur/Desktop//UCEC/UCEC_genes10.csv', what = 'character', sep = ',')

PI=0
v=0
for(i in seq(from=1, to=length(genes), by=1))#no of cancers
{
  if (beta[beta$Gene==genes[i],][4][,1]<0.05)
  {
    v=v+1
    k=beta[beta$Gene==genes[i],][2][,1]
    if (k>0)
    {
      b=1*(select(data1, genes[i])[,1]>median(select(data1, genes[i])[,1]))}
    
    if (k<0)
    {b=1*(select(data1, genes[i])[,1]<median(select(data1, genes[i])[,1]))}
    
    PI=PI+b
  }
}
data1$PI=PI

cutoff= 1#years
marker1= 1*(data1$PI>7)

Mayo4.2= survivalROC(Stime=data1$OS.time/365,  
                     status=data1$vital_status,      
                     marker = marker1,     
                     predict.time =  cutoff, method="KM")
plot(Mayo4.2$FP, Mayo4.2$TP, type="l", xlim=c(0,1), ylim=c(0,1.1),   
     xlab=paste( "FP", "", "AUROC = ",round(Mayo4.2$AUC,3)), 
     ylab="TP",main="Mayoscore 4, Method = KM 
 Year = 1", col = "red")
abline(0,1)
Mayo4.2