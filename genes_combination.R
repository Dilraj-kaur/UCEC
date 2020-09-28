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

#cancer1=scan(file = '/Users/macbook/Desktop/APOP_PAN_CANCER/top_ten/genelists/cancers.csv', what = 'character', sep = ',')

#out_file='/Users/macbook/Desktop/pancan_voting.csv'
#write.table(cbind("Cancer","HR","p-value","wald-p","logrank-p","C","%95 CI lower","%95 CI upper","min(PI)","max(PI)","cutoff"),
#            file=out_file,row.names=F,col.names=F,sep = ',')
#for(k in seq(from=1, to=length(cancer1), by=1))
#{





#cancer='THCA'
beta_file='/Users/dilrajkaur/Desktop/UCEC/newdata/UCEC_HRmedian_raw.csv'
tcga_datafile="/Users/dilrajkaur/Desktop/UCEC/newdata/UCEC-final.csv"

beta<-read.csv(file=beta_file,header =TRUE, sep = ",", dec = ".")
data1<-read.csv(file=tcga_datafile,header =TRUE, sep = ",", dec = ".")
genes <- scan(file = '/Users/dilrajkaur/Desktop/UCEC/newdata/input_genes/UCEC_genes15.network.csv', what = 'character', sep = ',')

######################################## voting based #######################################

gsize=10
cutoff=5

out_file='/Users/dilrajkaur/Desktop/UCEC/newdata/results/network.out'+gsize+'_'+cutoff+'genes.csv'
write.table(cbind("genes","genesize","HR","p-value","C","%95 CI lower","%95 CI upper","logrank-p","g1","g2","cutoff"),
            file=out_file,row.names=F,col.names=F,sep = ',')




f=combn(genes,gsize,FUN = list)

for(i in seq(from=1, to=length(f), by=1))#no of cancers
{
  g1=f[[i]]
  PI=0
  
  # if (('BNIP3L' %in% g1)&('TNFRSF12A' %in% g1)&('PSEN1' %in% g1))
  {
    
  for(m in seq(from=1, to=length(g1), by=1))
  {
    g=g1[m]
    k=beta[beta$Gene==g,][2][,1]
    if (k>0)
    {b=1*(select(data1, g)[,1]>median(select(data1, g)[,1]))}
  
    if (k<0)
    {b=1*(select(data1, g)[,1]<median(select(data1, g)[,1]))}
  
    PI=PI+b
  }
  
  #PI=PI+1*(data1$age_at_diagnosis>65)
  
  data1$PI=PI
  vals=(data1$PI>cutoff)
  
  surv_object <- Surv(time = data1$OS.time/30, event = data1$vital_status)
  fit1 <- survfit(surv_object~vals); 
  fit1.coxph <- coxph(surv_object~vals)
  first <- coef(summary(fit1.coxph))
  a=summary(fit1.coxph)
  
  if (first[5]<0.05)
  {write.table(cbind(paste(g1,collapse ='-'),gsize,first[2],first[5],a$concordance[1],a$conf.int[3],a$conf.int[4],a$logtest[3],fit1$n[1],fit1$n[2],cutoff),
               file=out_file,row.names=F,col.names=F,sep = ',',append = T);#output file
  }
  # {print(c(cat(g1),first[2],first[5],a$concordance[1],a$conf.int[3],a$conf.int[4],a$logtest[3]))}
  }
  

}



