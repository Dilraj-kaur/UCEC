library(survival)
library(survminer)
library(dplyr)

in_file="/Users/macbook/Desktop/APOP_PAN_CANCER/TCGA_raw.files/CESC-final.csv"
out_file="/Users/macbook/Desktop/APOP_PAN_CANCER/HR_APOP/CESC_HRmedian_raw.csv"
clin_feat=13

########################

LUSC_C <- read.csv(file=in_file,header =TRUE, sep = ",", dec = ".");
write.table(cbind("Gene","Beta","HR","P-value","GP1","GP2","Hr-Inv-lst","Concordance","Std_Error"),
            file=out_file,row.names=F,col.names=F,sep = ',');

for(i in seq(from=clin_feat+1, to=length(LUSC_C), by=1))#no of clinical features
{
  surv_object <- Surv(time = LUSC_C$OS.time, event = LUSC_C$vital_status)
  
  
  fit1 <- survfit(surv_object~(LUSC_C[,i])>(median(LUSC_C[,i])), data=LUSC_C); 
  summary(fit1);
  fit1.coxph <- coxph(surv_object ~ (LUSC_C[,i])>(median(LUSC_C[,i])), data = LUSC_C)
  # summary(fit1.coxph);
  first <- coef(summary(fit1.coxph))
  
  
  
  if((!is.na(first[5]))&&(!is.na(first[2])))
  {write.table(cbind(colnames(LUSC_C[i]),first[1],first[2],first[5],fit1$n[1],fit1$n[2],1/first[2],fit1.coxph$concordance[6],fit1.coxph$concordance[7]),
               file=out_file,row.names=F,col.names=F,sep = ',',append = T);#output file
  }
}
