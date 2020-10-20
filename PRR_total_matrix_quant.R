library(survival)
library(survminer)
library(dplyr)

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <<- apply(df_sorted, 1, mean)
  
  index_to_mean <<- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

LUSC_C <- read.csv(file='/Users/dilrajkaur/Desktop/TCGA-survival/LIHC/LIHC-final.csv',header =TRUE, sep = ",", dec = ".");

LUNG1 <- LUSC_C[,15:length(LUSC_C)]
LUNG1 <- LUNG1[, which(as.numeric(colSums(LUNG1!= 0)) > (nrow(LUNG1)/2))] # removing genes with 0 values in more than 50 % data.
LUSC_C<-cbind(LUSC_C[,1:14],LUNG1)
LUNG3q <- quantile_normalisation(t(LUSC_C[,15:length(LUSC_C)]))
LUNG3q<- as.data.frame(LUNG3q)
LUNG3f<-as.data.frame(t(LUNG3q))
LUSC_C[,15:length(LUSC_C)]<- LUNG3f
LUSC_C<-subset(LUSC_C,vital_status!="[NA]")
LUSC_C$vital_status<-as.numeric(ifelse(LUSC_C$vital_status=="0",0,1))

target<- as.data.frame.numeric(df_mean)    # For storing the mean values of the genes
rownames(target)<-colnames(LUNG3f)
write.table(target,file="/Users/dilrajkaur/Desktop/TCGA-survival/LIHC/LIHC-target_quantile.csv",row.names=T,col.names=F,sep = ',');
write.table(LUSC_C,file="//Users/dilrajkaur/Desktop/TCGA-survival/LIHC/LIHC-final_quantile.csv",row.names=F,col.names=T,sep = ',');
# LUSC_C[,19:length(LUSC_C)] <- quantile_normalisation(LUSC_C[,19:length(LUSC_C)])


write.table(cbind("Gene","Beta","HR","P-value","GP1","GP2","Hr-Inv-lst","Concordance","Std_Error"),
            file="/Users/dilrajkaur/Desktop/TCGA-survival/LIHC/LIHC-TOTAL-final-HRmedian-quantright.csv",row.names=F,col.names=F,sep = ',');

for(i in seq(from=15, to=length(LUSC_C), by=1))
{
  surv_object <- Surv(time = LUSC_C$OS.time, event = LUSC_C$vital_status)
  
  #survival analysis: fits cox ph model to find HR for mean cut
  fit1 <- survfit(surv_object~(LUSC_C[,i])>(median(LUSC_C[,i])), data=LUSC_C); 
  summary(fit1);
  fit1.coxph <- coxph(surv_object ~ (LUSC_C[,i])>(median(LUSC_C[,i])), data = LUSC_C)
  # summary(fit1.coxph);
  first <- coef(summary(fit1.coxph))
  #as.matrix(first)
  
  #check whether the pvalue is significant and HR is more than 20 (only bclxl had HR=20.8)
  if((!is.na(first[5]))&&(!is.na(first[2])))
  {write.table(cbind(colnames(LUSC_C[i]),first[1],first[2],first[5],fit1$n[1],fit1$n[2],1/first[2],fit1.coxph$concordance[6],fit1.coxph$concordance[7]),
               file="/Users/dilrajkaur/Desktop/TCGA-survival/LIHC/LIHC-TOTAL-final-HRmedian-quantright.csv",row.names=F,col.names=F,sep = ',',append = T);#output file
  }
}

