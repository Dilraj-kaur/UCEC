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



beta_file='/Users/dilrajkaur/Desktop/UCEC/UCEC_HRmedian_raw.csv'
#tcga_datafile="/Users/dilrajkaur/Desktop/UCEC/UCEC.RDS"
tcga_datafile2="/Users/dilrajkaur/Desktop/UCEC/UCEC-final.csv"

beta<-read.csv(file=beta_file,header =TRUE, sep = ",", dec = ".")
data1<-read.csv(file=tcga_datafile2,header =TRUE, sep = ",", dec = ".")
#data1<-readRDS(tcga_datafile)
genes <- scan(file = '///Users/dilrajkaur/Desktop/UCEC/newdata/input_genes/UCEC_genes4.csv', what = 'character', sep = ',')

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
data2=data1
######################################## clinical ############################

#table(data1$clinical_stage)


#data2<- subset(data1, ((race != '[NA]')&(race != '[Unknown]')&(race != '[Not Evaluated]')))

#example for stage 3,4 vs stage 1,2
surv_object <- Surv(time = data2$OS.time/30, event = data2$vital_status)


#survival analysis: fits cox ph model to find HR for mean cut
fit1 <- survfit(surv_object~((data2$clinical_stage!='Stage I')&(data2$clinical_stage!='Stage IA')&(data2$clinical_stage!='Stage IB')&(data2$clinical_stage!='Stage IC')&(data2$clinical_stage!='Stage II')&(data2$clinical_stage!='Stage IIB')))
#summary(fit1);
ggsurvplot(fit1, data=data2,xlab = "Time (Months)")


fit1.coxph <- coxph(surv_object~((data2$clinical_stage!='Stage I')&(data2$clinical_stage!='Stage IA')&(data2$clinical_stage!='Stage IB')&(data2$clinical_stage!='Stage IC')&(data2$clinical_stage!='Stage II')&(data2$clinical_stage!='Stage IIB')))

# summary(fit1.coxph);
first <- coef(summary(fit1.coxph))
#first
a=summary(fit1.coxph)
a$concordance[1]

#no of samples,HR,p-val,C,%95 CI L,%95CI U,logrank-p
print(c(nrow(data2),first[2],first[5],a$concordance[1],a$conf.int[3],a$conf.int[4],a$logtest[3]))

gg<-ggsurvplot(fit1, data=data2,xlab = "Time (Months)",pval = TRUE,legend.labs = c("Low risk", "High risk "),axes.offset=FALSE)
gg
ggsave('/Users/dilrajkaur/Desktop/UCEC/images/stage3_4.jpg', plot = print(gg),height = 5 ,width = 7, dpi = 1000)



##################################### subgroup- classification #########################
data2=data1
table(data2$menopause_status)
data2=subset(data1,((menopause_status=="Peri (6-12 months since last menstrual period)")|(menopause_status=="Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)")))

surv_object <- Surv(time = data2$OS.time/30, event = data2$vital_status)

#survival analysis: fits cox ph model to find HR for mean cut
fit1 <- survfit(surv_object~data2$PI>2); 
#summary(fit1);
ggsurvplot(fit1, data=data2)


fit1.coxph <- coxph(surv_object~data2$PI>2)
# summary(fit1.coxph);
first <- coef(summary(fit1.coxph))
#first
a=summary(fit1.coxph)

print(c(nrow(data2),first[2],first[5],a$concordance[1],a$conf.int[3],a$conf.int[4],a$logtest[3]))

gg<-ggsurvplot(fit1, data=data2,xlab = "Time (Months)",pval = TRUE,legend.labs = c("Low risk", "High risk"),axes.offset=FALSE)
gg
ggsave('/Users/dilrajkaur/Desktop/UCEC/newdata/images/menopause.subgroup.jpg', plot = print(gg),height = 5 ,width = 7, dpi = 1000)

###################### multivariate ##########################
data2=data1
#data2<-subset(data1,(ethnicity!='[NA]' )& (ethnicity!='[Not Evaluated]'))
data2<-subset(data2,(residual_tumor!='[NA]')&(residual_tumor!='[Unknown]'))
#data2<-subset(data2,(tumor_status!='[NA]')&(tumor_status!='[Unknown]'))
data2<-subset(data2,(menopause_status!='[NA]')&(menopause_status!='[Unknown]'))
data2<-subset(data2,(peritoneal_washing!='[NA]')&(peritoneal_washing!='[Not Evaluated]'))

data2 <- mutate(data2, PI = ifelse((PI > 2), "2", "1"))


data2 <- mutate(data2, clinical_stage = ifelse(((clinical_stage=='Stage III')|(clinical_stage=='Stage IIIA')|(clinical_stage=='Stage IIIB')|(clinical_stage=='Stage IIIC')|(clinical_stage=='Stage IIIC1')|(clinical_stage=='Stage IIIC2')|(clinical_stage=='Stage IV')|(clinical_stage=='Stage IVA')|(clinical_stage=='Stage IVB')), "2", "1"))
data2 <- mutate(data2, residual_tumor = ifelse(((residual_tumor=='R1')|(residual_tumor=='R2')), "2", "1"))
data2<- mutate(data2, histologic_diagnosis=ifelse((histologic_diagnosis=='Mixed serous and endometrioid')|(histologic_diagnosis=='Serous endometrial adenocarcinoma'),"2" ,"1"))
data2<-mutate(data2,menopause_status= ifelse(((menopause_status=='Peri(6-12 months since last menstrual period)')|(menopause_status=='Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)')),"2","1"))
#data2 <- mutate(data2, tumor_status = ifelse(((tumor_status=='WITH TUMOR')), "2", "1"))
data2 <- mutate(data2, peritoneal_washing = ifelse(((peritoneal_washing=='positive')), "2", "1"))
data2 <- mutate(data2, neoplasm_histologic_grade= ifelse(((neoplasm_histologic_grade=='G3')|(neoplasm_histologic_grade=='High Grade')), "2", "1"))



data2$clinical_stage<-factor(data2$clinical_stage,labels=c("Stage1,2","Stage3,4"),levels=c("1","2"))
data2$histologic_diagnosis<-factor(data2$histologic_diagnosis,labels =c("EEC","SEA,MSE"),levels=c("1","2"))
data2$residual_tumor<-factor(data2$residual_tumor,labels=c("R0","R1,R2"),levels = c("1","2"))
data2$menopause_status<-factor(data2$menopause_status,labels =c("Post","Pre"),levels = c("1","2"))
data2$neoplasm_histologic_grade<-factor(data2$neoplasm_histologic_grade,labels =c("Low grade","High grade"),levels = c("1","2"))
data2$peritoneal_washing<-factor(data2$peritoneal_washing,labels=c("negative","positive"),levels = c("1","2"))
#data2$tumor_status<-factor(data2$tumor_status,labels=c("Tumor free","With tumor"),levels = c("1","2"))
data2$Gene_voting_model<-factor(data2$PI,labels=c("Low risk","High risk"),levels = c("1","2"))


surv_object <- Surv(time = data2$OS.time/30, event = data2$vital_status)
fit.coxph <- coxph(surv_object ~ Gene_voting_model+clinical_stage+histologic_diagnosis+residual_tumor+menopause_status+neoplasm_histologic_grade+peritoneal_washing, 
                   data = data2)
gg<-ggforest(fit.coxph, data = data2)
ggsave('/Users/dilrajkaur/Desktop/UCEC/newdata/images/4.genes.final.multivariate.jpg', plot = print(gg),height = 6 ,width = 8, dpi = 1000)



###################### KM plot options #########################

ggsurv<-ggsurvplot(
  fit1,                            # survfit object with calculated statistics.
  data = data1,                   # data used to fit survival curves.
  risk.table = TRUE,               # show risk table.
  #title=" \t \t BCL2+BCLXL-BAX-BAK ",
  conf.int = FALSE,                 # show confidence intervals for
  # palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,150),                  # present narrower X axis, but not affect
  #ylim=c(0.1,1),
  xlab = "Time (Months)",         # customize X axis label.
  # break.time.by = 20,             # break X axis in time intervals by 500.
  ggtheme = theme_light(),         # customize plot and risk table with a theme.
  risk.table.height = 0.25,         # the height of the risk table
  risk.table.y.text = FALSE,        # show bars instead of names in text annotations
  # # in legend of risk table.
  conf.int.style = "step",         # customize style of confidence intervals
  legend=c(0.8,0.8),
  legend.labs = c("High Risk", "Low Risk"),  # change legend labels.
  # ncensor.plot = TRUE,               # plot the number of censored subjects at time t
  # ncensor.plot.height = 0.25,
  # pval = TRUE,                       # show p-value of log-rank test.
  axes.offset=FALSE,
)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("label", x = 30, y = 0.25, # x and y coordinates of the text
                    label = "HR = 2.94, p-val= 6.54e-10 \n 95%CI(2.09 - 4.15)
Wald p= 7e-10 & logrank, p= 9e-11", size = 4)
ggsurv


#fit1.coxph$concordance[6]