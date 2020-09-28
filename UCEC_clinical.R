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



beta_file='//Users/dilrajkaur/Desktop/UCEC/UCEC_HRmedian_raw.csv'
#tcga_datafile="/Users/dilrajkaur/Desktop/UCEC/UCEC.RDS"

beta<-read.csv(file=beta_file,header =TRUE, sep = ",", dec = ".")
data1<-read.csv(file='/Users/dilrajkaur/Desktop/UCEC/UCEC-final.csv', header = TRUE,sep=",")
genes <- scan(file = '/Users/dilrajkaur/Desktop/UCEC/UCEC_genes10.csv', what = 'character', sep = ',')

PI=0
v=0
for(i in seq(from=1, to=length(genes), by=1))
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

table(data2$tumor_status)


data2<- subset(data1, (pregnancies_full_term_count!= '[NA]'& pregnancies_full_term_count!='[Unknown]'))

#example for stage 3,4 vs stage 1,2
surv_object <- Surv(time = data2$OS.time/30, event = data2$vital_status)


#survival analysis: fits cox ph model to find HR for mean cut
fit1 <- survfit(surv_object~data2$pregnancies_full_term_count=='4+')

#fit1 <- survfit(surv_object~((data2$menopause_status!='Stage I')&(data2$clinical_stage!='Stage IA')&(data2$clinical_stage!='Stage IB')&(data2$clinical_stage!='Stage IC')&(data2$clinical_stage!='Stage II')&(data2$clinical_stage!='Stage IIB')))
#summary(fit1);
ggsurvplot(fit1, data=data2,xlab = "Time (Months)")

fit1.coxph<-coxph(surv_object~data2$pregnancies_full_term_count=='4+')
#fit1.coxph <- coxph(surv_object~((data2$clinical_stage!='Stage I')&(data2$clinical_stage!='Stage IA')&(data2$clinical_stage!='Stage IB')&(data2$clinical_stage!='Stage IC')&(data2$clinical_stage!='Stage II')&(data2$clinical_stage!='Stage IIB')))

# summary(fit1.coxph);
first <- coef(summary(fit1.coxph))
#first
a=summary(fit1.coxph)
a$concordance[1]

#no of samples,HR,p-val,C,%95 CI L,%95CI U,logrank-p
print(c(nrow(data2),first[2],first[5],a$concordance[1],a$conf.int[3],a$conf.int[4],a$logtest[3]))

gg<-ggsurvplot(fit1, data=data2,xlab = "Time (Months)",pval = TRUE,legend.labs = c("NOT HISPANIC OR LATINO", "HISPANIC OR LATINO"),axes.offset=FALSE)
gg
ggsave('/Users/dilrajkaur/Desktop/UCEC/ETHNICITY.jpg', plot = print(gg),height = 5 ,width = 7, dpi = 1000)


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