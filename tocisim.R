#individual level data provided in BMJ Veiga paper; 
ordinal<-c(rep(4,67),rep(5,41),rep(6,21))
toci<-c(rep(1,39),rep(0,28),rep(1,15),rep(0,26),rep(1,11),rep(0,10))
death<-c(rep(1,3),rep(0,36),rep(1,3),rep(0,25),rep(1,3),rep(0,12),rep(0,26),rep(1,8),rep(0,3),rep(1,3),rep(0,7))
data<-data.frame(ordinal,toci,death)
data$death<-factor(data$death)
data$toci<-factor(data$toci)
logit<-glm(death~toci+ordinal,data=data,family="binomial")
summary(logit)
exp(cbind(OR = coef(logit), confint(logit)))
pred<-data.frame(ordinal=c(rep(5,762-521),rep(6,521)),toci=1)

pred$toci_y<-predict(logit,pred,type="response")
pred$toci<-0
pred$toci_n<-predict(logit,pred,type="response")

#severity distribution
trials<-data.frame(trial=c(rep("COVACTA",4),rep("EMPACTA",4),rep("BACC",4),rep("CORIMUNO",4),rep("REMAP",4),rep("Veiga",4)),severity=rep(3:6,6),percent=c(3.4,27.9,30.4,38.3,9.3,64.2,26.5,0,15.6,79.8,4.1,0.5,0,100,0,0,0,0.4,69.8,29.8,0,51.9,31.8,16.3))
ggplot(trials,aes(x=severity,y=percent,fill=-severity))+geom_bar(stat="identity")+facet_wrap(~trial)+ggtitle("Tocilizumab Trials by Percentage of Patients in Baseline Severity Categories")+labs(x="Severity",y="Percent of Patients in Category")+theme(legend.position = "none",plot.title=element_text(hjust=0.5))

#wilcoxon rank sum test for difference in severity distributions
study<-c(rep("REMAP",755),rep("Veiga",129))
rank<-c(rep(4,3),rep(5,527),rep(6,225),rep(4,67),rep(5,41),rep(6,21))
boxplot(rank~study)
wilcox.test(rank~study)
mean(rank[study=="REMAP"])
mean(rank[study=="Veiga"])
755*129
72876/97395
#w-stat=72876;75% of REMAP-Veiga pairings will have REMAP patient more severe than Veiga


#data from appendix of distributions of LOSs reconstructed by using plot digitizer online; sampled evenly from numbers in stated ranges
library(dplyr)


ward<-data.frame(LOS=c(sample(c(1:4),1144,replace=T),sample(c(5:9),680,replace=T),sample(c(10:14),331,replace=T),sample(c(15:19),164,replace=T),sample(c(20:24),113,replace=T),sample(c(25:29),92,replace=T),sample(c(30:34),72,replace=T),sample(c(35:39),55,replace=T),sample(c(40:44),48,replace=T),sample(c(45:49),38,replace=T),sample(c(50:90),102,replace=T)))
library(fitdistrplus)
descdist(ward$LOS)
fitdist(ward$LOS,"gamma")
plot(density(ward$LOS))
0.93560655/0.07642482
1/0.07642482
mean(ward$LOS)
sd(ward$LOS)
#ward stay gamma=shape 0.93560655, scale=12.23952

icu<-data.frame(LOS=c(sample(c(1:4),166,replace=T),sample(c(5:9),122,replace=T),sample(c(10:14),36,replace=T),sample(c(15:19),27,replace=T),sample(c(20:24),10,replace=T),sample(c(25:29),7,replace=T),sample(c(30:34),5,replace=T),sample(c(35:39),2,replace=T),sample(c(40:44),3,replace=T),sample(c(45:49),1,replace=T)))
fitdist(icu$LOS,"gamma")
#icu stay gamma=shape 1.3103809, scale=1/0.1661165

imv<-data.frame(LOS=c(sample(c(1:4),84,replace=T),sample(c(5:9),97,replace=T),sample(c(10:14),92,replace=T),sample(c(15:19),68,replace=T),sample(c(20:24),54,replace=T),sample(c(25:29),46,replace=T),sample(c(30:34),30,replace=T),sample(c(35:39),18,replace=T),sample(c(40:44),9,replace=T),sample(c(45:49),8,replace=T),sample(c(50:90),42,replace=T)))
fitdist(imv$LOS,"gamma")
#imv shape=1.30491738 scale=1/0.06521814

hist(icu$LOS)
hist(imv$LOS)

#no tocilizumab
tracking<-data.frame(run=1:10000,in_hosp=NA,icu_noivm=NA,icu_ivm=NA,hosp_mort=NA,no_ivm_mort=NA,ivm_mort=NA,mort=NA)

for (f in 1:10000){
  icu_prob<-rnorm(1,mean=0.162,sd=0.0049)
  icu<-ceiling(icu_prob*5701)
  hosp<-5701-icu
  icu_ivm_prob<-rnorm(1,mean=0.7,sd=0.0151)
  icu_ivm<-ceiling(icu_ivm_prob*icu)
  icu_no_ivm<-icu-icu_ivm
  hosp_mort_prob<-rnorm(1,mean=.129,sd=.0049)
  hosp_mort<-ceiling(hosp*hosp_mort_prob)
  icu_ivm_mort_prob<-rnorm(1,mean=.292,sd=.0179)
  icu_ivm_mort<-ceiling(icu_ivm_mort_prob*icu_ivm)
  icu_no_ivm_mort_prob<-rnorm(1,mean=.195,sd=.0238)
  icu_no_ivm_mort<-ceiling(icu_no_ivm_mort_prob*icu_no_ivm)
  totalmort<-hosp_mort+icu_ivm_mort+icu_no_ivm_mort
  tracking$in_hosp[f]<-hosp
  tracking$icu_noivm[f]<-icu_no_ivm
  tracking$icu_ivm[f]<-icu_ivm
  tracking$hosp_mort[f]<-hosp_mort
  tracking$no_ivm_mort[f]<-icu_no_ivm_mort
  tracking$ivm_mort[f]<-icu_ivm_mort
  tracking$mort[f]<-totalmort
}
tracking$hosp_los<-rgamma(10000,shape=0.96461310,scale=12.6665413)
tracking$no_ivm_los<-rgamma(10000,shape=1.3754566,scale=14.6/1.3754566)
tracking$ivm_los<-rgamma(10000,shape=1.30959774,scale=29.7/1.30959774)
tracking$tot_hosp_los<-tracking$hosp_los*tracking$in_hosp
tracking$tot_no_ivm_los<-tracking$no_ivm_los*tracking$icu_noivm
tracking$tot_ivm_los<-tracking$ivm_los*tracking$icu_ivm
tracking$hosp_cost<-tracking$tot_hosp_los*1104
tracking$no_ivm_cost<-(tracking$tot_no_ivm_los*.42*1104)+(tracking$tot_no_ivm_los*.58*5652)
tracking$ivm_cost<-(tracking$tot_ivm_los*.27*1104)+(tracking$tot_ivm_los*.73*5652)
tracking$cost<-tracking$hosp_cost+tracking$no_ivm_cost+tracking$ivm_cost
mean(tracking$cost)

tracktoci<-tracking
tracktoci$icu_ivm_prob<-tracktoci$icu_ivm/(tracktoci$icu_noivm+tracktoci$icu_ivm)
tracktoci$icu_ivm_cont<-log(tracktoci$icu_ivm_prob/(1-tracktoci$icu_ivm_prob))
tracktoci$icu_ivm_treat<-(tracktoci$icu_ivm_cont+df.mv$rposterior(n=10000,tau.sample = F))
tracktoci$icu_ivm_prob<-exp(tracktoci$icu_ivm_treat)/(1+exp(tracktoci$icu_ivm_treat))
tracktoci$icu_total<-tracktoci$icu_noivm+tracktoci$icu_ivm
tracktoci$icu_ivm<-ceiling(tracktoci$icu_total*tracktoci$icu_ivm_prob)
tracktoci$icu_noivm<-tracktoci$icu_total-tracktoci$icu_ivm

tracktoci$hosp_mort_prob_old<-tracktoci$hosp_mort/tracktoci$in_hosp
tracktoci$hosp_mort_cont<-log(tracktoci$hosp_mort_prob_old/(1-tracktoci$hosp_mort_prob_old))
tracktoci$hosp_mort_treat<-(tracktoci$hosp_mort_cont+df.mort$rposterior(n=10000,tau.sample=F))
tracktoci$hosp_mort_prob<-exp(tracktoci$hosp_mort_treat)/(1+exp(tracktoci$hosp_mort_treat))
tracktoci$hosp_mort<-ceiling(tracktoci$hosp_mort_prob*tracktoci$in_hosp)

tracktoci$no_ivm_prob_old<-tracktoci$no_ivm_mort/tracktoci$icu_noivm
tracktoci$no_ivm_mort_cont<-log(tracktoci$no_ivm_prob_old/(1-tracktoci$no_ivm_prob_old))
tracktoci$no_ivm_mort_treat<-(tracktoci$no_ivm_mort_cont+df.mort$rposterior(n=10000,tau.sample=F))
tracktoci$no_ivm_mort_prob<-exp(tracktoci$no_ivm_mort_treat)/(1+exp(tracktoci$no_ivm_mort_treat))
tracktoci$no_ivm_mort<-ceiling(tracktoci$no_ivm_mort_prob*tracktoci$icu_noivm)

tracktoci$ivm_prob_old<-tracktoci$ivm_mort/tracktoci$icu_ivm
tracktoci$ivm_mort_cont<-log(tracktoci$ivm_prob_old/(1-tracktoci$ivm_prob_old))
tracktoci$ivm_mort_treat<-(tracktoci$ivm_mort_cont+df.mort$rposterior(n=10000,tau.sample=F))
tracktoci$ivm_mort_prob<-exp(tracktoci$ivm_mort_treat)/(1+exp(tracktoci$ivm_mort_treat))
tracktoci$ivm_mort<-ceiling(tracktoci$ivm_mort_prob*tracktoci$icu_ivm)
tracktoci$mort<-tracktoci$hosp_mort+tracktoci$no_ivm_mort+tracktoci$ivm_mort

tracktoci$toci_doses<-rnorm(10000,.288,.0113)
tracktoci$toci_cost<-2793*tracktoci$toci_doses*3600+2793*(1-tracktoci$toci_doses)*1800

tracktoci$tot_hosp_los<-tracktoci$hosp_los*tracktoci$in_hosp
tracktoci$tot_no_ivm_los<-tracktoci$no_ivm_los*tracktoci$icu_noivm
tracktoci$tot_ivm_los<-tracktoci$ivm_los*tracktoci$icu_ivm
tracktoci$hosp_cost<-tracktoci$tot_hosp_los*1104
tracktoci$no_ivm_cost<-(tracktoci$tot_no_ivm_los*.42*1104)+(tracktoci$tot_no_ivm_los*.58*5652)
tracktoci$ivm_cost<-(tracktoci$tot_ivm_los*.27*1104)+(tracktoci$tot_ivm_los*.73*5652)
tracktoci$cost<-tracktoci$hosp_cost+tracktoci$no_ivm_cost+tracktoci$ivm_cost+tracktoci$toci_cost

compare<-data.frame(run=1:10000,in_hosp=tracking$in_hosp,diff_icu_ivm=tracktoci$icu_ivm-tracking$icu_ivm,diff_mort=tracktoci$mort-tracking$mort,diff_cost=tracktoci$cost-tracking$cost,tocicost=tracktoci$toci_cost)
mean(compare$diff_cost)
quantile(compare$diff_cost,c(0.025,0.975))

#sensitivity analysis
#only 1 dose of toci
tracktoci$toci_cost<-2793*1800
tracktoci$cost<-tracktoci$hosp_cost+tracktoci$no_ivm_cost+tracktoci$ivm_cost+tracktoci$toci_cost
mean(tracktoci$cost-tracking$cost)
quantile(tracktoci$cost-tracking$cost,c(0.025,0.975))

#prevented mech vent patients go to hosp_only, not icu noivm
tracktoci$mech_vent_prev<-tracking$icu_ivm-tracktoci$icu_ivm
tracktoci$in_hosp<-tracktoci$in_hosp+tracktoci$mech_vent_prev
tracktoci$icu_noivm<-5701-tracktoci$in_hosp-tracktoci$icu_ivm
tracktoci$tot_hosp_los<-tracktoci$hosp_los*tracktoci$in_hosp
tracktoci$tot_no_ivm_los<-tracktoci$no_ivm_los*tracktoci$icu_noivm
tracktoci$tot_ivm_los<-tracktoci$ivm_los*tracktoci$icu_ivm
tracktoci$hosp_cost<-tracktoci$tot_hosp_los*1104
tracktoci$no_ivm_cost<-(tracktoci$tot_no_ivm_los*.42*1104)+(tracktoci$tot_no_ivm_los*.58*5652)
tracktoci$ivm_cost<-(tracktoci$tot_ivm_los*.27*1104)+(tracktoci$tot_ivm_los*.73*5652)

tracktoci$toci_cost<-2793*tracktoci$toci_doses*3600+2793*(1-tracktoci$toci_doses)*1800
tracktoci$cost<-tracktoci$hosp_cost+tracktoci$no_ivm_cost+tracktoci$ivm_cost+tracktoci$toci_cost
mean(tracktoci$cost)
mean(tracktoci$cost-tracking$cost)
