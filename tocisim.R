#data from appendix of distributions of LOSs reconstructed by using plot digitizer online; sampled evenly from numbers in stated ranges
library(dplyr)


ward<-data.frame(LOS=c(sample(c(1:4),1144,replace=T),sample(c(5:9),680,replace=T),sample(c(10:14),331,replace=T),sample(c(15:19),164,replace=T),sample(c(20:24),113,replace=T),sample(c(25:29),92,replace=T),sample(c(30:34),72,replace=T),sample(c(35:39),55,replace=T),sample(c(40:44),48,replace=T),sample(c(45:49),38,replace=T),sample(c(50:90),102,replace=T)))
library(fitdistrplus)
descdist(ward$LOS)
fitdist(ward$LOS,"gamma")
plot(density(ward$LOS))
0.9142624/0.0740240
1/0.0740240
mean(ward$LOS)
sd(ward$LOS)
#ward stay gamma=shape 0.96461310, scale=12.6665413

icu<-data.frame(LOS=c(sample(c(1:4),166,replace=T),sample(c(5:9),122,replace=T),sample(c(10:14),36,replace=T),sample(c(15:19),27,replace=T),sample(c(20:24),10,replace=T),sample(c(25:29),7,replace=T),sample(c(30:34),5,replace=T),sample(c(35:39),2,replace=T),sample(c(40:44),3,replace=T),sample(c(45:49),1,replace=T)))
fitdist(icu$LOS,"gamma")
#icu stay gamma=shape 1.3754566, scale=5.76

imv<-data.frame(LOS=c(sample(c(1:4),84,replace=T),sample(c(5:9),97,replace=T),sample(c(10:14),92,replace=T),sample(c(15:19),68,replace=T),sample(c(20:24),54,replace=T),sample(c(25:29),46,replace=T),sample(c(30:34),30,replace=T),sample(c(35:39),18,replace=T),sample(c(40:44),9,replace=T),sample(c(45:49),8,replace=T),sample(c(50:90),42,replace=T)))
fitdist(imv$LOS,"gamma")
#imv shape=1.30959774 scale=15.2632188

hist(icu$LOS)
hist(imv$LOS)

tracking<-data.frame(run=1:1000,total_chg=NA,avg_chg=NA,toci_cost=NA)

for (f in 1:1000){
  #probability of icu admission=.217 with sd=.0151
  sim<-data.frame(patient=1:5652,icuprob=rnorm(5652,mean=0.217,sd=0.0151))
  sim$icu_y_n<-NA
  for (i in 1:5652){
    sim$icu_y_n[i]<-sample(c(1,0),size=1,prob=c(sim$icuprob[i],1-sim$icuprob[i]))
  }
  sim$ivm<-NA
  for (i in 1:5652){
    sim$ivm[i]<-ifelse(sim$icu_y_n[i]==1,rnorm(1,mean=.434,sd=.033),0)
  }
  sim$ivm_y_n<-NA
  for (i in 1:5652){
    sim$ivm_y_n[i]<-ifelse(sim$ivm[i]>0,sample(c(1,0),size=1,prob=c(sim$ivm[i],1-sim$ivm[i])),0)
  }
  
  sim$hosp_only<-ifelse(sim$icu_y_n==0,1,0)
  sim$icu_no_ivm<-ifelse(sim$icu_y_n==1&sim$ivm_y_n==0,1,0)
  sim$icu_ivm<-ifelse(sim$icu_y_n==1&sim$ivm_y_n==1,1,0)
  sim$mort<-NA
  
  for (i in 1:5652){
    sim$mort[i]<-ifelse(sim$hosp_only[i]==1,rnorm(1,mean=.237,sd=.0069),ifelse(sim$icu_no_ivm[i]==1,rnorm(1,mean=.297,sd=.0223),rnorm(1,mean=.452,sd=.0198)))
  }
  
  sim$mort_y_n<-NA
  for (i in 1:5652){
    sim$mort_y_n[i]<-sample(c(1,0),size=1,prob=c(sim$mort[i],1-sim$mort[i]))
  }
  
  sim$los_hosp<-NA
  sim$los_icu<-NA
  hosp<-sim%>%filter(hosp_only==1)
  icu_no_ivm<-sim%>%filter(icu_no_ivm==1)
  icu_ivm<-sim%>%filter(icu_ivm==1)
  
  hosp$los_hosp<-ceiling(rgamma(nrow(hosp),shape=0.96461310,scale=12.6665413))
  hosp$los_icu<-0
  hosp$total_los<-hosp$los_hosp+hosp$los_icu
  mean(hosp$los_hosp)
  
  #58% of total LOS was in ICU, 42% on ward
  #scale=mean/shape
  icu_no_ivm$los_hosp<-ceiling(rgamma(nrow(icu_no_ivm),shape=1.3754566,scale=14.6/1.3754566))
  icu_no_ivm$los_icu<-icu_no_ivm$los_hosp*0.58
  icu_no_ivm$los_hosp<-icu_no_ivm$los_hosp-icu_no_ivm$los_icu
  icu_no_ivm$total_los<-icu_no_ivm$los_hosp+icu_no_ivm$los_icu
  mean(icu_no_ivm$total_los)
  
  #27% of total LOS on ward, 73% in ICU
  #scale=mean/shape
  icu_ivm$los_hosp<-ceiling(rgamma(nrow(icu_ivm),shape=1.30959774,scale=29.7/1.30959774))
  icu_ivm$los_icu<-icu_ivm$los_hosp*.73
  icu_ivm$los_hosp<-icu_ivm$los_hosp-icu_ivm$los_icu
  icu_ivm$total_los<-icu_ivm$los_hosp+icu_ivm$los_icu
  mean(icu_ivm$total_los)
  sd(icu_ivm$total_los)
  
  #averages all align closely with the averages stated in the paper, as do distributions
  sim<-rbind(hosp,icu_no_ivm,icu_ivm)
  sim$cost_of_stay<-sim$los_hosp*1104+sim$los_icu*5652
  
  #tocilizumab simulation
  icu_ivm$toci_risk<-NA
  icu_ivm$logodds.control<-log(icu_ivm$ivm/(1-icu_ivm$ivm))
  icu_ivm$rposterior<-df.ex$rposterior(n=nrow(icu_ivm),tau.sample=F)
  icu_ivm$logodds.treat<-icu_ivm$logodds.control+icu_ivm$rposterior
  icu_ivm$prob.treat<-exp(icu_ivm$logodds.treat)/(1+exp(icu_ivm$logodds.treat))
  icu_ivm$new_ivm<-NA
  for (i in 1:nrow(icu_ivm)){
    icu_ivm$new_ivm[i]<-sample(c(1,0),size=1,prob=c(icu_ivm$prob.treat[i],1-icu_ivm$prob.treat[i]))
  }
  new_icu_ivm<-subset(icu_ivm,new_ivm==1)
  new_icu_ivm$los_hosp<-ceiling(rgamma(nrow(new_icu_ivm),shape=1.30959774,scale=29.7/1.30959774))
  new_icu_ivm$los_icu<-new_icu_ivm$los_hosp*.73
  new_icu_ivm$los_hosp<-new_icu_ivm$los_hosp-new_icu_ivm$los_icu
  new_icu_ivm$total_los<-new_icu_ivm$los_hosp+new_icu_ivm$los_icu
  
  new_icu_no_ivm<-subset(icu_ivm,new_ivm==0)
  new_icu_no_ivm$los_hosp<-ceiling(rgamma(nrow(new_icu_no_ivm),shape=1.3754566,scale=14.6/1.3754566))
  new_icu_no_ivm$los_icu<-new_icu_no_ivm$los_hosp*0.58
  new_icu_no_ivm$los_hosp<-new_icu_no_ivm$los_hosp-new_icu_no_ivm$los_icu
  new_icu_no_ivm$total_los<-new_icu_no_ivm$los_hosp+new_icu_no_ivm$los_icu
  new_icu_ivm$cost_of_stay<-new_icu_ivm$los_hosp*1104+new_icu_ivm$los_icu*5652
  new_icu_no_ivm$cost_of_stay<-new_icu_no_ivm$los_hosp*1104+new_icu_no_ivm$los_icu*5652
  
  colnames(icu_ivm)
  allnew<-rbind(new_icu_ivm,new_icu_no_ivm)
  sim$stay_cost_toci<-ifelse(sim$icu_ivm==1,allnew$cost_of_stay[match(sim$patient,allnew$patient)],sim$cost_of_stay)
  
  sim$cost_of_toci<-ifelse(sim$icu_ivm==1,1800,0)
  sim$total_cost_toci<-sim$stay_cost_toci+sim$cost_of_toci
  median(sim$total_cost_toci)
  mean(sim$cost_of_stay)
  mean(sim$total_cost_toci)
  
  sim$change<-sim$total_cost_toci-sim$cost_of_stay
  
  tracking$total_chg[f]<-sum(sim$change)
  tracking$avg_chg[f]<-mean(sim$change)
  tracking$toci_cost[f]<-sum(sim$cost_of_toci)}
mean(tracking$total_chg)
mean(tracking$toci_cost)
quantile(tracking$total_chg,c(0.025,0.975))
quantile(tracking$toci_cost,c(0.025,0.975))

#mort using pooled effects from all trials
trackmort<-data.frame(run=1:1000,prev_mort=NA,new_mort=NA,diff=NA)
mort<-icu_ivm
mort$logodds.control<-log(mort$mort/(1-mort$mort))
for (f in 1:1000){
  mort$rposterior<-df.ex$rposterior(n=nrow(mort),tau.sample=F)
  mort$logodds.treat<-mort$logodds.control+mort$rposterior
  mort$prob.treat<-exp(mort$logodds.treat)/(1+exp(mort$logodds.treat))
  mort$new_mort<-NA
  for (i in 1:nrow(mort)){
    mort$new_mort[i]<-sample(c(1,0),size=1,prob=c(mort$prob.treat[i],1-mort$prob.treat[i]))
  }
  trackmort$prev_mort[f]<-sum(mort$mort)
  trackmort$new_mort[f]<-sum(mort$new_mort)
  trackmort$diff[f]<-trackmort$new_mort[f]-trackmort$prev_mort[f]
}
mean(trackmort$diff)
quantile(trackmort$diff,c(0.025,0.975))

#mort using REMAP-CAP only, with largest variation; slightly skewed distribution
trackmort<-data.frame(run=1:1000,prev_mort=NA,new_mort=NA,diff=NA)

for (f in 1:1000){
  mort$rposterior<-rnorm(n=nrow(mort),mean=0.91,sd=(1.18-0.91)/1.96)
  mort$prob.treat<-mort$mort*mort$rposterior
  mort$new_mort<-NA
  for (i in 1:nrow(mort)){
    mort$new_mort[i]<-sample(c(1,0),size=1,prob=c(mort$prob.treat[i],1-mort$prob.treat[i]))
  }
  trackmort$prev_mort[f]<-sum(mort$mort)
  trackmort$new_mort[f]<-sum(mort$new_mort)
  trackmort$diff[f]<-trackmort$new_mort[f]-trackmort$prev_mort[f]
}
mean(trackmort$diff)
quantile(trackmort$diff,c(0.025,0.975))
