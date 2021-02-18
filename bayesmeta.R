library(bayesmeta)
library(metafor)
#mortality tocilizumab
df<-data.frame(Author=c("Salvarini","COVACTA","EMPACTA","BACC","CORIMUNO","REMAP-CAP","Veiga","RECOVERY"),
               Ee=c(1,58,26,9,7,98,14,596),Ne=c(63,294,249,161,63,350,65,2022),
               Ec=c(1,28,11,3,8,142,6,694),Nc=c(60,144,128,82,67,397,64,2094))

df.es<-escalc(measure="RR",ai=Ee,n1i=Ne,ci=Ec,n2i=Nc,slab=Author,data=df)

td<-TurnerEtAlPrior(outcome="all-cause mortality",comparator1="pharmacological",comparator2="placebo / control")

df.ex<-bayesmeta(y=df.es$yi,
                           sigma=sqrt(df.es$vi),
                           label=df.es$Author,mu.prior.mean=0,mu.prior.sd=4,
                 tau.prior = td$dprior)

forestplot(df.ex,shrinkage = F)

#probability of mean<0
1-(1-df.ex$pposterior(mu=0))
#mean mortality rate of hospitalized patients is 20%
#https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2773971
prob.control<-0.2
logodds.control<-log(prob.control/(1-prob.control))
logodds.treat<-(logodds.control+df.ex$rposterior(n=10000,tau.sample=F))
prob.treat<-exp(logodds.treat)/(1+exp(logodds.treat))
riskdiff<-(prob.treat-prob.control)
nnt<-1/riskdiff
quantile(riskdiff,c(0.025,0.5,0.975))
1/quantile(riskdiff,c(0.025,0.5,0.975))

#mechvent tocilizumab
df<-data.frame(Author=c("COVACTA","EMPACTA","BACC","CORIMUNO","Veiga","Salvarani","RECOVERY"),
               Ee=c(51,20,11,5,7,6,215),Ne=c(183,249,161,63,65,60,1754),
               Ec=c(33,16,8,14,11,5,273),Nc=c(90,128,81,67,64,63,1800))
td<-TurnerEtAlPrior(outcome="cause-specific mortality / major morbidity event / composite (mortality or morbidity)",comparator1="pharmacological",comparator2="placebo / control")

df.es<-escalc(measure="RR",ai=Ee,n1i=Ne,ci=Ec,n2i=Nc,slab=Author,data=df)

df.ex<-bayesmeta(y=df.es$yi,
                           sigma=sqrt(df.es$vi),
                           label=df.es$Author,mu.prior.mean=0,mu.prior.sd=4,
                 tau.prior = td$dprior)

forestplot(df.ex,shrinkage = F)

#probability of mean<0
1-(1-df.ex$pposterior(mu=0))

#https://journals.lww.com/ccejournal/fulltext/2020/10000/timing_of_intubation_and_in_hospital_mortality_in.45.aspx
prob.control<-0.129
logodds.control<-log(prob.control/(1-prob.control))
logodds.treat<-(logodds.control+df.ex$rposterior(n=10000,tau.sample=F))
prob.treat<-exp(logodds.treat)/(1+exp(logodds.treat))
riskdiff<-(prob.treat-prob.control)

quantile(riskdiff,c(0.025,0.5,0.975))
1/quantile(riskdiff,c(0.025,0.5,0.975))
