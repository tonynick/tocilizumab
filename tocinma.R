#do for mechanical ventilation
data<-read.csv("C:/Users/tonynickonchuk/Desktop/TociBamlan/tocimort.csv",stringsAsFactors = F)
data<-read.csv("C:/Users/tonynickonchuk/Desktop/TociBamlan/tocimechvent.csv",stringsAsFactors = F)

network<-mtc.network(data)
plot(network)
#use for mortality
model<-mtc.model(network,type="consistency",likelihood="binom",
                 link="logit",linearModel="fixed",n.chain=3,
                 powerAdjust=NA,dic=T,hy.prior=mtc.hy.empirical.lor("mortality","pharma-control"))
#use for mech vent
model<-mtc.model(network,type="consistency",likelihood="binom",
                 link="logit",linearModel="fixed",n.chain=3,
                 powerAdjust=NA,dic=T,hy.prior=mtc.hy.empirical.lor("semi-objective","pharma-control"))
results<-mtc.run(model,n.adapt=500,n.iter=20000,thin=1)
gemtc::forest(relative.effect(results,"standardCare"),digits=2,use.description=T)
rank.prob<-rank.probability(results,preferredDirection = -1)
plot(rank.prob)
mtcresults = as.data.frame(round(relative.effect.table(results),2))