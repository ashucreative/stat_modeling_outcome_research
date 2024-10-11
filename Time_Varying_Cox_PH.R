#install.packages("survsim")
library(survsim)
N=100 #number of patients
set.seed(123)

?simple.surv.sim

df.tf<-simple.surv.sim(#baseline time fixed
  n=N, 
  foltime=500, # maximum time of follow-up
  dist.ev=c('llogistic'), # time to event
  anc.ev=c(0.68), #ancillary parameter
  beta0.ev=c(5.8), #β0 parameter for the time to event distribution.
  anc.cens=1.2, # real number containing the ancillary parameter for the time to censoring distribution
  beta0.cens=7.4, # 	real number containing the β0 parameter for the time to censoring distribution
  
  z=list(c("unif", 0.8, 1.2)),
  beta=list(c(-0.4),c(0)),
  
  x=list(c("bern", 0.5), # distribution for grp
         c("normal", 70, 13)))#distribution for age 

names(df.tf)[c(1,6,7)]<-c("id","grp","age")


# Generate time-varying covariates

set.seed(123)
nft<-sample(1:10,
              N,replace=T)#number of follow up time points
crp<-round(abs(rnorm(sum(nft)+N,
                     mean=100,sd=40)),1)
time<-NA
id<-NA
i=0
for(n in nft){
  i=i+1
  time.n<-sample(1:500,n)
  time.n<-c(0,sort(time.n))
  time<-c(time,time.n)
  id.n<-rep(i,n+1)
  id<-c(id,id.n)
}

df.td <- cbind(data.frame(id,time)[-1,],crp)
# Merge and reshape the two datasets

library(survival)
df<-tmerge(df.tf,df.tf,id=id,
             endpt=event(stop,status))

head(round(df))

# split into the simulated time intervals
df <- tmerge(df,df.td,id=id,
             crp=tdc(time,crp))
head(round(df),10)

# Fit Cox model
fit.tdc <- coxph(Surv(tstart,tstop,endpt)~
                   grp+age+crp+cluster(id),df)
fit.tdc


# CPH with time-varying coefficients using lung dataset
fit2 <- coxph(Surv(time, status) ~
                age +ph.karno+sex,
              data=lung)
zph <- cox.zph(fit2)
zph

# Plot results
plot(zph[2],lwd=2)
abline(0,0, col=1,lty=3,lwd=2)
abline(h= fit2$coef[2], col=3, lwd=2, lty=2)
  legend("bottomright",
  legend=c('Reference line for null effect',
  "Average hazard over time",
  "Time-varying hazard"),
  lty=c(3,2,1), col=c(1,3,1), lwd=2)

# Two methods to include time-varying effect
# Step function to explore time-varying coefficient using survSplit function

lung.split <- survSplit(Surv(time, status) ~ .,
                          data= lung, cut=c(180, 350),
                          episode= "tgroup", id="id")
head(lung.split[-c(1,4,6:8)])

# Cox reg model stratified by the tgroup variable
fit.split <- coxph(Surv(tstart, time, status) ~
                     age + ph.karno:strata(tgroup)+
                     sex,
                   data=lung.split)
fit.split
# check assumption
cox.zph(fit.split)

# These results show that that there is no correlation
# between transformed survival time and the scaled
# Schoenfeld residuals.

# User-defined parametric continuous function to describe the time-varying coefficient
fit.tt <- coxph(Surv(time, status) ~
                  age + ph.karno + tt(ph.karno)+ sex,
                data=lung,
                tt = function(x, t, ...) x * log(t+20))
fit.tt

# Coefficients for both the ph.karno and tt(ph.karno) are significant, implying effect of ph.karno varies with time.
# Plot results
zph.tt <- cox.zph(fit2,
                  transform=function(t) log(t+20))
plot(zph.tt[2])
abline(coef(fit.tt)[2:3], col=2)


# Investigating time-varying coefficient with timereg package
install.packages("timereg")
library(timereg)
fit.out <- timecox(Surv(time,status)~
                       age+sex+ph.karno,
                     data=lung,n.sim=500,
                     max.time=700)
summary(fit.out)

par(mfrow=c(2,2))
plot(fit.out)

# Setting sex and age time-fixed variables
fit.const <- timecox(Surv(time,status)~
                       const(age)+const(sex)+ph.karno,
                     data=lung,n.sim=500,
                     max.time=700)
coef(fit.const)
