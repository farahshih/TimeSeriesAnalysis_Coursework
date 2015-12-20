setwd("C:\\Farahshih\\SFO\\Study at Berkeley\\STAT 153\\Project")
library(TSA)

####################### Read data ######################
mydata<-read.csv("mydata.csv")
str(mydata)
date<-seq(as.Date("1996/1/1"),as.Date("2015/04/1"),by="month")
mydata$month<-strftime(date,"%Y/%m")
mydata$month


########### Explore Data & Prewhiten ###########
visit<-mydata$us_source
visit<-ts(visit,start=c(1996,1),end=c(2015,4),frequency=12)

stl_visit<-stl(visit,s.window = "periodic") ## earlier time fluctuate more violently
plot(stl_visit, main="Seasonal Decomposition of Time Series by Loess") #use loess

windows(pointsize=12)
Month<-c('J','F','M','A','M','J','J','A','S','O','N','D')
plot(visit,ylab='Number of visitors',xlab='Time', 
     main="Monthly Taiwanese visitors to the US")
points(visit,pch=Month)
abline(h=mean(visit))


######### Examine the relationship between CPI and visitors ##############
cpi<-mydata$CPI
cpi<-ts(cpi,start=c(1996,1),end=c(2015,4),frequency=12)
summary(cpi)
summary(visit)

plot(cpi,ylab='CPI',xlab='Time')
points(cpi,pch=Month) ; abline(h=mean(cpi))

stl_cpi<-stl(cpi,s.window = "periodic") 
plot(stl_cpi,main="Seasonal Decomposition of CPI by loess") 

visit_cpi_dif=ts.intersect(Visit=diff(diff(log(visit),12)),CPI=cpi)
plot(visit_cpi_dif,yax.flip=T,main="First and Seasonal Differenced Log(visit) and CPI") 

ccf(y=as.vector(visit_cpi_dif[,1]),x=as.vector(visit_cpi_dif[,2]),
    main='#visitors & CPI',ylab='CCF')

prewhiten(x=as.vector(visit_cpi_dif[,2]),y=as.vector(visit_cpi_dif[,1]),
          main='Sample CCF of Prewhitened Differenced Log(Visits) and CPI',ylab='CCF')


############## Model Specification ################
plot(log(visit),ylab='Number of visitors',xlab='Time', 
     main="Logarithm of Monthly Taiwanese visitors to the US")
points(log(visit),pch=Month)
abline(h=mean(log(visit)))  

## look at preintervention data
preint<-mydata$us_source[1:68]
preint<-ts(visit,start=c(1996,1),end=c(2001,8),frequency=12)
plot(preint)

# plot acf 
acf(as.vector(log(preint)),lag.max=36)  # show seasonality
acf(as.vector(diff(log(preint),lag=12)),lag.max=36) ## there is pattern (nonstationary)
acf(as.vector(diff(diff(log(preint)),lag=12)),lag.max=36,ci.type='ma',
    main="ACF of First and Seasonal Difference of logarithm of visitors (pre-intervention)") 
pacf(as.vector(diff(diff(log(preint)),lag=12)),lag.max=36,
     main="PACF of First and Seasonal Difference of logarithm of visitors (pre-intervention)") 

plot(diff(diff(log(preint)),lag=12),xlab='Time') 
abline(h=0)


############## Model fitting ###############
## Model 1 
m1.visit<-arimax(log(visit),order=c(2,1,1),seasonal=list(order=c(2,1,1),period=12),
                 xtransf = data.frame(I911=1*(seq(visit)==69),I911=1*(seq(visit)==69)),
                 transfer = list(c(0,0),c(1,0)),
                 xreg=data.frame(CPI=cpi,SARS=c(rep(0,85),rep(1,5),rep(0,(232-90)))),
                 method='ML')
m1.visit 

## Model 2
m2.visit<-arimax(log(visit),order=c(2,1,1),seasonal=list(order=c(2,1,1),period=12),
                 xtransf = data.frame(I911=1*(seq(visit)==69),I911=1*(seq(visit)==69)),
                 transfer = list(c(0,0),c(1,0)),
                 xreg=data.frame(CPI=cpi,SARS=c(rep(0,85),rep(1,5),rep(0,(232-90)))),
                 fixed = c(0,0,NA,NA,NA,NA,NA,NA,0,NA,NA),method='ML')
m2.visit 


############## Model Diagnostics ###############
detectAO(m2.visit)
detectIO(m2.visit)

m3.visit<-arimax(log(visit),order=c(2,1,1),seasonal=list(order=c(2,1,1),period=12),
                 xtransf = data.frame(I911=1*(seq(visit)==69),I911=1*(seq(visit)==69)),
                 transfer = list(c(0,0),c(1,0)),
                 xreg=data.frame(CPI=cpi,SARS=c(rep(0,85),rep(1,5),rep(0,(232-90)))),
                 io=c(26,62,88,89,122,158),fixed=c(0,0,NA,NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA,0,NA,NA),method='ML')
m3.visit 
detectIO(m3.visit)
detectAO(m3.visit)

mod<-m3.visit 
tsdiag(mod,gof=30,omit.initial=T)
LB.test(mod)
qqnorm(residuals(mod)); qqline(residuals(mod))
hist(residuals(mod),xlab='Standardized Residuals')
shapiro.test(residuals(mod))

periodogram(residuals(mod))
spec(residuals(mod),span=16,main="Smoothed Periodogram")

plot(log(visit),ylab='log(visit)',main="Logs of Number of Visitors and Fitted Values")
points(fitted(mod))


### checking ARCH 
McLeod.Li.test(mod,main="McLeod.Li Test Statistics for m3.visit model")
acf(as.vector(residuals(m3.visit))^2,na.action=na.omit,
    main="Sample ACF of the Squared residuals(m3.visit)")
pacf(as.vector(residuals(m3.visit))^2,na.action=na.omit,
     main="Sample PACF of the Squared residuals(m3.visit)")

acf(abs(as.vector(residuals(m3.visit))),na.action=na.omit,
    main="Sample ACF of the Absolute residuals(m3.visit)")
pacf(abs(as.vector(residuals(m3.visit))),na.action=na.omit,
    main="Sample PACF of the Absolute residuals(m3.visit)")


############ Frequency Analysis #############
sp=periodogram(visit,ylab="Visitors Periodogram")
abline(h=0)
which(sp$spec==max(sp$spec))
1/sp$freq[20] ; sp$freq[20]
1/sp$freq[40] ; sp$freq[40]
1/sp$freq[1]
axis(1,at=c(0.0833,0.1667))

sp2=spec(visit,log='no',sub='',xlab='Frequency',ylab='Smoothed Sample Spectral Density')

k=kernel('daniell',m=2)  # reduce the variability
sp3=spec(visit,kernel=k,log='no',sub='',xlab='Frequency',
         ylab='Smoothed Sample Spectral Density')


########## Computing AIC of different models for preintervention data ############
aic_tab2<-data.frame(matrix(nrow = 18,ncol = 18))
aic_tab2[1,]<-c(NA,NA,rep(1,4),rep(2,4),rep(3,4),rep(4,4))
aic_tab2[2,]<-c(NA,NA,rep(c(1:4),4))
aic_tab2[,1]<-c(NA,NA,rep(1,4),rep(2,4),rep(3,4),rep(4,4))
aic_tab2[,2]<-c(NA,NA,rep(c(1:4),4))
aic_tab2

# i=row, j=column, p=row, q=column
for (i in 7:18){
    for (j in 3:18){
        p<-aic_tab2[1,j] ; q<-aic_tab2[i,1]
        P<-aic_tab2[2,j] ; Q<-aic_tab2[i,2]
        possibleError<-tryCatch(
            model<-arima(log(preint),order=c(p,1,q),seasonal=list(order=c(P,1,Q),period=12)),
            error=function(e) e)
        if (inherits(possibleError,"error")) next
        if (inherits(possibleError,"warning")) next
        aic_tab2[i,j]<-model$aic
    }
}

aic_tab2[is.na(aic_tab2)]<-0
which.min(apply(aic_tab2,MARGIN=2,min)) #column
which.min(apply(aic_tab2,MARGIN=1,min)) #row


