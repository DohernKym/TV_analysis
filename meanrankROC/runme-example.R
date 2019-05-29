#########################################################################################################
# Aasthaa Bansal & Patrick Heagerty 
# University of Washington 
# May 5, 2018
#
# Code used in: 
# 'A tutorial on evaluating the time-varying discrimination accuracy of survival models used in dynamic decision-making', MDM, 2018
######################################################################################################### 
# http://faculty.washington.edu/abansal/software.html #

rm(list=ls())

library(survival)
library(survivalROC)
library(risksetROC)

# Download .R files provided as supplementary materials or download the meanrankROC package from
# github: https://github.com/aasthaa/meanrankROC_package
setwd("F:\\_0.Paper Work\\_4.Statistics\\통계관련자료\\ROC\\Rcode_online_supp")

source("MeanRank.q")
source("NNE-estimate.q")
source("NNE-CrossValidation.q")
source("interpolate.q")

source("dynamicTP.q")
source("NNE-estimate_TPR.q")

source("dynamicIntegrateAUC.R")

#Load in the datasets. Note: The PBC data is freely available in R.
data(pbc)
names(pbc)[names(pbc)=="fudays"] <- "time"
bDat <- pbc[1:312,] #baseline data
bDat$deathEver <- bDat$status
bDat$deathEver[which(bDat$status==1)] <- 0 #censor at transplant
bDat$deathEver[which(bDat$status==2)] <- 1 #death

head(bDat)

# Build dataset with time-dependent covariates
str(pbc)
length(unique(pbc$id))

(pbc2 <- tmerge(pbc, pbc, id=id, death = event(time, status))) # set range
summary(pbc2)
summary(pbcseq)
?tmerge
(pbc2 <- tmerge(pbc2, pbcseq, id=id, ascites = tdc(day, ascites), hepato = tdc(day, hepato),
   spiders = tdc(day, spiders), edema = tdc(day, edema), chol = tdc(day, chol), 
   bili = tdc(day, bili), albumin = tdc(day, albumin),
   protime = tdc(day, protime), alk.phos = tdc(day, alk.phos),
   ast = tdc(day, ast), platelet = tdc(day, platelet), stage = tdc(day, stage)))
colnames(pbc2)
length(unique(pbc2$id))
attr(pbc2, "tcount")
dim(pbc2)
pbc2 <- subset(pbc2, id>=1 & id<=312)
dim(pbc2)
length(unique(pbc2$id))
dim(bDat)

# According to documentation, some baseline values for protime and age 
  # in pbc were found to be incorrect. Correct values in pbcseq
bDat[1:5,]
subset(pbc2, tstart==0)[1:5,]
for(i in 1:312){
   if(pbc2$protime[which(pbc2$id==i & pbc2$tstart==0)] != bDat$protime[i])
      pbc2$protime[which(pbc2$id==i & pbc2$tstart==0)] <- bDat$protime[i]
      
   if(pbc2$age[which(pbc2$id==i & pbc2$tstart==0)] != bDat$age[i])
      pbc2$age[which(pbc2$id==i & pbc2$tstart==0)] <- bDat$age[i]
}

pbc2$deathEver <- pbc2$status
pbc2$deathEver[which(pbc2$status==1)] <- 0
pbc2$deathEver[which(pbc2$status==2)] <- 1

pbc2$death[which(pbc2$death==1)] <- 0
pbc2$death[which(pbc2$death==2)] <- 1

# Use 10-fold CV to get baseline scores
set.seed(49)

samples <- floor(runif(nrow(bDat), 1,11))
sampSizes <- sapply(seq(1:10), function(s){length(which(samples==s))} )
sampSizes # Check that no subsets with 0 subjects

while(min(sampSizes)==0) {
   samples <- floor(runif(nTrain, 1,11))
   sampSizes <- sapply(seq(1:10), function(s){length(which(samples==s))} )
}

### 10-fold cross-validation to get predicted baseline scores
score4Baseline_cv <- score5Baseline_cv <- rep(NA,nrow(bDat))

for(s in 1:10) {
   bDat_train <- bDat[-which(samples==s),]
   bDat_test <- bDat[which(samples==s),]
          
   mod <- coxph(Surv(time=time, event= deathEver) ~ log(bili) + log(protime) + edema + albumin + age,  
      data=bDat_train )         
   riskVals <- predict(mod, type="risk", newdata= bDat_test)
   score5Baseline_cv[which(samples==s)] <- riskVals
       
   mod <- coxph(Surv(time=time, event= deathEver) ~ log(protime) + edema + albumin + age, data=bDat_train )         
   riskVals <- predict(mod, type="risk", newdata= bDat_test)
   score4Baseline_cv[which(samples==s)] <- riskVals
}
bDat$score4baseline <- score4Baseline_cv
bDat$score5baseline <- score5Baseline_cv

# Fit model on all baseline data, use for prediction of time-varying score
coxMod5baseline <- coxph(Surv(time=time, event= deathEver) ~ log(bili) + log(protime) + edema + albumin + age, data= bDat)
pbc2$score5tv <- predict(coxMod5baseline, type="risk", newdata= pbc2)

coxMod4baseline <- coxph(Surv(time=time, event= deathEver) ~ log(protime) + edema + albumin + age, data= bDat)
pbc2$score4tv <- predict(coxMod4baseline, type="risk", newdata= pbc2)

##### Table 3
head(pbcseq)
## landmarkTimes
units <- 365.25
landmarkTimes <- c(1, 4, 6)*units

## A. AUC_I/D
tableAUC_ID <- matrix(nrow=2, ncol=length(landmarkTimes))
tableAUC_TV_ID <- matrix(nrow=2, ncol=length(landmarkTimes))

# Baseline risk scores
scores <- c("score4baseline","score5baseline")
for(i in 1:length(scores)) {
   currVar <- eval(parse(text=paste("bDat$", scores[i],sep="")))  

   mmm <- MeanRank( survival.time= bDat$time, survival.status= bDat$deathEver, marker= currVar )

   bandwidths <- 0.05 + c(1:80)/200
   IMSEs <- vector(length=length(bandwidths))
   for(j in 1:length(bandwidths)) {
      nnnC <- nne.CrossValidate( x= mmm$time, y= mmm$mean.rank, lambda=bandwidths[j] )  #Cross-validated bandwidth
      IMSEs[j] <- nnnC$IMSE
   }

   currLambdaOS <- mean(bandwidths[which(IMSEs==min(IMSEs, na.rm=T))])
   
   nnn <- nne( x= mmm$time, y= mmm$mean.rank, lambda=currLambdaOS, nControls=mmm$nControls )  #Fixed bandwidth
   tableAUC_ID[i,] <- sapply(landmarkTimes, function(x){ interpolate( x = nnn$x, y=nnn$nne, target=x ) } )
}
rownames(tableAUC_ID) <- scores
colnames(tableAUC_ID) <- landmarkTimes/units
round(tableAUC_ID, 3)

# Updated (time-varying) risk scores
scores <- c("score4tv", "score5tv")
for(i in 1:length(scores)) {
   currVar <- eval(parse(text=paste("pbc2$", scores[i],sep="")))  

   mmm <- MeanRank(survival.time= pbc2$tstop, survival.status= pbc2$death, marker= currVar, start= pbc2$tstart)

   bandwidths <- 0.05 + c(1:80)/200
   IMSEs <- vector(length=length(bandwidths))
   for(j in 1:length(bandwidths)) {
      nnnC <- nne.CrossValidate( x= mmm$time, y= mmm$mean.rank, lambda=bandwidths[j] )  # Cross-validated bandwidth
      IMSEs[j] <- nnnC$IMSE
   }

   currLambdaOS <- mean(bandwidths[which(IMSEs==min(IMSEs, na.rm=T))])
   
   nnn <- nne( x= mmm$time, y= mmm$mean.rank, lambda=currLambdaOS, nControls=mmm$nControls)  # Fixed bandwidth
   tableAUC_TV_ID[i,] <- sapply(landmarkTimes, function(x){interpolate(x = nnn$x, y=nnn$nne, target=x )})
}
rownames(tableAUC_TV_ID) <- scores
colnames(tableAUC_TV_ID) <- landmarkTimes/units
round(tableAUC_TV_ID, 3)

# B. c-index
round(dynamicIntegrateAUC(survival.time=bDat$time, survival.status=bDat$deathEver, marker=bDat$score4baseline, cutoffTime = units*10), 3)
round(dynamicIntegrateAUC(survival.time=bDat$time, survival.status=bDat$deathEver, marker=bDat$score5baseline, cutoffTime = units*10), 3)

round(dynamicIntegrateAUC(survival.time= pbc2$tstop, survival.status= pbc2$death, start=pbc2$tstart, marker=pbc2$score4tv, cutoffTime = units*10), 3)
round(dynamicIntegrateAUC(survival.time= pbc2$tstop, survival.status= pbc2$death, start=pbc2$tstart, marker=pbc2$score5tv, cutoffTime = units*10), 3)

fit1 <- coxph(Surv(time, deathEver) ~ score4baseline, data = bDat)
summary(fit1)
fit2 <- coxph(Surv(time, deathEver) ~ score5baseline, data = bDat)
summary(fit2)
fit3 <- coxph(Surv(tstop, death) ~ score4tv, data = pbc2)
summary(fit3)
fit4 <- coxph(Surv(tstop, death) ~ score5tv, data = pbc2)
summary(fit4)

# C. Sequential C/D AUCs on subsetted data at each timepoint and one year ahead to mimic landmark analysis 
# landmarkTimes ####   
units <- 365.25

landmarkTimes <- c(1, 4, 6)*units
tableAUC_CD <- matrix(nrow=4, ncol=length(landmarkTimes))

timeWindow <- 1

for(j in 1:length(landmarkTimes)) {
   currData <- subset(bDat, time >= (landmarkTimes[j]))
   currDataTV <- subset(pbc2, tstart <= (landmarkTimes[j]) & tstop > (landmarkTimes[j]))

   nobs <- nrow(currData)     
   out1 <- survivalROC( currData$time, currData$deathEver, marker= currData$score4baseline,
              predict.time = (landmarkTimes[j] + timeWindow*units), method = "NNE", span = 0.04 * nobs^(-0.2)  )
   tableAUC_CD[1,j] <- out1$AUC

   out1 <- survivalROC( currData$time, currData$deathEver, marker= currData$score5baseline,
              predict.time = (landmarkTimes[j] + timeWindow*units), method = "NNE", span = 0.04 * nobs^(-0.2)  )
   tableAUC_CD[2,j] <- out1$AUC


   nobs <- nrow(currDataTV)
   out1 <- survivalROC( currDataTV$time, currDataTV$deathEver, marker= currDataTV$score4tv, 
              predict.time = (landmarkTimes[j] + timeWindow*units), method = "NNE", span = 0.04 * nobs^(-0.2)  )
   tableAUC_CD[3,j] <- out1$AUC    
            
   out1 <- survivalROC( currDataTV$time, currDataTV$deathEver, marker= currDataTV$score5tv, 
              predict.time = (landmarkTimes[j] + timeWindow*units), method = "NNE", span = 0.04 * nobs^(-0.2)  )
   tableAUC_CD[4,j] <- out1$AUC
}
rownames(tableAUC_CD) <- c("score4baseline","score5baseline", "score4tv", "score5tv")
colnames(tableAUC_CD) <- landmarkTimes/units
round(tableAUC_CD, 3)

####Bootstrap 95% CIs
nBoot <- 500

## A. Bootstrap CIs - Baseline markers/scores
markers <- c("score4baseline","score5baseline")
set.seed(49)
Cindex_bstrap <- matrix(nrow=nBoot, ncol=length(markers))
bstrapRes <- list()

for(b in 1:nBoot) {
   currData <- bDat[sample(x=seq(1,nrow(bDat)), size=nrow(bDat), replace = T),]   
   kmfit <- survfit(Surv(time, deathEver) ~ 1, data= currData)
   
   currDataLM1 <- currData[which(currData$time >= (landmarkTimes[1])), ]
   currDataLM2 <- currData[which(currData$time >= (landmarkTimes[2])), ]
   currDataLM3 <- currData[which(currData$time >= (landmarkTimes[3])), ]

   aucID_scores <- NULL
   aucCD_scores <- matrix(nrow=length(scores), ncol=length(landmarkTimes))

   for(i in 1:length(markers)) {
      currVar <- eval(parse(text=paste("currData$",markers[i],sep="")))

      ### AUC I/D
      mmm <- MeanRank( survival.time= currData$time, survival.status= currData$deathEver, marker= currVar )
   	  nnn <- nne( x= mmm$time, y= mmm$mean.rank, lambda=0.3 )  #Fixed bandwidth
      aucID_scores <- rbind(aucID_scores, 
            sapply(landmarkTimes, function(t){ interpolate( x = nnn$x, y=nnn$nne, target=t ) } ))
   
      ### C-index
      Cindex_bstrap[b,i] <- dynamicIntegrateAUC(survival.time=currData$time, survival.status= currData$deathEver,
         marker=currVar, cutoffTime = units*10)
   
      ### AUC C/D landmark  
      if(markers[i]=="score4baseline" | markers[i]=="score5baseline") {
      	   
         currDataLM <- currDataLM1
         currVecLM <- eval(parse(text=paste("currDataLM$", markers[i], sep="")))
         out1 <- survivalROC( currDataLM$time, currDataLM$deathEver, marker=currVecLM,
                   predict.time = (landmarkTimes[1] + timeWindow*units), method = "NNE", span = 0.04 * nobs^(-0.2)  )

         currDataLM <- currDataLM2
         currVecLM <- eval(parse(text=paste("currDataLM$", markers[i], sep="")))
         out2 <- survivalROC( currDataLM$time, currDataLM$deathEver, marker=currVecLM,
                   predict.time = (landmarkTimes[2] + timeWindow*units), method = "NNE", span = 0.04 * nobs^(-0.2)  )

         currDataLM <- currDataLM3
         currVecLM <- eval(parse(text=paste("currDataLM$", markers[i], sep="")))
         out3 <- survivalROC( currDataLM$time, currDataLM$deathEver, marker=currVecLM,
                   predict.time = (landmarkTimes[3] + timeWindow*units), method = "NNE", span = 0.04 * nobs^(-0.2)  )

         aucCD_scores[i,] <- c(out1$AUC, out2$AUC, out3$AUC)
      }
   }
   bstrapRes[[b]] <- list(aucID_scores=aucID_scores, aucCD_scores=aucCD_scores)
}

# Get CIs for c-indices
Cindex_CIs <- round(apply(Cindex_bstrap, 2, quantile, probs=c(0.025,0.975)),3)
colnames(Cindex_CIs) <- markers
Cindex_CIs

# Get CIs for AUCs
AUC_ID_CIs <- NULL
AUC_CD_CIs <- NULL

for(t in 1:length(landmarkTimes)) {
   AUC_ID <- NULL
   AUC_CD <- NULL
   for(b in 1:nBoot) {
      AUC_ID <- cbind( AUC_ID, bstrapRes[[b]]$aucID_scores[,t] )
      AUC_CD <- cbind( AUC_CD, bstrapRes[[b]]$aucCD_scores[,t] ) 
   }
   AUC_ID_CI_raw <- round( apply(AUC_ID, 1, quantile, probs=c(0.025,0.975)), 3 )
   AUC_CD_CI_raw <- round( apply(AUC_CD, 1, quantile, probs=c(0.025,0.975)), 3 )
   
   AUC_ID_CIs <- cbind(AUC_ID_CIs, sapply(seq(2), function(x) paste("(", AUC_ID_CI_raw[1,x], ", ", AUC_ID_CI_raw[2,x], ")", sep="" ) ) )
   AUC_CD_CIs <- cbind(AUC_CD_CIs, sapply(seq(2), function(x) paste("(", AUC_CD_CI_raw[1,x], ", ", AUC_CD_CI_raw[2,x], ")", sep="" ) ) )
}
rownames(AUC_CD_CIs) <- rownames(AUC_ID_CIs) <- c("4-cov model", "5-cov model")
colnames(AUC_CD_CIs) <- colnames(AUC_ID_CIs) <- landmarkTimes/units
AUC_ID_CIs
AUC_CD_CIs

## B. Bootstrap CIs - Time-varying scores
markers <- c("score4tv","score5tv")

set.seed(49)
Cindex_bstrapTV <- matrix(nrow=nBoot, ncol=length(markers))
bstrapResTV <- list()

for(b in 1:nBoot) {
   #sample individuals
   subjs <- unique(pbc2$id)
   currSubjs <- sample(x=subjs, size=length(subjs), replace = T)
   currData <- NULL
   for(j in 1:length(currSubjs)) 
      currData <- rbind(currData, pbc2[which(pbc2$id==currSubjs[j]),]) 
   
   kmfit <- survfit(Surv(time=tstart, time2=tstop, event=death) ~ 1, data=currData)
      
   currDataLM1 <- currData[which(currData$tstart <= (landmarkTimes[1]) & currData$tstop > (landmarkTimes[1])), ]
   currDataLM2 <- currData[which(currData$tstart <= (landmarkTimes[2]) & currData$tstop > (landmarkTimes[2])), ]
   currDataLM3 <- currData[which(currData$tstart <= (landmarkTimes[3]) & currData$tstop > (landmarkTimes[3])), ]


   aucID_scores <- NULL
   aucCD_scores <- matrix(nrow=length(scores), ncol=length(landmarkTimes))

   for(i in 1:length(markers)) {
      currVar <- eval(parse(text=paste("currData$", markers[i],sep="")))

      ### AUC I/D
      mmm <- MeanRank( survival.time= currData$tstop, survival.status= currData$death, start=currData$tstart, marker= currVar )
   	  nnn <- nne( x= mmm$time, y= mmm$mean.rank, lambda=0.3, nControls=mmm$nControls )  #Fixed bandwidth
      aucID_scores <- rbind(aucID_scores, 
            sapply(landmarkTimes, function(t){ interpolate( x = nnn$x, y=nnn$nne, target=t ) } ))

      ### C-index
      Cindex_bstrapTV[b,i] <- dynamicIntegrateAUC(survival.time=currData$tstop, survival.status=currData$death,
         start=currData$tstart, marker= currVar, cutoffTime = units*10)   
   
      ### AUC C/D landmark (for the scores only)       	   
         currDataLM <- currDataLM1
         currVecLM <- eval(parse(text=paste("currDataLM$", markers[i], sep="")))
         out1 <- survivalROC( currDataLM$time, currDataLM$deathEver, marker=currVecLM,
                   predict.time = (landmarkTimes[1] + timeWindow*units), method = "NNE", span = 0.04 * nobs^(-0.2)  )

         currDataLM <- currDataLM2
         currVecLM <- eval(parse(text=paste("currDataLM$", markers[i], sep="")))
         out2 <- survivalROC( currDataLM$time, currDataLM$deathEver, marker=currVecLM,
                   predict.time = (landmarkTimes[2] + timeWindow*units), method = "NNE", span = 0.04 * nobs^(-0.2)  )

         currDataLM <- currDataLM3
         currVecLM <- eval(parse(text=paste("currDataLM$", markers[i], sep="")))
         out3 <- survivalROC( currDataLM$time, currDataLM$deathEver, marker=currVecLM,
                   predict.time = (landmarkTimes[3] + timeWindow*units), method = "NNE", span = 0.04 * nobs^(-0.2)  )

         aucCD_scores[i,] <- c(out1$AUC, out2$AUC, out3$AUC)

   }
   bstrapResTV[[b]] <- list(aucID_scores, aucCD_scores)
}

# Get CIs for c-indices
Cindex_CIs <- round(apply(Cindex_bstrapTV, 2, quantile, probs=c(0.025,0.975)), 3)
colnames(Cindex_CIs) <- markers
Cindex_CIs

# Get CIs for AUCs
AUC_CD_CIs <- NULL
AUC_ID_CIs <- NULL

for(t in 1:length(landmarkTimes)) {
   AUC_CD <- NULL
   AUC_ID <- NULL
   for(b in 1:nBoot) {
      AUC_ID <- cbind( AUC_ID, bstrapResTV[[b]][[1]][,t] )
      AUC_CD <- cbind( AUC_CD, bstrapResTV[[b]][[2]][,t] ) 
   }
   AUC_CD_CI_raw <- round( apply(AUC_CD, 1, quantile, probs=c(0.025,0.975)), 3 )
   AUC_ID_CI_raw <- round( apply(AUC_ID, 1, quantile, probs=c(0.025,0.975)), 3 )
   
   AUC_CD_CIs <- cbind(AUC_CD_CIs, sapply(seq(2), function(x) paste("(", AUC_CD_CI_raw[1,x], ", ", AUC_CD_CI_raw[2,x], ")", sep="" ) ) )
   AUC_ID_CIs <- cbind(AUC_ID_CIs, sapply(seq(2), function(x) paste("(", AUC_ID_CI_raw[1,x], ", ", AUC_ID_CI_raw[2,x], ")", sep="" ) ) )
}
rownames(AUC_CD_CIs) <- rownames(AUC_ID_CIs) <- c("4-cov model", "5-cov model")
colnames(AUC_CD_CIs) <- colnames(AUC_ID_CIs) <- landmarkTimes/units
AUC_ID_CIs
AUC_CD_CIs

# C-index difference
getCindexBstrapCI <- function(nBoot, inData, markerVarName1, markerVarName2, timeVarName, eventVarName, cutoffTime) {
   set.seed(49)
   resultStar <- vector(length=nBoot)
   for(i in 1:nBoot) {
   	  datStar <- inData[sample(seq(nrow(inData)), nrow(inData), replace=T), ]
   	  
   	  markerVar1 <- eval(parse(text=paste("datStar$", markerVarName1, sep="")))
   	  markerVar2 <- eval(parse(text=paste("datStar$", markerVarName2, sep="")))
   	  timeVar <- eval(parse(text=paste("datStar$", timeVarName, sep="")))
   	  eventVar <- eval(parse(text=paste("datStar$", eventVarName, sep="")))

      kmfit <- survfit(Surv(timeVar, eventVar) ~ 1)


      ### Marker 1
      mmm <- MeanRank( survival.time= timeVar, survival.status= eventVar, marker= markerVar1 )

      #Get overlap between survival function and mmm
      meanRanks <-  mmm$mean.rank[which(mmm$time <= cutoffTime)]
      survTimes <- mmm$time[mmm$time <= cutoffTime]
      timeMatch <- match(survTimes, kmfit$time)

      S_t <- kmfit$surv[timeMatch]

      #Calculate weights for c-index
      f_t <- c( 0, (S_t[-length(S_t)] - S_t[-1]) )
      S_tao <- S_t[length(S_t)]
      weights <- (2*f_t*S_t)/(1-S_tao^2)
      
      Cindex1 <- sum(meanRanks * weights) #C-index


      ### Marker 2
      mmm <- MeanRank( survival.time= timeVar, survival.status= eventVar, marker= markerVar2 )

      #Get overlap between survival function and mmm
      meanRanks <-  mmm$mean.rank[which(mmm$time <= cutoffTime)]
      survTimes <- mmm$time[mmm$time <= cutoffTime]
      timeMatch <- match(survTimes, kmfit$time)

      S_t <- kmfit$surv[timeMatch]

      #Calculate weights for c-index
      f_t <- c( 0, (S_t[-length(S_t)] - S_t[-1]) )
      S_tao <- S_t[length(S_t)]
      weights <- (2*f_t*S_t)/(1-S_tao^2)
      
      Cindex2 <- sum(meanRanks * weights) #C-index
      
      resultStar[i] <- Cindex1 - Cindex2
      
   }
   return( quantile(resultStar, probs=c(0.025, 0.975)) )
}

getCindexBstrapCI(nBoot=500, inData=bDat, markerVarName1="score5baseline", markerVarName2="score4baseline",
   timeVarName="time", eventVarName="deathEver", cutoffTime=units*10)

### FIGURES

# Figure 2
# pdf("Figure2_meanRank_pbc_baseline.pdf", height=4, width=9)
# par(mfrow=c(1,2))
# par(ps=10)

# AUC I/D
mmmBaseline <- MeanRank( survival.time= bDat$time, survival.status= bDat$deathEver, marker= bDat$score5baseline )
print(length(mmmBaseline$time))
nnnBaseline <- nne( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.2, nControls= mmmBaseline$nControls )  #Fixed bandwidth
plot( mmmBaseline$time, mmmBaseline$mean.rank, xlab="Time (years)", ylab=expression(AUC^"I/D"*"(t)"), ylim=c(0.4,1), col="blue", pch=21, cex=.8, axes=F, xlim=c(0,11)*units)
axis(1, at=seq(0,10,by=2)*units, labels=seq(0,10,by=2))
axis(2)
box()
abline(h=0.5, lty=2)
lines( nnnBaseline$x, nnnBaseline$nne, col="blue", lwd=2 )

mmmBaseline <- MeanRank( survival.time= bDat$time, survival.status= bDat$deathEver, marker= bDat$score4baseline )
print(length(mmmBaseline$time))
nnnBaseline <- nne( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.2, nControls= mmmBaseline$nControls )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="orange", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="orange", lwd=2 )

# ROC I/D (TPR)
fpr <- 0.1

mmmBaseline <- dynamicTP( p=fpr, survival.time= bDat$time, survival.status= bDat$deathEver, marker= bDat$score5baseline )
print(length(mmmBaseline$time))
nnnBaseline <- nne_TPR( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda= 0.3, nControls= mmmBaseline$nControls,
   nCases= mmmBaseline$nCases, p=fpr, survival.time= bDat$time, survival.status= bDat$deathEver, marker= bDat$score5baseline )  #Fixed bandwidth
plot( mmmBaseline$time, mmmBaseline$mean.rank, xlab="Time (years)", ylab=expression(ROC[t]^"I/D"*"(FPF=10%)"), ylim=c(0,1), col="blue", pch=21, cex=.8, axes=F, xlim=c(0,11)*units)
axis(1, at=seq(0,10,by=2)*units, labels=seq(0,10,by=2))
axis(2)
box()
abline(h=0.5, lty=2)
lines( nnnBaseline$x, nnnBaseline$nne, col="blue", lwd=2 )

mmmBaseline <- dynamicTP( p=fpr, survival.time= bDat$time, survival.status= bDat$deathEver, marker= bDat$score4baseline )
print(length(mmmBaseline$time))
nnnBaseline <- nne_TPR( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.3, nControls= mmmBaseline$nControls,
   nCases= mmmBaseline$nCases, p=fpr, survival.time= bDat$time, survival.status= bDat$deathEver, marker= bDat$score4baseline )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="orange", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="orange", lwd=2 )

legend(x=1.5*units, y=0.9, legend=c("4 covariates", "5 covariates"), col=c("orange","blue"), lty=1, lwd=2, horiz=T)

dev.off()

# Figure 3
mmm <- MeanRank( survival.time= pbc2$tstop, survival.status= pbc2$death, marker= currVar, start= pbc2$tstart )

pdf("Figure3_meanRank_pbc_TV.pdf", height=4, width=9)
par(mfrow=c(1,2))
par(ps=10)

# AUC I/D
mmmTV <- MeanRank( survival.time= pbc2$tstop, survival.status= pbc2$death, marker= pbc2$score5tv, start= pbc2$tstart )
print(length(mmmTV$time))
nnn <- nne( x= mmmTV$time, y= mmmTV$mean.rank, lambda=0.2, nControls=mmmTV$nControls )  #Fixed bandwidth
plot( mmmTV$time, mmmTV$mean.rank, xlab="Time (years)", ylab=expression(AUC^"I/D"*"(t)"), ylim=c(0.4,1), col="blue", pch=21, cex=.8, axes=F, xlim=c(0,11)*units)
axis(1, at=seq(0,10,by=2)*units, labels=seq(0,10,by=2))
axis(2)
box()
abline(h=0.5, lty=2)
lines( nnn$x, nnn$nne, col="blue", lwd=2, lty=1 )

mmmTV <- MeanRank( survival.time= pbc2$tstop, survival.status= pbc2$death, marker= pbc2$score4tv, start= pbc2$tstart )
print(length(mmmTV$time))
nnn <- nne( x= mmmTV$time, y= mmmTV$mean.rank, lambda=0.2, nControls=mmmTV$nControls )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="orange", pch=21, cex=.8)
lines( nnn$x, nnn$nne, col="orange", lwd=2, lty=1 )

# ROC I/D (TPR)
mmmTV <- dynamicTP( p=fpr, survival.time =pbc2$tstop, survival.status= pbc2$death, marker= pbc2$score5tv, start= pbc2$tstart )
print(length(mmmTV$time))
nnn <- nne_TPR( x= mmmTV$time, y= mmmTV$mean.rank, lambda= 0.3, nControls= mmmTV$nControls,
   nCases= mmmTV$nCases, p=fpr, survival.time =pbc2$tstop, survival.status= pbc2$death, marker= pbc2$score5tv, start= pbc2$tstart )  #Fixed bandwidth
plot( mmmTV$time, mmmTV$mean.rank, xlab="Time (years)", ylab=expression(ROC[t]^"I/D"*"(FPF=10%)"), ylim=c(0,1), col="blue", pch=21, cex=.8, axes=F, xlim=c(0,11)*units)
axis(1, at=seq(0,10,by=2)*units, labels=seq(0,10,by=2))
axis(2)
box()
abline(h=0.5, lty=2)
lines( nnn$x, nnn$nne, col="blue", lwd=2 )

mmmTV <- dynamicTP( p=fpr, survival.time =pbc2$tstop, survival.status= pbc2$death, marker= pbc2$score4tv, start= pbc2$tstart )
nnn <- nne_TPR( x= mmmTV$time, y= mmmTV$mean.rank, lambda= 0.3, nControls= mmmTV$nControls,
   nCases= mmmTV$nCases, p=fpr, survival.time =pbc2$tstop, survival.status= pbc2$death, marker= pbc2$score4tv, start= pbc2$tstart )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="orange", pch=21, cex=.8 )
lines( nnn$x, nnn$nne, col="orange", lwd=2 )

legend(x=2.5*units, y=0.3, legend=c("4 covariates", "5 covariates"), col=c("orange","blue"), lty=1, lwd=2, horiz=T)

dev.off()

## Figure 4: With CIs for baseline and updated risk scores from 5-covariate model using the incident/dynamic approach
pdf("Figure4_meanRank_pbc_baseline_TV_CI.pdf", height=4, width=9)
par(mfrow=c(1,2))
par(ps=10)

# ROC I/D
mmmBaseline <- MeanRank( survival.time= bDat$time, survival.status= bDat$deathEver, marker= bDat$score5baseline )
nnnBaseline <- nne( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.2, nControls= mmmBaseline$nControls )  #Fixed bandwidth
plot( mmmBaseline$time, mmmBaseline$mean.rank, xlab="Time (years)", ylab=expression(AUC^"I/D"*"(t)"), ylim=c(0.4,1), col="lightblue", pch=21, cex=.8, axes=F, xlim=c(0,11)*units)
axis(1, at=seq(0,10,by=2)*units, labels=seq(0,10,by=2))
axis(2)
box()
abline(h=0.5, lty=2)
lines( nnnBaseline$x, nnnBaseline$nne, col="blue", lwd=2 )
lines( nnnBaseline$x, nnnBaseline$nne + 1.96*sqrt(nnnBaseline$var), col="blue", lty=2 )
lines( nnnBaseline$x, nnnBaseline$nne - 1.96*sqrt(nnnBaseline$var), col="blue", lty=2 )

mmmTV <- MeanRank( survival.time= pbc2$tstop, survival.status= pbc2$death, marker= pbc2$score5tv, start= pbc2$tstart )
nnn <- nne( x= mmmTV$time, y= mmmTV$mean.rank, lambda=0.2, nControls=mmmTV$nControls )  #Fixed bandwidth
plot( mmmTV$time, mmmTV$mean.rank, xlab="Time (years)", ylab=expression(AUC^"I/D"*"(t)"), ylim=c(0.4,1), col="lightblue", pch=21, cex=.8, axes=F, xlim=c(0,11)*units)
axis(1, at=seq(0,10,by=2)*units, labels=seq(0,10,by=2))
axis(2)
box()
abline(h=0.5, lty=2)
lines( nnn$x, nnn$nne, col="blue", lwd=2, lty=1 )
lines( nnn$x, nnn$nne + 1.96*sqrt(nnn$var), col="blue", lty=2 )
lines( nnn$x, nnn$nne - 1.96*sqrt(nnn$var), col="blue", lty=2 )

