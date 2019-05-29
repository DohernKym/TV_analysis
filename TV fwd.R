setwd("F:\\_0.Paper Work\\_0.Ongoing\\.Git Work\\TV analysis")

rm(list=ls())

library(tidyverse)
library(data.table)

library(survival)
library(survivalROC)

library(nlme) 
library(lme4) 

# Download .R files provided as supplementary materials or download the meanrankROC package from
# github: https://github.com/aasthaa/meanrankROC_package

## loading "meanrankROC" R package, reference : http://faculty.washington.edu/abansal/software.html
source("meanrankROC/MeanRank.q")
source("meanrankROC/NNE-estimate.q")
source("meanrankROC/NNE-CrossValidation.q")
source("meanrankROC/interpolate.q")
source("meanrankROC/dynamicTP.q")
source("meanrankROC/NNE-estimate_TPR.q")
source("meanrankROC/dynamicIntegrateAUC.R")

# read longitudinal data
df <- read.csv("F:\\_0.Paper Work\\_0.Ongoing\\_Longitudinal\\sepsis.csv", header = T, sep = ",", na.strings = c("", " ", "NA", "#N/A"), stringsAsFactors = T) %>%
  dplyr::mutate(sex=factor(sex, levels = c("male", "female"), labels = c("Male", "Female")),
                inh=factor(inh, levels = c(0, 1), labels = c("No Inhaltion", "Inhalation")),
                sepsis=factor(sepsis, levels = c(0, 1), labels = c("No Sepsis", "Sepsis")),
                mortality=factor(death, levels = c(0, 1), labels = c("Survivors", "Non-survivors"))) %>%
  group_by(id) %>% 
  dplyr::mutate(days = max(hd)+1)%>% ungroup %>%
  dplyr::select(-plt, -wbc, -pt, -tb, -cr, -crp, -glc, -lct, -pct, -ld, -pf, everything()) %>%
  filter(rowSums(is.na(.))!=11) # 검사값이 없는 열 제외.
length(unique(df$id))

# tranform to wide for each subject
dfw <- df %>% setDT(.) %>%
  dcast(., id+age+sex+tbsa+inh+apachIV+han.score+ABSI+days+mortality~hd,
        value.var = c("plt")) %>% dplyr::select(1:10) %>% 
  mutate(mortality = as.numeric(mortality=="Non-survivors"))

# the creation of start,stop datasets which have multiple intervals for each subject, 
# along with the covariate values that apply over that interval. 
(dfw1 <- tmerge(dfw, dfw, id=id, deathEvent = event(days, mortality))) # set range
all.equal(dfw1$mortality, dfw1$deathEvent)
(dfl <- tmerge(dfw1, df, id=id, wbc = tdc(hd, wbc), 
                plt = tdc(hd, plt), glc = tdc(hd, glc), cr = tdc(hd, cr), 
                crp = tdc(hd, crp), lct = tdc(hd, lct), pct = tdc(hd, pct),
                tb = tdc(hd, tb), pt = tdc(hd, pt), pf = tdc(hd, pf),
                ld = tdc(hd, ld)))
dfl1 <- dfl %>% mutate(plt = -1*plt, pf = -1*pf) 
  # plt와 pf는 증가할 수록 사망할 가능성이 낮아 지는 값으로 원래 값으로 하먄 AUC값이 0.5보다 작게 나옴.
  # 작게 나온 AUC를 1-AUC로 계산하는 것이 맞는지, -1을 곱해서 숫자가 커질수록 사망 이벤트가 높아지는 걸로 바꾸어서 
  # 보는 것이 맞는지 잘 모르겠습니다.

attr(dfl, "tcount")

b_df <- dfl1 %>% filter(tstart==0) # baseline data

nrow(dfl1[dfl1$tstart==0,])
length(unique(b_df$id))

str(dfl)
summary(dfl)

### 전체적으로 생존에 영향을 미치는 여부 보기
  # formula가 맞는지 모르겠습니다.
glmm.fwd <- glmer(death ~ age+tbsa+inh
                  +tb*hd+wbc*hd+pt*hd+tb*hd+cr*hd+crp*hd
                  +glc*hd+lct*hd+ld*hd+pf*hd+plt*hd
                  + (1+hd|id), data=df, family = binomial)

summary(glmm.fwd)
print(glmm.fwd, correlation=T)
  ## 결과 해석은 어떻게 해야 하는지?

glmm.fwd.no <- glmer(death ~ age+tbsa+inh
                  +tb+wbc+pt+tb+cr+crp+glc+lct+ld+pf+plt+hd
                  + (1+hd|id), data=df, family = binomial)

summary(glmm.fwd.no)
print(glmm.fwd.no, correlation=T)


### 각각의 검사가 생존과 사망에 따라 시간에 따라 차이가 있는지 분석 
nlme_tb <- lme(tb ~ hd*death,  random = ~hd|id, data = df, na.action = na.omit)
nlme_lct <- lme(lct ~ hd*death,  random = ~hd|id, data = df, na.action = na.omit)
nlme_pt <- lme(pt ~ hd*death,  random = ~hd|id, data = df, na.action = na.omit)
nlme_tb <- lme(tb ~ hd*death,  random = ~hd|id, data = df, na.action = na.omit)
nlme_wbc <- lme(wbc ~ hd*death,  random = ~hd|id, data = df, na.action = na.omit)
nlme_crp <- lme(crp ~ hd*death,  random = ~hd|id, data = df, na.action = na.omit)
nlme_cr <- lme(cr ~ hd*death,  random = ~hd|id, data = df, na.action = na.omit)
nlme_ld <- lme(ld ~ hd*death,  random = ~hd|id, data = df, na.action = na.omit)
nlme_glc <- lme(glc ~ hd*death,  random = ~hd|id, data = df, na.action = na.omit)
nlme_pf <- lme(pf ~ hd*death,  random = ~hd|id, data = df, na.action = na.omit)
nlme_plt <- lme(plt ~ hd*death,  random = ~hd|id, data = df, na.action = na.omit)

## 결과해석을 어떻게 해야 하는지? 두군에서 시간에따라 의미 있는 차이가 있다고 이야기 할 수 있는지?
summary(nlme_tb)
summary(nlme_lct)
summary(nlme_pt)
summary(nlme_tb)
summary(nlme_wbc)
summary(nlme_crp)
summary(nlme_cr)
summary(nlme_ld)
summary(nlme_glc)
summary(nlme_pf)
summary(nlme_plt)

## 평균 분포 plot

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, 
##   and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  datac$nci <- ciMult*sqrt(datac$sd/datac$N + datac$sd^2/(2*(datac$N-1)))
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  return(datac)
}

#### plt ####
plot_plt <- summarySE(df, measurevar="plt", groupvars=c("hd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_plt, mapping = aes(x = as.numeric(as.character(hd)), y = plt, 
                         group = mortality,
                         color = mortality)) +
  geom_errorbar(aes(ymin=plt-nci, ymax=plt+nci), width=0.7) +
  geom_line() +
  geom_point() +
  xlab("Day before Outcome") +
  ylab("plt") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  scale_y_continuous(limits = c(0, 550), breaks= c(0:550*50)) + # Set tick every 4
  scale_x_continuous(limits = c(0, 60), breaks = c(0:60*10)) +
  theme_bw() +
  theme(legend.position=c(0.75, 0.88)) + 
  ggtitle("Trajectory of WBC") +
  theme(plot.title = element_text(hjust = 0.5))

#### pt ####
plot_pt <- summarySE(df, measurevar="pt", groupvars=c("hd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_pt, mapping = aes(x = as.numeric(as.character(hd)), y = pt, 
                         group = mortality,
                         color = mortality)) +
  geom_errorbar(aes(ymin=pt-nci, ymax=pt+nci), width=0.7) +
  geom_line() +
  geom_point() +
  xlab("Day before Outcome") +
  ylab("pt") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  ggtitle("Tracjetory of pt before Last Outcome") +
  scale_y_continuous(limits = c(10, 33), breaks= c(0:20*5)) + # Set tick every 4
  scale_x_continuous(limits = c(0, 60), breaks = c(0:60*10)) +
  theme(legend.position=c(0.75, 0.88)) + 
  ggtitle("Trajectory of Prothrombin time") +
  theme(plot.title = element_text(hjust = 0.5))

#### tb ####
plot_tb <- summarySE(df, measurevar="tb", groupvars=c("hd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_tb, mapping = aes(x = as.numeric(as.character(hd)), y = tb, 
                         group = mortality,
                         color = mortality)) +
  geom_errorbar(aes(ymin=tb-nci, ymax=tb+nci), width=0.7) +
  geom_line() +
  geom_point() +
  xlab("Day before Outcome") +
  ylab("tb") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  ggtitle("Tracjetory of tb before Last Outcome") +
  scale_y_continuous(limits = c(0, 8), breaks= c(0:8*1)) + # Set tick every 4
  scale_x_continuous(limits = c(0, 60), breaks = c(0:60*10)) +
  theme(legend.position=c(0.75, 0.88)) + 
  ggtitle("Trajectory of Bilirubin") +
  theme(plot.title = element_text(hjust = 0.5))

#### lct ####
plot_lct <- summarySE(df, measurevar="lct", groupvars=c("hd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_lct, mapping = aes(x = as.numeric(as.character(hd)), y = lct, 
                         group = mortality,
                         color = mortality)) +
  geom_errorbar(aes(ymin=lct-nci, ymax=lct+nci), width=0.7) +
  geom_line() +
  geom_point() +
  xlab("Day before Outcome") +
  ylab("lct") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  ggtitle("Tracjetory of lct before Last Outcome") +
  scale_y_continuous(limits = c(0, 7), breaks= c(0:8*2)) + # Set tick every 4
  scale_x_continuous(limits = c(0, 60), breaks = c(0:60*10)) +
  theme(legend.position=c(0.75, 0.88)) + 
  ggtitle("Trajectory of WBC") +
  theme(plot.title = element_text(hjust = 0.5))

#### cr ####
plot_cr <- summarySE(df, measurevar="cr", groupvars=c("hd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_cr, mapping = aes(x = as.numeric(as.character(hd)), y = cr, 
                         group = mortality,
                         color = mortality)) +
  geom_errorbar(aes(ymin=cr-nci, ymax=cr+nci), width=0.7) +
  geom_line() +
  geom_point() +
  xlab("Day before Outcome") +
  ylab("cr") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  ggtitle("Tracjetory of cr before Last Outcome") +
  scale_y_continuous(limits = c(0, 2), breaks= c(0:2*1)) + # Set tick every 4
  scale_x_continuous(limits = c(0, 60), breaks = c(0:60*10)) +
  theme(legend.position=c(0.77, 0.88)) + 
  ggtitle("Trajectory of WBC") +
  theme(plot.title = element_text(hjust = 0.5))

#### glc ####
plot_glc <- summarySE(df, measurevar="glc", groupvars=c("hd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_glc, mapping = aes(x = as.numeric(as.character(hd)), y = glc, 
                         group = mortality,
                         color = mortality)) +
  geom_errorbar(aes(ymin=glc-nci, ymax=glc+nci), width=0.7) +
  geom_line() +
  geom_point() +
  xlab("Day before Outcome") +
  ylab("glc") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  ggtitle("Tracjetory of glc before Last Outcome") +
  scale_y_continuous(limits = c(100, 300), breaks= c(0:400*50)) + # Set tick every 4
  scale_x_continuous(limits = c(0, 60), breaks = c(0:60*10)) +
  theme(legend.position=c(0.75, 0.88)) + 
  ggtitle("Trajectory of WBC") +
  theme(plot.title = element_text(hjust = 0.5))


#### crp ####
plot_crp <- summarySE(df, measurevar="crp", groupvars=c("hd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_crp, mapping = aes(x = as.numeric(as.character(hd)), y = crp, 
                         group = mortality,
                         color = mortality)) +
  geom_errorbar(aes(ymin=crp-nci, ymax=crp+nci), width=0.7) +
  geom_line() +
  geom_point() +
  xlab("Day before Outcome") +
  ylab("crp") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  ggtitle("Tracjetory of crp before Last Outcome") +
  scale_y_continuous(limits = c(0, 250), breaks= c(0:400*50)) + # Set tick every 4
  scale_x_continuous(limits = c(0, 60), breaks = c(0:60*10)) +
  theme(legend.position=c(0.77, 0.88)) + 
  ggtitle("Trajectory of WBC") +
  theme(plot.title = element_text(hjust = 0.5))

#### wbc ####
plot_wbc <- summarySE(df, measurevar="wbc", groupvars=c("hd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_wbc, mapping = aes(x = as.numeric(as.character(hd)), y = wbc, 
                         group = mortality,
                         color = mortality)) +
  geom_errorbar(aes(ymin=wbc-nci, ymax=wbc+nci), width=0.7) +
  geom_line() +
  geom_point() +
  xlab("Day before Outcome") +
  ylab("wbc") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  ggtitle("Trajectory of Biomarker") +
  scale_y_continuous(limits = c(5, 32), breaks= c(0:40*5)) + # Set tick every 4
  scale_x_continuous(limits = c(0, 60), breaks = c(0:60*10)) +
  theme(legend.position=c(0.78, 0.9)) + 
  ggtitle("Trajectory of WBC") +
  theme(plot.title = element_text(hjust = 0.5))

#### pf ####
plot_pf <- summarySE(df, measurevar="pf", groupvars=c("hd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_pf, mapping = aes(x = as.numeric(as.character(hd)), y = pf, 
                               group = mortality,
                               color = mortality)) +
  geom_errorbar(aes(ymin=pf-nci, ymax=pf+nci), width=0.7) +
  geom_line() +
  geom_point() +
  xlab("Day before Outcome") +
  ylab("PF ratio") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  ggtitle("Tracjetory of crp before Last Outcome") +
  scale_y_continuous(limits = c(150, 400), breaks= c(0:400*50)) + # Set tick every 4
  scale_x_continuous(limits = c(0, 61), breaks = c(0:60*10)) +
  theme(legend.position=c(0.75, 0.88)) + 
  ggtitle("Trajectory of PF ratio") +
  theme(plot.title = element_text(hjust = 0.5))

#### ld ####
plot_ld <- summarySE(df, measurevar="ld", groupvars=c("hd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_ld, mapping = aes(x = as.numeric(as.character(hd)), y = ld, 
                               group = mortality,
                               color = mortality)) +
  geom_errorbar(aes(ymin=ld-nci, ymax=ld+nci), width=0.7) +
  geom_line() +
  geom_point() +
  xlab("Day before Outcome") +
  ylab("ld") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  ggtitle("Tracjetory of ld before Last Outcome") +
  scale_y_continuous(limits = c(100, 800), breaks= c(0:400*50)) + # Set tick every 4
  scale_x_continuous(limits = c(0, 60), breaks = c(0:60*10)) +
  theme(legend.position=c(0.77, 0.88)) + 
  ggtitle("Trajectory of LD") +
  theme(plot.title = element_text(hjust = 0.5))

##### AUC by incident cases/dynamic controls

## landmarkTimes for various time
units <- 5
landmarkTimes <- c(1:12)*units

## A. AUC_I/D
tableAUC_ID <- matrix(nrow=10, ncol=length(landmarkTimes))
tableAUC_TV_ID <- matrix(nrow=10, ncol=length(landmarkTimes))

# Baseline biomarkers at admission
  # 입원시 처음 측정한 값만으로 측정
scores <- c("plt", "lct", "wbc", "tb", "pf", "pt", "glc", "ld", "cr", "crp")
for(i in 1:length(scores)) {
  currVar <- eval(parse(text=paste("b_df$", scores[i],sep="")))  
  
  mmm <- MeanRank(survival.time= b_df$days, survival.status=b_df$mortality, marker= currVar)
  
  bandwidths <- 0.05 + c(1:80)/200
  IMSEs <- vector(length=length(bandwidths))
  for(j in 1:length(bandwidths)) {
    nnnC <- nne.CrossValidate(x= mmm$time, y= mmm$mean.rank, lambda=bandwidths[j])  #Cross-validated bandwidth
    IMSEs[j] <- nnnC$IMSE
  }
  
  currLambdaOS <- mean(bandwidths[which(IMSEs==min(IMSEs, na.rm=T))])
  
  nnn <- nne( x= mmm$time, y= mmm$mean.rank, lambda=currLambdaOS, nControls=mmm$nControls )  #Fixed bandwidth
  tableAUC_ID[i,] <- sapply(landmarkTimes, function(x){interpolate(x = nnn$x, y=nnn$nne, target=x)})
}
rownames(tableAUC_ID) <- scores
colnames(tableAUC_ID) <- landmarkTimes/units
round(tableAUC_ID, 3)

# Updated time-varying biomarkers 
# longitudinal하게 측정된 값을 이용한 경우
scores <- c("plt", "lct", "wbc", "tb", "pf", "pt", "glc", "ld", "cr", "crp")
for(i in 1:length(scores)) {
  currVar <- eval(parse(text=paste("dfl1$", scores[i],sep="")))  
  
  mmm <- MeanRank(survival.time= dfl1$tstop, survival.status= dfl1$deathEvent, marker= currVar, start= dfl1$tstart)
  
  bandwidths <- 0.05 + c(1:80)/200
  IMSEs <- vector(length=length(bandwidths))
  for(j in 1:length(bandwidths)) {
    nnnC <- nne.CrossValidate(x= mmm$time, y= mmm$mean.rank, lambda=bandwidths[j])  # Cross-validated bandwidth
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
  # c-index의 정확한 의미가 무엇인지 잘 모르겠습니다. AUC랑 무슨 차이가 있는지 code가 의미하는 것이 무엇인지?)
  # cutofftime을 설정하는 기준은 무엇인지?

# Baseline biomarkers at admission(입원시 처음 측정한 값만으로 측정)
round(dynamicIntegrateAUC(survival.time=b_df$days, survival.status=b_df$mortality, marker=b_df$plt, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time=b_df$days, survival.status=b_df$mortality, marker=b_df$lct, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time=b_df$days, survival.status=b_df$mortality, marker=b_df$wbc, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time=b_df$days, survival.status=b_df$mortality, marker=b_df$tb, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time=b_df$days, survival.status=b_df$mortality, marker=b_df$pf, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time=b_df$days, survival.status=b_df$mortality, marker=b_df$pt, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time=b_df$days, survival.status=b_df$mortality, marker=b_df$glc, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time=b_df$days, survival.status=b_df$mortality, marker=b_df$ld, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time=b_df$days, survival.status=b_df$mortality, marker=b_df$cr, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time=b_df$days, survival.status=b_df$mortality, marker=b_df$crp, cutoffTime = units*13), 3)

# Updated time-varying biomarkers (longitudinal하게 측정된 값을 이용한 경우)
round(dynamicIntegrateAUC(survival.time= dfl1$tstop, survival.status= dfl1$deathEvent, start=dfl1$tstart, marker=dfl1$plt, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time= dfl1$tstop, survival.status= dfl1$deathEvent, start=dfl1$tstart, marker=dfl1$lct, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time= dfl1$tstop, survival.status= dfl1$deathEvent, start=dfl1$tstart, marker=dfl1$wbc, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time= dfl1$tstop, survival.status= dfl1$deathEvent, start=dfl1$tstart, marker=dfl1$tb, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time= dfl1$tstop, survival.status= dfl1$deathEvent, start=dfl1$tstart, marker=dfl1$pf, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time= dfl1$tstop, survival.status= dfl1$deathEvent, start=dfl1$tstart, marker=dfl1$pt, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time= dfl1$tstop, survival.status= dfl1$deathEvent, start=dfl1$tstart, marker=dfl1$glc, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time= dfl1$tstop, survival.status= dfl1$deathEvent, start=dfl1$tstart, marker=dfl1$ld, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time= dfl1$tstop, survival.status= dfl1$deathEvent, start=dfl1$tstart, marker=dfl1$cr, cutoffTime = units*13), 3)
round(dynamicIntegrateAUC(survival.time= dfl1$tstop, survival.status= dfl1$deathEvent, start=dfl1$tstart, marker=dfl1$crp, cutoffTime = units*13), 3)

#### Bootstrap 95% CIs
nBoot <- 500

## A. Bootstrap CIs - Baseline markers
markers <- c("plt", "lct", "wbc", "tb", "pf", "pt", "glc", "ld", "cr", "crp")
set.seed(49)
Cindex_bstrap <- matrix(nrow=nBoot, ncol=length(markers))
bstrapRes <- list()

for(b in 1:nBoot) {
  currData <- b_df[sample(x=seq(1,nrow(b_df)), size=nrow(b_df), replace = T),]   
 
  kmfit <- survfit(Surv(days, mortality) ~ 1, data= currData)
  
  aucID_scores <- NULL

  for(i in 1:length(markers)) {
    currVar <- eval(parse(text=paste("currData$",markers[i],sep="")))
    
    ### AUC I/D
    mmm <- MeanRank(survival.time= currData$days, survival.status= currData$mortality, marker= currVar )
    nnn <- nne(x= mmm$time, y= mmm$mean.rank, lambda=0.3)  # Fixed bandwidth
    aucID_scores <- rbind(aucID_scores, 
                          sapply(landmarkTimes, function(t){interpolate( x = nnn$x, y=nnn$nne, target=t )}))
    
    ### C-index
    Cindex_bstrap[b,i] <- dynamicIntegrateAUC(survival.time=currData$days, survival.status= currData$mortality,
                                              marker=currVar, cutoffTime = units*13)
 }
  bstrapRes[[b]] <- list(aucID_scores=aucID_scores)
}

# Get CIs for c-indices
Cindex_CIs <- round(apply(Cindex_bstrap, 2, quantile, probs=c(0.025,0.975)),3)
colnames(Cindex_CIs) <- markers
Cindex_CIs

# Get CIs for AUCs
AUC_ID_CIs <- NULL

for(t in 1:length(landmarkTimes)) {
  AUC_ID <- NULL
 
  for(b in 1:nBoot) {
    AUC_ID <- cbind( AUC_ID, bstrapRes[[b]]$aucID_scores[,t] )
   }
  AUC_ID_CI_raw <- round( apply(AUC_ID, 1, quantile, probs=c(0.025,0.975)), 3 )
 
  AUC_ID_CIs <- cbind(AUC_ID_CIs, sapply(seq(10), function(x) paste("(", AUC_ID_CI_raw[1,x], ", ", AUC_ID_CI_raw[2,x], ")", sep="" ) ) )
}
rownames(AUC_ID_CIs) <- c("plt", "lct", "wbc", "tb", "pf", "pt", "glc", "ld", "cr", "crp")
colnames(AUC_ID_CIs) <- landmarkTimes/units
AUC_ID_CIs


#### baseline AUC plot...rownames(tableAUC_ID) <- scores
a<-round(tableAUC_ID, 3)
b <- AUC_ID_CIs
c<-data.frame(marker = c("plt", "lct", "wbc", "tb", "pf", "pt", "glc", "ld", "cr", "crp"),
              matrix(b, nrow=10)) 
c
?matrix
b

me(AUC_ID_CIs)
b<-data.frame(a)
d<- c %>% 
  mutate(str_extract_all("[0-9]+\\.[0-9]+"))%>%
  unnest(X1) %>% 
  group_by(marker) %>%
  mutate(x1 = paste0("up", row_number())) %>%
  spread(up, x1) %>% ungroup()
d
utate()

str(b)

## AUC plot 1 ##
Au1 <- left_join(data.frame(tableROC) %>% 
                   mutate(marker = rownames(.)) %>%
                   gather(., hd.bwd, auc, 1:13, factor_key=TRUE), 
                 data.frame(tableCiUp) %>% 
                   mutate(marker = rownames(.)) %>%
                   gather(., hd.bwd, up, 1:13, factor_key=TRUE), 
                 by=c("marker", "hd.bwd")) %>%
  left_join(., data.frame(tableCiDn) %>%
              mutate(marker = rownames(.)) %>%
              gather(., hd.bwd, dn, 1:13, factor_key=TRUE), 
            by=c("marker", "hd.bwd")) %>%
  mutate(hd.bwd=-1*as.numeric(str_extract(hd.bwd,"[0-9]+"))) %>%
  filter(marker=="WBC" | marker=="Bilirubin" | marker=="Glucose" | marker=="LD" | marker=="CRP") %>%
  mutate(marker=factor(marker, levels = c("CRP", "Bilirubin", "LD", "Glucose", "WBC")))

theme_set(theme_bw())
ggplot(Au1, mapping = aes(x = as.numeric(as.character(hd.bwd)), y = auc, 
                          group = marker,
                          color = marker)) +
  geom_errorbar(aes(ymin=dn, ymax=up), width=5, position = position_dodge(2)) +
  geom_line(position=position_dodge(2)) +
  geom_point(position=position_dodge(2)) +
  xlab("Days before Outcome") +
  ylab("AUC") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  scale_y_continuous(limits = c(0.3, 1), breaks= c(0:10*0.1)) + # Set tick every 4
  scale_x_continuous(limits = c(-62, 2), breaks = c(-60:0*10)) +
  theme_bw() +
  theme(legend.position="right") + 
  ggtitle("AUC(95% CI) of Biomarkers") +
  theme(plot.title = element_text(hjust = 0.5))

## AUC plot 2 ##
Au2 <- left_join(data.frame(tableROC) %>% 
                   mutate(marker = rownames(.)) %>%
                   gather(., hd.bwd, auc, 1:13, factor_key=TRUE), 
                 data.frame(tableCiUp) %>% 
                   mutate(marker = rownames(.)) %>%
                   gather(., hd.bwd, up, 1:13, factor_key=TRUE), 
                 by=c("marker", "hd.bwd")) %>%
  left_join(., data.frame(tableCiDn) %>%
              mutate(marker = rownames(.)) %>%
              gather(., hd.bwd, dn, 1:13, factor_key=TRUE), 
            by=c("marker", "hd.bwd")) %>%
  mutate(hd.bwd=-1*as.numeric(str_extract(hd.bwd,"[0-9]+"))) %>%
  filter(marker=="Platelet" | marker=="Lactate" | marker=="Prothrombin time" | marker=="Creatinine" | marker=="PF ratio") %>%
  mutate(marker=factor(marker, levels = c("Platelet", "Prothrombin time", "Lactate", "Creatinine", "PF ratio")))

str(Au2)
theme_set(theme_bw())
ggplot(Au2, mapping = aes(x = as.numeric(as.character(hd.bwd)), y = auc, 
                          group = marker,
                          color = marker)) +
  geom_errorbar(aes(ymin=dn, ymax=up), width=5, position = position_dodge(2)) +
  geom_line(position=position_dodge(2)) +
  geom_point(position=position_dodge(2)) +
  xlab("Days before Outcome") +
  ylab("AUC") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  scale_y_continuous(limits = c(0.3, 1), breaks= c(0:10*0.1)) + # Set tick every 4
  scale_x_continuous(limits = c(-62, 2), breaks = c(-60:0*10)) +
  theme_bw() +
  theme(legend.position="right") + 
  ggtitle("AUC(95% CI) of Biomarkers") +
  theme(plot.title = element_text(hjust = 0.5))









## B. Bootstrap CIs - Time-varying scores ####
markers <- c("plt", "lct", "wbc", "tb", "pf", "pt", "glc", "ld", "cr", "crp")

set.seed(49)
Cindex_bstrapTV <- matrix(nrow=nBoot, ncol=length(markers))
bstrapResTV <- list()

for(b in 1:nBoot) {
  #sample individuals
  subjs <- unique(dfl1$id)
  currSubjs <- sample(x=subjs, size=length(subjs), replace = T)
  currData <- NULL
  for(j in 1:length(currSubjs)) 
    currData <- rbind(currData, dfl1[which(dfl1$id==currSubjs[j]),]) 
  
  kmfit <- survfit(Surv(time=tstart, time2=tstop, event=deathEvent) ~ 1, data=currData)
  
  aucID_scores <- NULL

  for(i in 1:length(markers)) {
    currVar <- eval(parse(text=paste("currData$", markers[i],sep="")))
    
    ### AUC I/D TV
    mmm <- MeanRank( survival.time= currData$tstop, survival.status= currData$deathEvent, start=currData$tstart, marker= currVar )
    nnn <- nne( x= mmm$time, y= mmm$mean.rank, lambda=0.3, nControls=mmm$nControls )  #Fixed bandwidth
    aucID_scores <- rbind(aucID_scores, 
                          sapply(landmarkTimes, function(t){ interpolate( x = nnn$x, y=nnn$nne, target=t ) } ))
    
    ### C-index TV
    Cindex_bstrapTV[b,i] <- dynamicIntegrateAUC(survival.time=currData$tstop, survival.status=currData$deathEvent,
                                                start=currData$tstart, marker= currVar, cutoffTime = units*13)
  }
  bstrapResTV[[b]] <- list(aucID_scores)
}

# Get CIs for c-indices
Cindex_CIsTV <- round(apply(Cindex_bstrapTV, 2, quantile, probs=c(0.025,0.975)), 3)
colnames(Cindex_CIsTV) <- markers
Cindex_CIsTV

# Get CIs for AUCs
AUC_ID_CIsTV <- NULL

for(t in 1:length(landmarkTimes)) {
  AUC_ID <- NULL
  for(b in 1:nBoot) {
    AUC_ID <- cbind( AUC_ID, bstrapResTV[[b]][[1]][,t] )
   }
  AUC_ID_CI_raw <- round( apply(AUC_ID, 1, quantile, probs=c(0.025,0.975)), 3 )
  
  AUC_ID_CIsTV <- cbind(AUC_ID_CIsTV, sapply(seq(10), function(x) paste("(", AUC_ID_CI_raw[1,x], ", ", AUC_ID_CI_raw[2,x], ")", sep="" ) ) )
}
rownames(AUC_ID_CIsTV) <- c("plt", "lct", "wbc", "tb", "pf", "pt", "glc", "ld", "cr", "crp")
colnames(AUC_ID_CIsTV) <- landmarkTimes/units
AUC_ID_CIsTV

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

getCindexBstrapCI(nBoot=500, inData=b_df, markerVarName1="plt", markerVarName2="glc",
                  timeVarName="days", eventVarName="mortality", cutoffTime=units*13)


### FIGURES

# Figure 2
# pdf("Figure2_meanRank_pbc_baseline.pdf", height=4, width=9)
# par(mfrow=c(1,2))
# par(ps=10)

# AUC I/D for baseline biomarkers
mmmBaseline <- MeanRank(survival.time= b_df$days, survival.status= b_df$mortality,
                         marker= b_df$plt)
print(length(mmmBaseline$time))
nnnBaseline <- nne(x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.2, nControls= mmmBaseline$nControls)  #Fixed bandwidth
plot(mmmBaseline$time, mmmBaseline$mean.rank, xlab="Time (Days)", ylab=expression(AUC^"I/D"*"(t)"), 
      ylim=c(0,1), col="blue", pch=21, cex=.8, axes=F, xlim=c(0,60))
axis(1, at=seq(0,60,by=5), labels=seq(0,60,by=5))
axis(2)
box()
abline(h=0.5, lty=2)
lines(nnnBaseline$x, nnnBaseline$nne, col="blue", lwd=2 )

mmmBaseline <- MeanRank( survival.time= b_df$days, survival.status= b_df$mortality,
                         marker= b_df$lct )
print(length(mmmBaseline$time))
nnnBaseline <- nne( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.2, nControls= mmmBaseline$nControls )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="orange", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="orange", lwd=2 )

mmmBaseline <- MeanRank( survival.time= b_df$days, survival.status= b_df$mortality, 
                         marker= b_df$wbc )
print(length(mmmBaseline$time))
nnnBaseline <- nne( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.2, nControls= mmmBaseline$nControls )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="red", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="red", lwd=2 )

mmmBaseline <- MeanRank( survival.time= b_df$days, survival.status= b_df$mortality,
                         marker= b_df$tb )
print(length(mmmBaseline$time))
nnnBaseline <- nne( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.2, nControls= mmmBaseline$nControls )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="green", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="green", lwd=2 )

mmmBaseline <- MeanRank( survival.time= b_df$days, survival.status= b_df$mortality,
                         marker= b_df$pf )
print(length(mmmBaseline$time))
nnnBaseline <- nne( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.2, nControls= mmmBaseline$nControls )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="purple", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="purple", lwd=2 )

mmmBaseline <- MeanRank( survival.time= b_df$days, survival.status= b_df$mortality,
                         marker= b_df$pt )
print(length(mmmBaseline$time))
nnnBaseline <- nne( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.2, nControls= mmmBaseline$nControls )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="saddlebrown", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="saddlebrown", lwd=2 ,lty=2)

mmmBaseline <- MeanRank( survival.time= b_df$days, survival.status= b_df$mortality, 
                         marker= b_df$glc )
print(length(mmmBaseline$time))
nnnBaseline <- nne( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.2, nControls= mmmBaseline$nControls )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="seagreen3", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="seagreen3", lwd=2 ,lty=2)

mmmBaseline <- MeanRank( survival.time= b_df$days, survival.status= b_df$mortality, 
                         marker= b_df$ld )
print(length(mmmBaseline$time))
nnnBaseline <- nne( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.2, nControls= mmmBaseline$nControls )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="gray", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="gray", lwd=2 ,lty=2)

mmmBaseline <- MeanRank( survival.time= b_df$days, survival.status= b_df$mortality, 
                         marker= b_df$cr )
print(length(mmmBaseline$time))
nnnBaseline <- nne( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.2, nControls= mmmBaseline$nControls )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="orangered", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="orangered", lwd=2 ,lty=2)

mmmBaseline <- MeanRank( survival.time= b_df$days, survival.status= b_df$mortality,
                         marker= b_df$crp )
print(length(mmmBaseline$time))
nnnBaseline <- nne( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.2, nControls= mmmBaseline$nControls )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="purple4", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="purple4", lwd=2 ,lty=2)

# ROC I/D (TPR) -> 이 그래프가 의미하는 것은 무엇인지요? AUC의 변화만 보여주고 이것은 빼도 되는거 아닌지 궁금합니다.
fpr <- 0.1

mmmBaseline <- dynamicTP( p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$plt )
print(length(mmmBaseline$time))
nnnBaseline <- nne_TPR( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda= 0.3, nControls= mmmBaseline$nControls,
                        nCases= mmmBaseline$nCases, p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$score5baseline )  #Fixed bandwidth
plot( mmmBaseline$time, mmmBaseline$mean.rank, xlab="Time (years)", ylab=expression(ROC[t]^"I/D"*"(FPF=10%)"), ylim=c(0,1), 
      col="blue", pch=21, cex=.8, axes=F, xlim=c(0,200))
axis(1, at=seq(0,200,by=5), labels=seq(0,200,by=5))
axis(2)
box()
abline(h=0.5, lty=2)
lines( nnnBaseline$x, nnnBaseline$nne, col="blue", lwd=2 )

mmmBaseline <- dynamicTP( p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$lct )
print(length(mmmBaseline$time))
nnnBaseline <- nne_TPR( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.3, nControls= mmmBaseline$nControls,
                        nCases= mmmBaseline$nCases, p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$score4baseline )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="orange", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="orange", lwd=2 )

mmmBaseline <- dynamicTP( p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$wbc )
print(length(mmmBaseline$time))
nnnBaseline <- nne_TPR( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.3, nControls= mmmBaseline$nControls,
                        nCases= mmmBaseline$nCases, p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$score4baseline )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="red", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="red", lwd=2 )

mmmBaseline <- dynamicTP( p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$tb )
print(length(mmmBaseline$time))
nnnBaseline <- nne_TPR( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.3, nControls= mmmBaseline$nControls,
                        nCases= mmmBaseline$nCases, p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$score4baseline )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="green", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="green", lwd=2 )

mmmBaseline <- dynamicTP( p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$pf )
print(length(mmmBaseline$time))
nnnBaseline <- nne_TPR( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.3, nControls= mmmBaseline$nControls,
                        nCases= mmmBaseline$nCases, p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$score4baseline )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="purple", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="purple", lwd=2 )

mmmBaseline <- dynamicTP( p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$pt )
print(length(mmmBaseline$time))
nnnBaseline <- nne_TPR( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.3, nControls= mmmBaseline$nControls,
                        nCases= mmmBaseline$nCases, p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$score4baseline )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="cyan", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="cyan", lwd=2 )

mmmBaseline <- dynamicTP( p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$glc )
print(length(mmmBaseline$time))
nnnBaseline <- nne_TPR( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.3, nControls= mmmBaseline$nControls,
                        nCases= mmmBaseline$nCases, p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$score4baseline )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="darkorchid1", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="darkorchid1", lwd=2 )

mmmBaseline <- dynamicTP( p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$ld )
print(length(mmmBaseline$time))
nnnBaseline <- nne_TPR( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.3, nControls= mmmBaseline$nControls,
                        nCases= mmmBaseline$nCases, p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$score4baseline )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="darkseagreen4", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="darkseagreen4", lwd=2 )

mmmBaseline <- dynamicTP( p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$cr )
print(length(mmmBaseline$time))
nnnBaseline <- nne_TPR( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.3, nControls= mmmBaseline$nControls,
                        nCases= mmmBaseline$nCases, p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$score4baseline )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="deeppink4", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="deeppink4", lwd=2 )

mmmBaseline <- dynamicTP( p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$crp )
print(length(mmmBaseline$time))
nnnBaseline <- nne_TPR( x= mmmBaseline$time, y= mmmBaseline$mean.rank, lambda=0.3, nControls= mmmBaseline$nControls,
                        nCases= mmmBaseline$nCases, p=fpr, survival.time= b_df$days, survival.status= b_df$mortality, marker= b_df$score4baseline )  #Fixed bandwidth
points( mmmBaseline$time, mmmBaseline$mean.rank, col="gray0", pch=21, cex=.8 )
lines( nnnBaseline$x, nnnBaseline$nne, col="gray0", lwd=2 )

legend(x=1.5*units, y=0.9, legend=c("4 covariates", "5 covariates"), col=c("orange","blue"), lty=1, lwd=2, horiz=T)


# Figure
mmm <- MeanRank(survival.time= dfl1$tstop, survival.status= dfl1$deathEvent, marker= currVar, start= dfl1$tstart)

# pdf("Figure3_meanRank_pbc_TV.pdf", height=4, width=9)
# par(mfrow=c(1,2))
# par(ps=10)

# AUC I/D for updated time varing biomarkers
mmmTV <- MeanRank(survival.time= dfl1$tstop, survival.status= dfl1$deathEvent, 
                  marker= dfl1$plt, start= dfl1$tstart)
print(length(mmmTV$time))
nnn <- nne( x= mmmTV$time, y= mmmTV$mean.rank, lambda=0.2, nControls=mmmTV$nControls )  #Fixed bandwidth
plot( mmmTV$time, mmmTV$mean.rank, xlab="Time (years)", ylab=expression(AUC["I/D"]*"(t)"), ylim=c(0.4,1), 
      col="blue", pch=21, cex=.8, axes=F, xlim=c(0,200))
axis(1, at=seq(0,200,by=5), labels=seq(0,200,by=5))
axis(2)
box()
abline(h=0.5, lty=2)
lines( nnn$x, nnn$nne, col="blue", lwd=2, lty=1 )

mmmTV <- MeanRank( survival.time= dfl1$tstop, survival.status= dfl1$deathEvent,
                   marker= dfl1$lct, start= dfl1$tstart )
print(length(mmmTV$time))
nnn <- nne( x= mmmTV$time, y= mmmTV$mean.rank, lambda=0.2, nControls=mmmTV$nControls )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="orange", pch=21, cex=.8)
lines( nnn$x, nnn$nne, col="orange", lwd=2, lty=1 )

mmmTV <- MeanRank( survival.time= dfl1$tstop, survival.status= dfl1$deathEvent,
                   marker= dfl1$wbc, start= dfl1$tstart )
print(length(mmmTV$time))
nnn <- nne( x= mmmTV$time, y= mmmTV$mean.rank, lambda=0.2, nControls=mmmTV$nControls )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="red", pch=21, cex=.8)
lines( nnn$x, nnn$nne, col="red", lwd=2, lty=1 )

mmmTV <- MeanRank( survival.time= dfl1$tstop, survival.status= dfl1$deathEvent,
                   marker= dfl1$tb, start= dfl1$tstart )
print(length(mmmTV$time))
nnn <- nne( x= mmmTV$time, y= mmmTV$mean.rank, lambda=0.2, nControls=mmmTV$nControls )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="green", pch=21, cex=.8)
lines( nnn$x, nnn$nne, col="green", lwd=2, lty=1 )

mmmTV <- MeanRank( survival.time= dfl1$tstop, survival.status= dfl1$deathEvent,
                   marker= dfl1$pf, start= dfl1$tstart )
print(length(mmmTV$time))
nnn <- nne( x= mmmTV$time, y= mmmTV$mean.rank, lambda=0.2, nControls=mmmTV$nControls )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="purple", pch=21, cex=.8)
lines( nnn$x, nnn$nne, col="purple", lwd=2, lty=1 )

mmmTV <- MeanRank( survival.time= dfl1$tstop, survival.status= dfl1$deathEvent,
                   marker= dfl1$pt, start= dfl1$tstart )
print(length(mmmTV$time))
nnn <- nne( x= mmmTV$time, y= mmmTV$mean.rank, lambda=0.2, nControls=mmmTV$nControls )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="saddlebrown", pch=21, cex=.8)
lines( nnn$x, nnn$nne, col="saddlebrown", lwd=2, lty=2 )

mmmTV <- MeanRank( survival.time= dfl1$tstop, survival.status= dfl1$deathEvent,
                   marker= dfl1$glc, start= dfl1$tstart )
print(length(mmmTV$time))
nnn <- nne( x= mmmTV$time, y= mmmTV$mean.rank, lambda=0.2, nControls=mmmTV$nControls )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="seagreen3", pch=21, cex=.8)
lines( nnn$x, nnn$nne, col="seagreen3", lwd=2, lty=2 )

mmmTV <- MeanRank( survival.time= dfl1$tstop, survival.status= dfl1$deathEvent, 
                   marker= dfl1$ld, start= dfl1$tstart )
print(length(mmmTV$time))
nnn <- nne( x= mmmTV$time, y= mmmTV$mean.rank, lambda=0.2, nControls=mmmTV$nControls )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="plum", pch=21, cex=.8)
lines( nnn$x, nnn$nne, col="plum", lwd=2, lty=2 )

mmmTV <- MeanRank( survival.time= dfl1$tstop, survival.status= dfl1$deathEvent,
                   marker= dfl1$cr, start= dfl1$tstart )
print(length(mmmTV$time))
nnn <- nne( x= mmmTV$time, y= mmmTV$mean.rank, lambda=0.2, nControls=mmmTV$nControls )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="orangered", pch=21, cex=.8)
lines( nnn$x, nnn$nne, col="orangered", lwd=2, lty=2 )

mmmTV <- MeanRank( survival.time= dfl1$tstop, survival.status= dfl1$deathEvent,
                   marker= dfl1$crp, start= dfl1$tstart )
print(length(mmmTV$time))
nnn <- nne( x= mmmTV$time, y= mmmTV$mean.rank, lambda=0.2, nControls=mmmTV$nControls )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="purple4", pch=21, cex=.8)
lines( nnn$x, nnn$nne, col="purple4", lwd=2, lty=2 )

# ROC I/D (TPR)
mmmTV <- dynamicTP( p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$plt, start= dfl1$tstart )
print(length(mmmTV$time))
nnn <- nne_TPR( x= mmmTV$time, y= mmmTV$mean.rank, lambda= 0.3, nControls= mmmTV$nControls,
                nCases= mmmTV$nCases, p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$plt, start= dfl1$tstart )  #Fixed bandwidth
plot( mmmTV$time, mmmTV$mean.rank, xlab="Time (years)", ylab=expression(ROC[t]^"I/D"*"(FPF=10%)"), ylim=c(0,1), col="blue", pch=21, cex=.8, axes=F, xlim=c(0,11)*units)
axis(1, at=seq(0,60,by=5), labels=seq(0,60,by=5))
axis(2)
box()
abline(h=0.5, lty=2)
lines( nnn$x, nnn$nne, col="blue", lwd=2 )

mmmTV <- dynamicTP( p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$lct, start= dfl1$tstart )
nnn <- nne_TPR( x= mmmTV$time, y= mmmTV$mean.rank, lambda= 0.3, nControls= mmmTV$nControls,
                nCases= mmmTV$nCases, p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$lct, start= dfl1$tstart )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="orange", pch=21, cex=.8 )
lines( nnn$x, nnn$nne, col="orange", lwd=2 )

mmmTV <- dynamicTP( p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$wbc, start= dfl1$tstart )
nnn <- nne_TPR( x= mmmTV$time, y= mmmTV$mean.rank, lambda= 0.3, nControls= mmmTV$nControls,
                nCases= mmmTV$nCases, p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$wbc, start= dfl1$tstart )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="red", pch=21, cex=.8 )
lines( nnn$x, nnn$nne, col="red", lwd=2 )

mmmTV <- dynamicTP( p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$tb, start= dfl1$tstart )
nnn <- nne_TPR( x= mmmTV$time, y= mmmTV$mean.rank, lambda= 0.3, nControls= mmmTV$nControls,
                nCases= mmmTV$nCases, p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$tb, start= dfl1$tstart )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="green", pch=21, cex=.8 )
lines( nnn$x, nnn$nne, col="green", lwd=2 )

mmmTV <- dynamicTP( p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$pf, start= dfl1$tstart )
nnn <- nne_TPR( x= mmmTV$time, y= mmmTV$mean.rank, lambda= 0.3, nControls= mmmTV$nControls,
                nCases= mmmTV$nCases, p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$pf, start= dfl1$tstart )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="purple", pch=21, cex=.8 )
lines( nnn$x, nnn$nne, col="purple", lwd=2 )

mmmTV <- dynamicTP( p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$pt, start= dfl1$tstart )
nnn <- nne_TPR( x= mmmTV$time, y= mmmTV$mean.rank, lambda= 0.3, nControls= mmmTV$nControls,
                nCases= mmmTV$nCases, p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$pt, start= dfl1$tstart )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="saddlebrown", pch=21, cex=.8 )
lines( nnn$x, nnn$nne, col="saddlebrown", lwd=2, lty=2 )

mmmTV <- dynamicTP( p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$glc, start= dfl1$tstart )
nnn <- nne_TPR( x= mmmTV$time, y= mmmTV$mean.rank, lambda= 0.3, nControls= mmmTV$nControls,
                nCases= mmmTV$nCases, p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$glc, start= dfl1$tstart )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="seagreen3", pch=21, cex=.8 )
lines( nnn$x, nnn$nne, col="seagreen3", lwd=2, lty=2 )

mmmTV <- dynamicTP( p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$ld, start= dfl1$tstart )
nnn <- nne_TPR( x= mmmTV$time, y= mmmTV$mean.rank, lambda= 0.3, nControls= mmmTV$nControls,
                nCases= mmmTV$nCases, p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$ld, start= dfl1$tstart )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="plum", pch=21, cex=.8 )
lines( nnn$x, nnn$nne, col="plum", lwd=2, lty=2)

mmmTV <- dynamicTP( p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$cr, start= dfl1$tstart )
nnn <- nne_TPR( x= mmmTV$time, y= mmmTV$mean.rank, lambda= 0.3, nControls= mmmTV$nControls,
                nCases= mmmTV$nCases, p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$cr, start= dfl1$tstart )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="orangered", pch=21, cex=.8 )
lines( nnn$x, nnn$nne, col="orangered", lwd=2, lty=2 )

mmmTV <- dynamicTP( p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$crp, start= dfl1$tstart )
nnn <- nne_TPR( x= mmmTV$time, y= mmmTV$mean.rank, lambda= 0.3, nControls= mmmTV$nControls,
                nCases= mmmTV$nCases, p=fpr, survival.time =dfl1$tstop, survival.status= dfl1$deathEvent, marker= dfl1$crp, start= dfl1$tstart )  #Fixed bandwidth
points( mmmTV$time, mmmTV$mean.rank, col="purple4", pch=21, cex=.8 )
lines( nnn$x, nnn$nne, col="purple4", lwd=2, lty=2 )

# legend(x=2.5*units, y=0.3, legend=c("4 covariates", "5 covariates"), col=c("orange","blue"), lty=1, lwd=2, horiz=T)

## save.image(file = "TV fwd.RData")

