setwd("F:\\00___Paper_Work\\01_TV_analysis")

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

### 전체적으로 생존에 영향을 미치는 여부 보기
# formula가 맞는지 모르겠습니다.
glmm.bwd <- glmer(death ~ age+tbsa+inh+tb*hd.bwd+wbc*hd.bwd+pt*hd.bwd+tb*hd.bwd+cr*hd.bwd+crp*hd.bwd+glc*hd.bwd+lct*hd.bwd+ld*hd.bwd+pf*hd.bwd
                  + (1+hd.bwd|id), data=df, family = binomial)

summary(glmm.bwd)
print(glmm.bwd, correlation=T)
## 결과 해석은 어떻게 해야 하는지?

### 각각의 검사가 생존과 사망에 따라 시간에 따라 차이가 있는지 분석 
nlme_tb <- lme(tb ~ hd.bwd*death,  random = ~hd.bwd|id, data = df, na.action = na.omit)
nlme_lct <- lme(lct ~ hd.bwd*death,  random = ~hd.bwd|id, data = df, na.action = na.omit)
nlme_pt <- lme(pt ~ hd.bwd*death,  random = ~hd.bwd|id, data = df, na.action = na.omit)
nlme_tb <- lme(tb ~ hd.bwd*death,  random = ~hd.bwd|id, data = df, na.action = na.omit)
nlme_wbc <- lme(wbc ~ hd.bwd*death,  random = ~hd.bwd|id, data = df, na.action = na.omit)
nlme_crp <- lme(crp ~ hd.bwd*death,  random = ~hd.bwd|id, data = df, na.action = na.omit)
nlme_cr <- lme(cr ~ hd.bwd*death,  random = ~hd.bwd|id, data = df, na.action = na.omit)
nlme_ld <- lme(ld ~ hd.bwd*death,  random = ~hd.bwd|id, data = df, na.action = na.omit)
nlme_glc <- lme(glc ~ hd.bwd*death,  random = ~hd.bwd|id, data = df, na.action = na.omit)
nlme_pf <- lme(pf ~ hd.bwd*death,  random = ~hd.bwd|id, data = df, na.action = na.omit)

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

## Biomarker 평균 변화 plot ####

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

#### plt ##
plot_plt <- summarySE(df, measurevar="plt", groupvars=c("hd.bwd","mortality"), na.rm = T)
 
theme_set(theme_bw())
ggplot(plot_plt, mapping = aes(x = as.numeric(as.character(hd.bwd)), y = plt, 
                         group = mortality,
                         color = mortality)) +
  geom_errorbar(aes(ymin=plt-nci, ymax=plt+nci), width=0.7) +
  geom_line() +
  geom_point() +
  xlab("Days before Outcome") +
  ylab("Platelet") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  scale_y_continuous(limits = c(0, 450), breaks= c(0:550*50)) + # Set tick every 4
  scale_x_continuous(limits = c(-61, 1), breaks = c(-60:0*10)) +
  theme_bw() +
  theme(legend.position=c(0.77, 0.88)) + 
  ggtitle("Trajectory of Platelet") +
  theme(plot.title = element_text(hjust = 0.5))

#### pt ##
plot_pt <- summarySE(df, measurevar="pt", groupvars=c("hd.bwd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_pt, mapping = aes(x = as.numeric(as.character(hd.bwd)), y = pt, 
                         group = mortality,
                         color = mortality)) +
  geom_errorbar(aes(ymin=pt-nci, ymax=pt+nci), width=0.7) +
  geom_line() +
  geom_point() +
  xlab("Days before Outcome") +
  ylab("Prothrombin Time") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  scale_y_continuous(limits = c(11, 24), breaks= c(0:20*2)) + # Set tick every 4
  scale_x_continuous(limits = c(-61, 1), breaks = c(-60:0*10)) +
  theme_bw() +
  theme(legend.position=c(0.77, 0.88)) + 
  ggtitle("Trajectory of Prothrombin Time") +
  theme(plot.title = element_text(hjust = 0.5))

#### tb ####
plot_tb <- summarySE(df, measurevar="tb", groupvars=c("hd.bwd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_tb, mapping = aes(x = as.numeric(as.character(hd.bwd)), y = tb, 
                         group = mortality,
                         color = mortality)) +
  geom_errorbar(aes(ymin=tb-nci, ymax=tb+nci), width=0.7) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  xlab("Days before Outcome") +
  ylab("Bilirubin") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  scale_y_continuous(limits = c(0, 5), breaks= c(0:8*1)) + # Set tick every 4
  scale_x_continuous(limits = c(-61, 1), breaks = c(-60:0*10)) +
  theme_bw() +
  theme(legend.position=c(0.77, 0.88)) + 
  ggtitle("Trajectory of Bilirubin") +
  theme(plot.title = element_text(hjust = 0.5))

#### lct ####
plt_lct <- summarySE(df, measurevar="lct", groupvars=c("hd.bwd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plt_lct, mapping = aes(x = as.numeric(as.character(hd.bwd)), y = lct, 
                         group = mortality,
                         color = mortality)) +
  geom_errorbar(aes(ymin=lct-nci, ymax=lct+nci), width=0.7) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  xlab("Days before Outcome") +
  ylab("Lactate") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  scale_y_continuous(limits = c(0.5, 5), breaks= c(0:25*1)) + 
  scale_x_continuous(limits = c(-61, 1), breaks = c(-60:0*10)) +
  theme_bw() +
  theme(legend.position=c(0.77, 0.88)) + 
  ggtitle("Trajectory of Lactate") +
  theme(plot.title = element_text(hjust = 0.5))

#### cr ####
plot_cr <- summarySE(df, measurevar="cr", groupvars=c("hd.bwd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_cr, mapping = aes(x = as.numeric(as.character(hd.bwd)), y = cr, 
                         group = mortality,
                         color = mortality)) +
  geom_errorbar(aes(ymin=cr-nci, ymax=cr+nci), width=0.7) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  xlab("Days before Outcome") +
  ylab("Creatinine") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  scale_y_continuous(limits = c(0.2, 2.5), breaks= c(0:5*0.5)) + # Set tick every 4
  scale_x_continuous(limits = c(-61, 1), breaks = c(-60:0*10)) +
  theme_bw() +
  theme(legend.position=c(0.77, 0.88)) + 
  ggtitle("Trajectory of Creatinine") +
  theme(plot.title = element_text(hjust = 0.5))

#### glc ####
plot_glc <- summarySE(df, measurevar="glc", groupvars=c("hd.bwd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_glc, mapping = aes(x = as.numeric(as.character(hd.bwd)), y = glc, 
                         group = mortality,
                         color = mortality)) +
  geom_errorbar(aes(ymin=glc-nci, ymax=glc+nci), width=0.7) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  xlab("Days before Outcome") +
  ylab("Glucose") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  scale_y_continuous(limits = c(120, 300), breaks= c(0:400*50)) +
  scale_x_continuous(limits = c(-61, 1), breaks = c(-60:0*10)) +
  theme_bw() +
  theme(legend.position=c(0.8, 0.85)) + 
  ggtitle("Trajectory of Glucose") +
  theme(plot.title = element_text(hjust = 0.5))

#### crp ####
plot_crp <- summarySE(df, measurevar="crp", groupvars=c("hd.bwd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_crp, mapping = aes(x = as.numeric(as.character(hd.bwd)), y = crp, 
                         group = mortality,
                         color = mortality)) +
  geom_errorbar(aes(ymin=crp-nci, ymax=crp+nci), width=0.7) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  xlab("Days before Outcome") +
  ylab("C-reative Protein") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  scale_y_continuous(limits = c(30, 250), breaks= c(0:400*50)) + # Set tick every 4
  scale_x_continuous(limits = c(-61, 1), breaks = c(-60:0*10)) +
  theme_bw() +
  theme(legend.position=c(0.8, 0.9)) +   
  ggtitle("Trajectory of C-reative Protein") +
  theme(plot.title = element_text(hjust = 0.5))

#### wbc ####
plot_wbc <- summarySE(df, measurevar="wbc", groupvars=c("hd.bwd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_wbc, mapping = aes(x = as.numeric(as.character(hd.bwd)), y = wbc, 
                         group = mortality,
                         color = mortality)) +
  # geom_smooth(method="lm",alpha=0.5,aes(group=mortality))+
  geom_errorbar(aes(ymin=wbc-nci, ymax=wbc+nci), width=0.7) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  xlab("Days before Outcome") +
  ylab("WBC") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  scale_y_continuous(limits = c(8, 20), breaks= c(0:40*2)) + # Set tick every 4
  scale_x_continuous(limits = c(-61, 1), breaks = c(-60:0*10)) +
  theme_bw() +
  theme(legend.position=c(0.78, 0.9)) + 
  ggtitle("Trajectory of WBC") +
  theme(plot.title = element_text(hjust = 0.5))

#### pf ####
plot_pf <- summarySE(df, measurevar="pf", groupvars=c("hd.bwd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_pf, mapping = aes(x = as.numeric(as.character(hd.bwd)), y = pf, 
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
  scale_y_continuous(limits = c(150, 350), breaks= c(0:400*50)) + # Set tick every 4
  scale_x_continuous(limits = c(-61, 1), breaks = c(-60:0*10)) +
  theme(legend.position=c(0.75, 0.88)) + 
  ggtitle("Trajectory of PF ratio") +
  theme(plot.title = element_text(hjust = 0.5))

#### ld ####
plot_ld <- summarySE(df, measurevar="ld", groupvars=c("hd.bwd","mortality"), na.rm = T)

theme_set(theme_bw())
ggplot(plot_ld, mapping = aes(x = as.numeric(as.character(hd.bwd)), y = ld, 
                              group = mortality,
                              color = mortality)) +
  geom_errorbar(aes(ymin=ld-nci, ymax=ld+nci), width=0.7) +
  geom_line() +
  geom_point() +
  xlab("Day before Outcome") +
  ylab("LD") +
  scale_colour_hue(name="Mortaltiy",    # Legend label, use darker colors
                   l=40) +  # Use darker colors, lightness=40
  ggtitle("Tracjetory of ld before Last Outcome") +
  scale_y_continuous(limits = c(100, 750), breaks= c(0:400*50)) + # Set tick every 4
  scale_x_continuous(limits = c(-61, 1), breaks = c(-60:0*10)) +
  theme(legend.position=c(0.77, 0.88)) + 
  ggtitle("Trajectory of LD") +
  theme(plot.title = element_text(hjust = 0.5))


## 이 부분은 각 시점에서 검사로 마지막날의 사망율을 예측하는 AUC 구하기 
  ## 시간에 따라서 변화하는 것이 아니고 그 시점에서 에측이므로 time dependent ROC가 아닌 일반 ROC 로 구함.

library(pROC)

## landmarkTimes for various time
units <- -5
landmarkTimes <- c(12:0)*units

tableROC <- matrix(nrow=10, ncol=length(landmarkTimes))
tableCiUp <- matrix(nrow=10, ncol=length(landmarkTimes))
tableCiDn <- matrix(nrow=10, ncol=length(landmarkTimes))

for(j in 1:length(landmarkTimes)) {
  currDataTV <- subset(df, hd.bwd == (landmarkTimes[j]))
  
  out1 <- auc(currDataTV$mortality ~ currDataTV$plt) 
  tableROC[1,j] <- ci(out1)[2] 
  tableCiUp[1,j] <- ci(out1)[1]
  tableCiDn[1,j] <- ci(out1)[3]
  out1 <- auc(currDataTV$mortality ~ currDataTV$lct) 
  tableROC[2,j] <- ci(out1)[2] 
  tableCiUp[2,j] <- ci(out1)[1]
  tableCiDn[2,j] <- ci(out1)[3]  
  out1 <- auc(currDataTV$mortality ~ currDataTV$wbc) 
  tableROC[3,j] <- ci(out1)[2] 
  tableCiUp[3,j] <- ci(out1)[1]
  tableCiDn[3,j] <- ci(out1)[3] 
  out1 <- auc(currDataTV$mortality ~ currDataTV$tb) 
  tableROC[4,j] <- ci(out1)[2] 
  tableCiUp[4,j] <- ci(out1)[1]
  tableCiDn[4,j] <- ci(out1)[3]
  out1 <- auc(currDataTV$mortality ~ currDataTV$pf) 
  tableROC[5,j] <- ci(out1)[2] 
  tableCiUp[5,j] <- ci(out1)[1]
  tableCiDn[5,j] <- ci(out1)[3]
  out1 <- auc(currDataTV$mortality ~ currDataTV$pt) 
  tableROC[6,j] <- ci(out1)[2] 
  tableCiUp[6,j] <- ci(out1)[1]
  tableCiDn[6,j] <- ci(out1)[3]
  out1 <- auc(currDataTV$mortality ~ currDataTV$ld) 
  tableROC[7,j] <- ci(out1)[2] 
  tableCiUp[7,j] <- ci(out1)[1]
  tableCiDn[7,j] <- ci(out1)[3]
  out1 <- auc(currDataTV$mortality ~ currDataTV$glc) 
  tableROC[8,j] <- ci(out1)[2] 
  tableCiUp[8,j] <- ci(out1)[1]
  tableCiDn[8,j] <- ci(out1)[3]
  out1 <- auc(currDataTV$mortality ~ currDataTV$cr) 
  tableROC[9,j] <- ci(out1)[2] 
  tableCiUp[9,j] <- ci(out1)[1]
  tableCiDn[9,j] <- ci(out1)[3]
  out1 <- auc(currDataTV$mortality ~ currDataTV$crp) 
  tableROC[10,j] <- ci(out1)[2] 
  tableCiUp[10,j] <- ci(out1)[1]
  tableCiDn[10,j] <- ci(out1)[3]
  } 
rownames(tableROC) <- rownames(tableCiUp) <- rownames(tableCiDn) <- 
  c("Platelet", "Lactate", "WBC", "Bilirubin", "PF ratio", "Prothrombin time", "LD", "Glucose", "Creatinine", "CRP")
colnames(tableROC) <- colnames(tableCiUp) <- colnames(tableCiDn) <- landmarkTimes
round(tableROC, 3)

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

# save.image(file = "TV bwd.RData")
