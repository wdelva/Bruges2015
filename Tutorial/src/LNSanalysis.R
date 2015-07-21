#Define the project path to files
pp<-"/Users/wimdelva/Dropbox (Personal)/FWO Age Mixing/LNS/Analysis"

#General path settings
dd<-paste0(pp, "/Datasets")
dod<-paste0(pp, "/Do Files")
gd<-paste0(pp, "/Graphics")

setwd(pp)

#Create path to specific datasets
master<-paste0(dd, "/lns participants and relationships 07.01.2015.rda")
lwmastercleaned<-paste0(dd, "/lns cleaned rel (long) round (wide) 07.01.2015.rda")
llmastercleaned<-paste0(dd, "/lns cleaned rel (long) round (long) 07.01.2015.rda")
wwmastercleaned<-paste0(dd, "/lns cleaned rel (wide) round (wide) 07.01.2015.rda")
llemastercleaned<-paste0(dd,"/lns cleaned excluded rel (long) round (long) 26.01.2015.rda")
w2emastercleaned<-paste0(dd,"/lns cleaned excluded rel (wide) round (two) 26.01.2015.rda")

#######
#Open libraries
#######
library(epicalc)
library(plyr)
library(dplyr)
library(gmodels)
library(nlme)
library(lme4)
library(MCMCglmm)
library(ggplot2)
library(scales)
library(Hmisc)
library(gridExtra)
library(mgcv)


#######
#Load dataset
#######
load(llemastercleaned)
load(w2emastercleaned)

#######
#Table 1: Exclusions
#See Exclusions R script
#######

#######
#Table 2: Relationship characteristics reported by men and women
#######

#Re-categorize place and residence
clle$place2 <- revalue(clle$place, c("Mozambique" = "Elsewhere"))

clle$relresp2 <- revalue(clle$relresp, c("Mozambique" = "Elsewhere", "Chizumulu" = "Elsewhere"))

clle$curresp2 <- revalue(clle$curresp, c("Mozambique" = "Elsewhere", "Chizumulu" = "Elsewhere"))



rcatvars<-select(clle, durII, firstsexII, lastsexII,
                 otherpartpII, sfII, otherpartII, alivep, pt, ongoing, 
                 evercond, insample, jobp, place2, relresp2, hiv, hivp, testp,
                 curresp2)

rnumvars<-select(clle, agepII)

#univariate tables
lapply(rcatvars, tab1)
lapply(rnumvars, summary)

#Create function for making custom-made 2x2 tables with missing values and column proportions
tblfxn <- function (x, y) {
  count <- table(x, y, useNA = "ifany")
  freq <- round(prop.table(count, 2)*100, 1)
  
  tbl <- cbind(count[,1], freq[,1])
  for (i in 2:ncol(count)) {
    tbl <- cbind(tbl, count[,i], freq[,i])
  }
  
  cnames <- character(0)
  for (i in 1:ncol(count)) {
    cnames <- c(cnames, colnames(count)[i], "%")
  }
  colnames(tbl) <- cnames
  print(tbl)
  
  chi <- chisq.test(x, y)
  print(chi)
}

lapply(rcatvars, tblfxn, y=clle$sex)
by(rnumvars, clle$sex, summary)

#######
#Table 3: Average age difference and SD around the age differences by participant age group and gender
#######
byagecat <- group_by(clle, agecat, sex)
table3 <- summarise(byagecat,
                    n = n(),
                    mean = mean(agedif, na.rm = TRUE),
                    sd = sd(agedif, na.rm = TRUE),
                    se = sd/sqrt(n))
table3


#######
#Table 4: Average bridgewidth for people with more than one partner
#######
bysex <- group_by(subset(cwe, totalpartners > 1), sex)

summary(subset(cwe, totalpartners > 1)$bridgewidth)

table4 <- summarise(bysex,
                    n = n(),
                    median = median(bridgewidth, na.rm= TRUE),
                    iqr = IQR(bridgewidth, na.rm =TRUE))
table4

#Now calculate the fraction of people who reported more than one partner and have zero bridgewidth
cwe$bwzero <- ifelse(cwe$bridgewidth == 0, "Yes", "No")
cwe$tpone <- ifelse(cwe$totalpartners > 1, "Yes", "No")

tabpct(cwe$bwzero, cwe$tpone)

#######
#Table 5: Mean and median of participant’s mean age differences in relationships for socio-demographic characteristics (by sex)
#######

men<-subset(cwe, sex=="Male")
women<-subset(cwe, sex=="Female")

pcatvarsm<-select(men, agecat, everpolyII, elecII, floorII, ownedII, builtII, 
                  tribe, jobII, wall, roof, pit, educ, relig, marstat, test)

pcatvarsw<-select(women,agecat, everpolyII, elecII, floorII, ownedII, builtII, 
                  tribe, jobII, wall, roof, pit, educ, relig, marstat, test)

lapply(pcatvarsm, tab1)
lapply(pcatvarsw, tab1)

lapply(pcatvarsm, function(x) {
  ddply(.data = men, .variables = .(x), 
        .fun = summarize, 
        meanadmean = mean(admean))
})

lapply(pcatvarsw, function(x) {
  ddply(.data = women, .variables = .(x), 
        .fun = summarize, 
        meanadmean = mean(admean))
})

#######
#Table 6: Sexual behaviours and relationship characteristics associated with AD
#AD-both binary and coninuous
#Stratified by gender
#Mixed effects to account for hierarchical data
#######
men<-subset(clle, sex=="Male")
women<-subset(clle, sex=="Female")

predsexm<-men[,c("firstsexII", "lastsexII", "otherpartpII", "sfII", "otherpartII",
                 "pt", "ongoing", "evercond","relresp2", "curresp2")]

predsexw<-women[,c("firstsexII", "lastsexII", "otherpartpII", "sfII", "otherpartII",
                   "pt", "ongoing", "evercond","relresp2", "curresp2")]

#The true likelihood can also be approximated using numerical integration. 
#Quadrature methods are common, and perhaps most common among these use the Gaussian quadrature rule, 
#frequently with the Gauss-Hermite weighting function. It is also common to incorporate adaptive 
#algorithms that adaptively vary the step size near points with high error. 
#The accuracy increases as the number of integration points increases. 
#Using a single integration point is equivalent to the so-called Laplace approximation. 
#Each additional integration point will increase the number of computations and thus the speed 
#to convergence, although it increases the accuracy. 
#The 'nAGG' specified in the models below refers to the points of integration. 

lapply(predsexm, function(x) {
  x.model<-glmer(adcat ~ x + (1 | id), data=men, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
  print(x.model, corr=FALSE)
  se <- sqrt(diag(vcov(x.model)))
  tab <- cbind(Est = fixef(x.model), LL = fixef(x.model) - 1.96 * se, 
               UL = fixef(x.model) + 1.96 *se)
  exp(tab)  
})

#Binary female models
lapply(predsexw, function(x) {
  x.model<-glmer(adcat ~ x + (1 | id), data=women, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
  print(x.model, corr=FALSE)
  se <- sqrt(diag(vcov(x.model)))
  tab <- cbind(Est = fixef(x.model), LL = fixef(x.model) - 1.96 * se, 
               UL = fixef(x.model) + 1.96 *se)
  exp(tab)
})

#Continuous male models
lapply(predsexm, function(x) {
  x.model<-lmer(agedif ~ x + (1 | id), data=men)
  print(summary(x.model), corr=FALSE)
  confint(x.model)
})

#Continuous female models
lapply(predsexw, function(x) {
  x.model<-lmer(agedif ~ x + (1 | id), data=women)
  print(summary(x.model), corr=FALSE)
  confint(x.model)
})


#######
#Table 7: What are large bridgewidths associated with?
#Zero-inflated negative binomial regression because excess zeros and data is over-dispersed
#The SD is a lot larger than the mean
#######
library(pscl)
library(MASS)

pred<-cwe[, c("sex", "age", "agecat", "everpolyII", "elecII", "ownedII",
              "notnyanja", "pit", "educ", "relig", "marstat", "nevercond",
              "nonspouse", "onceoff", "mcp", "pmcp", "resoff")]

lapply(pred, function(x) {
  x.model<-zeroinfl(bridgewidth ~ x | x, data = cwe, dist = "negbin", EM = TRUE)
  summary(x.model)
  tab<-exp(cbind(Estimate = coef(x.model), confint(x.model)))
})

########
#Table 8:Determining the best base model for ‘average age-difference’ and ‘age of participant’ in female models 
#Using General Additive Models
########

#Create a subset of the data for women
women <- subset(cwe, sex == "Female" & !is.na(hiv.1))

#Center age on 18
women$age <- women$age-18

#Model with just linear version of age-difference
m1 <- glm(hiv.1 ~ admean, data=women, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

#Model using thin plate regression splines
m2 <- gam(hiv.1~s(admean), data=women, family=binomial)
summary(m2)
plot(m2)
AIC(m2)

#Model using linear age and splines for AD
m3 <- gam(hiv.1 ~ s(admean) + age, data=women, family=binomial)
summary(m3)
plot(m3)
AIC(m3)

#Model using splines for age and splines for AD
m4 <- gam(hiv.1 ~ s(admean) + s(age), data=women, family=binomial)
summary(m4)
plot(m4)
AIC(m4)

#Model using tensor product smooth of age and AD
m5 <- gam(hiv.1 ~ te(admean, age), data=women, family=binomial)
summary(m5)
plot(m5)
AIC(m5)

#Visualizing the best model
vis.gam(m5, type = "response", plot.type = "contour")

dev.copy(pdf, paste0(gd, '/Figure meanADfemale.pdf'))
dev.off()

########
#Table 9:Determining the best base model for ‘average age-difference’ and ‘age of participant’ in male models 
#Using General Additive Models
########

#Create a subset of the data for women
men <- subset(cwe, sex == "Male" & !is.na(hiv.1))

#Center age on 18
men$age <- men$age-18

#Model with just linear version of age-difference
m1 <- glm(hiv.1 ~ admean, data=men, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

#Model using thin plate regression splines
m2 <- gam(hiv.1~s(admean), data=men, family=binomial)
summary(m2)
plot(m2)
AIC(m2)

#Model using linear age and splines for AD
m3 <- gam(hiv.1 ~ s(admean) + age, data=men, family=binomial)
summary(m3)
plot(m3)
AIC(m3)

#Model using splines for age and splines for AD
m4 <- gam(hiv.1 ~ s(admean) + s(age), data=men, family=binomial)
summary(m4)
plot(m4)
AIC(m4)

#Model using tensor product smooth of age and AD
m5 <- gam(hiv.1 ~ te(admean, age), data=men, family=binomial)
summary(m5)
plot(m5)
AIC(m5)

#Visualizing the best model
vis.gam(m5, type = "response", plot.type = "contour")

dev.copy(pdf, paste0(gd, '/Figure meanADmale.pdf'))
dev.off()

########
#Table 10:Determining the best base model for ‘maximum age-difference’ and ‘age of participant’ in female models 
#Using General Additive Models
########

#Create a subset of the data for women
women <- subset(cwe, sex == "Female" & !is.na(hiv.1))

#Center age on 18
women$age <- women$age-18

#Model with just linear version of age-difference
m1 <- glm(hiv.1 ~ admax, data=women, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

#Model using thin plate regression splines
m2 <- gam(hiv.1~s(admax), data=women, family=binomial)
summary(m2)
plot(m2)
AIC(m2)

#Model using linear age and splines for AD
m3 <- gam(hiv.1 ~ s(admax) + age, data=women, family=binomial)
summary(m3)
plot(m3)
AIC(m3)

#Model using splines for age and splines for AD
m4 <- gam(hiv.1 ~ s(admax) + s(age), data=women, family=binomial)
summary(m4)
plot(m4)
AIC(m4)

#Model using tensor product smooth of age and AD
m5 <- gam(hiv.1 ~ te(admax, age), data=women, family=binomial)
summary(m5)
plot(m5)
AIC(m5)

#Visualizing the best model
vis.gam(m5, type = "response", plot.type = "contour")

dev.copy(pdf, paste0(gd, '/Figure maxADfemale.pdf'))
dev.off()

########
#Table 11:Determining the best base model for ‘maximum age-difference’ and ‘age of participant’ in male models 
#Using General Additive Models
########

#Create a subset of the data for women
men <- subset(cwe, sex == "Male" & !is.na(hiv.1))

#Center age on 18
men$age <- men$age-18

#Model with just linear version of age-difference
m1 <- glm(hiv.1 ~ admax, data=men, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

#Model using thin plate regression splines
m2 <- gam(hiv.1~s(admax), data=men, family=binomial)
summary(m2)
plot(m2)
AIC(m2)

#Model using linear age and linear for AD
m3 <- gam(hiv.1 ~ admax + age, data=men, family=binomial)
summary(m3)
plot(m3)
AIC(m3)

#Model using splines for age and linear for AD
m4 <- gam(hiv.1 ~ admax + s(age), data=men, family=binomial)
summary(m4)
plot(m4)
AIC(m4)

#Model using tensor product smooth of age and AD
m5 <- gam(hiv.1 ~ te(admax, age), data=men, family=binomial)
summary(m5)
plot(m5)
AIC(m5)

#Visualizing the best model
vis.gam(m4, type = "response", plot.type = "contour")

dev.copy(pdf, paste0(gd, '/Figure maxADmale.pdf'))
dev.off()

########
#Table 12:Determining the best base model for ‘bridgewidth’ and ‘age of participant’ in female models 
#Using General Additive Models
########

#Create a subset of the data for women
women <- subset(cwe, sex == "Female" & !is.na(hiv.1))

#Center age on 18
women$age <- women$age-18

#Model with just linear version of bw
m1 <- glm(hiv.1 ~ bridgewidth, data=women, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

#Model using thin plate regression splines
m2 <- gam(hiv.1~s(bridgewidth), data=women, family=binomial)
summary(m2)
plot(m2)
AIC(m2)

#Model using linear age and linear for bw
m3 <- gam(hiv.1 ~ bridgewidth + age, data=women, family=binomial)
summary(m3)
plot(m3)
AIC(m3)

#Model using splines for age and linear for bw
m4 <- gam(hiv.1 ~ bridgewidth + s(age), data=women, family=binomial)
summary(m4)
plot(m4)
AIC(m4)

#Model using tensor product smooth of age and bw
m5 <- gam(hiv.1 ~ te(bridgewidth, age), data=women, family=binomial)
summary(m5)
plot(m5)
AIC(m5)

#Visualizing the best model
vis.gam(m5, type = "response", plot.type = "contour")

dev.copy(pdf, paste0(gd, '/Figure bwfemale.pdf'))
dev.off()


########
#Table 13:Determining the best base model for ‘bridgewidth’ and ‘age of participant’ in male models 
#Using General Additive Models
########

#Create a subset of the data for women
men <- subset(cwe, sex == "Male" & !is.na(hiv.1))

#Center age on 18
men$age <- men$age-18

#Model with just linear version of bw
m1 <- glm(hiv.1 ~ bridgewidth, data=men, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

#Model using thin plate regression splines
m2 <- gam(hiv.1~s(bridgewidth), data=men, family=binomial)
summary(m2)
plot(m2)
AIC(m2)

#Model using linear age and linear for bw
m3 <- gam(hiv.1 ~ bridgewidth + age, data=men, family=binomial)
summary(m3)
plot(m3)
AIC(m3)

#Model using splines for age and linear for bw
m4 <- gam(hiv.1 ~ bridgewidth + s(age), data=men, family=binomial)
summary(m4)
plot(m4)
AIC(m4)

#Model using tensor product smooth of age and bw
m5 <- gam(hiv.1 ~ te(bridgewidth, age), data=men, family=binomial)
summary(m5)
plot(m5)
AIC(m5)

#Visualizing the best model
vis.gam(m4, type = "response", plot.type = "contour")

dev.copy(pdf, paste0(gd, '/Figure bwmale.pdf'))
dev.off()

########
#Table 14:Determining the best base model for ‘ADR in past 3 years’ and ‘age of participant’ in female models 
#Using General Additive Models
########

# Create var to indicate if participant had an ADR in the past three years.
cwe$adcatever <- as.factor(ifelse((cwe$adcat.1 == "Age-disparate" & !is.na(cwe$adcat.1)) |
                          (cwe$adcat.2 == "Age-disparate" & !is.na(cwe$adcat.2)) |
                          (cwe$adcat.3 == "Age-disparate" & !is.na(cwe$adcat.3)) | 
                          (cwe$adcat.4 == "Age-disparate" & !is.na(cwe$adcat.4)) |
                          (cwe$adcat.5 == "Age-disparate" & !is.na(cwe$adcat.5)),
                        "Yes", "No"))

#Create a subset of the data for women
women <- subset(cwe, sex == "Female" & !is.na(hiv.1))

#Center age on 18
women$age <- women$age-18

#Model with just factor version of ADR
m1 <- glm(hiv.1 ~ adcatever, data=women, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

#Model using linear age
m2 <- gam(hiv.1~adcatever + age, data=women, family=binomial)
summary(m2)
plot(m2)
AIC(m2)

#Model using splines age 
m3 <- gam(hiv.1 ~ adcatever + s(age), data=women, family=binomial)
summary(m3)
plot(m3)
AIC(m3)

#Model using splines for age and ADR interaction
m4 <- gam(hiv.1 ~ s(age, by=adcatever), data=women, family=binomial)
summary(m4)
plot(m4)
AIC(m4)

#Visualizing the best model
vis.gam(m3, type = "response", plot.type = "contour")

dev.copy(pdf, paste0(gd, '/Figure ADRfemale.pdf'))
dev.off()


########
#Table 15:Determining the best base model for ‘ADR in past 3 years’ and ‘age of participant’ in male models 
#Using General Additive Models
########

#Create a subset of the data for women
men <- subset(cwe, sex == "Male" & !is.na(hiv.1))

#Center age on 18
men$age <- men$age-18

#Model with just factor ADR
m1 <- glm(hiv.1 ~ adcatever, data=men, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

#Model using linear age
m2 <- gam(hiv.1~adcatever + age, data=men, family=binomial)
summary(m2)
plot(m2)
AIC(m2)

#Model using splines age 
m3 <- gam(hiv.1 ~ adcatever + s(age), data=men, family=binomial)
summary(m3)
plot(m3)
AIC(m3)

#Model using splines for age and ADR interaction
m4 <- gam(hiv.1 ~ s(age, by=adcatever), data=men, family=binomial)
summary(m4)
plot(m4)
AIC(m4)

#Visualizing the best model
vis.gam(m3, type = "response", plot.type = "contour")

dev.copy(pdf, paste0(gd, '/Figure ADRmale.pdf'))
dev.off()

#####################################################
#Figures
#####################################################

#Set themes for graphs

theme<-theme(axis.title.x = element_text(face="bold"),
             axis.text.x  = element_text(angle=90, face="bold", colour="black"),
             axis.title.y = element_text(face="bold", size=10),
             axis.text.y = element_text(angle=90, face="bold", colour="black"),
             panel.grid.major = element_line(colour = 'black'), 
             panel.grid.minor = element_line(colour = NA), 
             panel.background = element_rect(fill = 'white'),
             strip.background = element_rect(fill = 'white'))



#######
#Figure 1: HIV prevalence by sex and age group
#######
clle$agecat2<-cut(clle$age,
                  breaks=c(-Inf, 19, 24, 29, 34, 39, 44, 50),
                  labels=c("18-19", "20-24", "25-29", "30-24", "35-39", "40-44", "45-50"))

#Need to convert the hiv variable to a numeric variable where hiv==1, so that I can use it to calculate the prev
clle$hivgraph<-as.numeric(clle$hiv)
clle$hivgraph[clle$hivgraph == 1] <- 0
clle$hivgraph[clle$hivgraph == 2] <- 1

#Create a small dataframe that has the prevalence, upper and lower 95% CI AD by age and sex
bysexage<-group_by(clle, sex, agecat2)
hivprev <- summarise(bysexage, 
                     n=n(),
                     r = sum(hivgraph, na.rm = TRUE),
                     prev = r/n,
                     lower = prop.test(r, n)$conf.int[1],
                     upper = prop.test(r, n)$conf.int[2])


fig1<-ggplot(na.omit(hivprev), aes(x=agecat2, y=prev, group = sex, color = sex)) + 
  geom_point(stat="identity", size=3) + 
  geom_errorbar(stat="identity", aes(ymin = lower, ymax = upper), width=0.15) +
  geom_line(stat="identity", aes(color=sex)) +
  scale_y_continuous(labels = percent) +
  scale_color_brewer(name="Gender", palette="Set1") +
  xlab("Age Group") + 
  ylab("Prevalence of HIV") +  
  theme

fig1

dev.copy(pdf,paste0(gd,'/Figure 1.pdf'))
dev.off()


#######
#Figure 2: Age-gaps by age and partner type
#######
sub<-select(clle, age, agedif, pt, sex)

fig2<-ggplot(na.omit(sub), aes(x=age, y=agedif, colour=pt)) + 
  geom_point(alpha=1/4) +
  stat_smooth(aes(fill = pt), method = "gam", formula = y ~ s(x), size = 1, se=TRUE) +
  scale_fill_brewer(name="Partner Type", palette="Set1") +
  scale_color_brewer(name="Partner Type", palette="Set1") +
  xlab("Age of Participant") + 
  ylab("Partner age difference") +
  scale_y_continuous(limits=c(-20,40))+
  facet_grid(. ~ sex) +
  theme

fig2

dev.copy(pdf,paste0(gd,'/Figure 2.pdf'))
dev.off()

#######
#Figure 3: Scatter plot of ages versus partners ages, with line for x=y and fitted
#Stratified by gender 
#Build two mixed effects models (for each gender) that takes into account heteroscedastic errors
#######
library(nlme)

#Subset datasets for men so don't have to use interaction term in model
men <- subset(select(clle, id, age, agepII, sex), sex == "Male")
women <-subset(select(clle, id, age, agepII, sex), sex == "Female")

#Centre age on 18
men$age2 <- men$age - 18
women$age2 <- women$age -18

#Now fit separate models for men and women
m1m <- lme(agepII ~ age2, data = men, random = ~1 | id, method = "ML")
m1w <- lme(agepII ~ age2, data = women, random = ~1 | id, method = "ML")

#Now add the varPower variance structure for each model
#The formula to calculate the weights for variance is |v|^(2*t)
#Age can't be at 0 in varPower formula because then it will evaluate to 0 variance for the first level (18 year olds)
m2m <- update(m1m, weights = varPower(value = 0.5, form = ~age2 + 1))
m2w <- update(m1w, weights = varPower(value = 0.5, form = ~age2 + 1))

anova(m1m, m2m)
anova(m1w, m2w)

#Now fit the final models for men and women with REML instead of ML
m3m <- update(m2m, method="REML")
m3w <- update(m2w, method="REML")

#Post estimation
plot(m3m)
plot(m3w)
summary(m3m)
summary(m3w)
cim <- intervals(m3m, which = "all")
ciw <- intervals(m3w, which = "all")

#add the predicted values to the dataset
#Level 0 mean population level predictions
men$pred <- predict(m3m, men, level = 0)
women$pred <- predict(m3w, women, level = 0)

#Create design matrix
#[-2] drops the response from the formula
matm<- model.matrix(formula(m3m) [-2], men)
matw<- model.matrix(formula(m3w) [-2], women)

predvarm <- diag(matm %*% vcov(m3m) %*% t(matm))
predvarw <- diag(matw %*% vcov(m3w) %*% t(matw))

men$se <-sqrt(predvarm)
women$se <- sqrt(predvarw)

#Need to calculate a vector of variances because the variances should be different for each age
#Create a vector of the variance weights. This is based upon the formula for varPower |v|^(2*t)
#(2*t) is the power
#|v| is the absolute value of each of the ages from the age vector
powerm <- (attributes(m3m$apVar)$Pars["varStruct.power"])
powerw <- (attributes(m3w$apVar)$Pars["varStruct.power"])
men$weights <- (1 + men$age2)^ powerm
women$weights <- (1 + women$age2)^powerw

#Extracting the residual variance for 18 year olds (0— baseline)
basevarm <- m3m$sigma^2
basevarw <- m3w$sigma^2

#Now create a vector that contains the residual variances for all ages
men$resvar <- basevarm * men$weights
women$resvar <- basevarw * women$weights

#Extract the variance of the random intercept
men$interceptvar <- getVarCov(m3m)[1,1]
women$interceptvar <- getVarCov(m3w)[1,1]

#Caluclate the prediction stardard errors
men$se2 <- sqrt(men$interceptvar + men$resvar)
women$se2 <- sqrt(women$interceptvar + women$resvar)

#Now plot
fig3a<-ggplot(men, aes(x = age, y = pred)) +
  geom_point(aes(x=age, y=agepII), position=position_jitter(width=0.5,height=0.5), alpha = 0.5) +
  geom_ribbon(aes(ymin=pred-2*se, ymax=pred+2*se, fill = "95% Confidence Interval"), alpha=0.5) +
  geom_ribbon(aes(ymin=pred-2*se2, ymax=pred+2*se2, fill = "95% Prediction Interval"), alpha=0.2) +
  geom_abline(intercept=0, slope=1, size = 1, aes(color = "One-to-one")) +
  geom_line(aes(color="Best fit"), size = 1) +
  scale_y_continuous(name="Partners age") +
  scale_fill_manual('Intervals', values = c('black', 'red')) +
  scale_color_manual('Lines', values = c('black', 'blue')) +
  ggtitle("Men") +
  xlab("Participant's age") +
  annotate("text", x = 28, y = 48, label = "Power is 0.24, (95% CI: 0.19-0.30)", size = 3)# +
  #theme

fig3a

#Now plot for the age-mixing tutorial in Bruges
fig3Bruges<-ggplot(men, aes(x = age, y = pred)) +
  geom_point(aes(x=age, y=agepII), position=position_jitter(width=0.5,height=0.5), alpha = 0.5) +
  #geom_ribbon(aes(ymin=pred-2*se, ymax=pred+2*se, fill = "95% Confidence Interval"), alpha=0.5) +
  geom_ribbon(aes(ymin=pred-2*se2, ymax=pred+2*se2, fill = "95% Prediction Interval"), alpha=0.2) +
  geom_abline(intercept=0, slope=1, size = 1, aes(color = "Diagonal (x = y)")) +
  geom_line(aes(color="Linear fit"), size = 1) +
  scale_y_continuous(name="Age woman") +
  scale_fill_manual('', values = c('black', 'red')) +
  scale_color_manual('', values = c('black', 'blue')) +
  #ggtitle("Men") +
  xlab("Age man") +
  #annotate("text", x = 28, y = 48, label = "Power is 0.24, (95% CI: 0.19-0.30)", size = 3) +
  theme_bw(base_size = 18, base_family = "")
#theme

fig3Bruges


dev.copy(pdf,paste0(gd,'/Figure3Bruges.pdf'), width=9.5)
dev.off()

fig3b<-ggplot(women, aes(x = age, y = pred)) +
  geom_point(aes(x=age, y=agepII), position=position_jitter(width=0.5,height=0.5), alpha = 0.5) +
  geom_ribbon(aes(ymin=pred-2*se, ymax=pred+2*se, fill = "95% Confidence Interval"), alpha=0.5) +
  geom_ribbon(aes(ymin=pred-2*se2, ymax=pred+2*se2, fill = "95% Prediction Interval"), alpha =0.2) +
  geom_abline(intercept=0, slope=1, size = 1, aes(color = "One-to-one")) +
  geom_line(aes(color="Best fit"), size = 1) +
  scale_y_continuous(name="Partners age") +
  scale_color_manual('Lines', values = c('black', 'blue')) +
  scale_fill_manual('Intervals', values = c('black', 'red')) +
  xlab("Participant's age") +
  ggtitle("Women") +
  annotate("text", x = 28, y = 70, label = "Power is 0.19, (95% CI:0.14-0.24)", size = 3) +
  theme

fig3b

dev.copy(pdf,paste0(gd,'/Figure 3b.pdf'))
dev.off()





##########################################################################################

#######
#Tables : Are averaage age gaps, maximum age gaps, bridgewidths, and having had an age-disparate relationship related to the participants' HIV status? 
#Stratify by gender of participant
#Univariate
#######

# Create var to indicate if participant had an ADR in the past three years.
cwe$adcatever <- ifelse((cwe$adcat.1 == "Age-disparate" & !is.na(cwe$adcat.1)) |
                          (cwe$adcat.2 == "Age-disparate" & !is.na(cwe$adcat.2)) |
                          (cwe$adcat.3 == "Age-disparate" & !is.na(cwe$adcat.3)) | 
                          (cwe$adcat.4 == "Age-disparate" & !is.na(cwe$adcat.4)) |
                          (cwe$adcat.5 == "Age-disparate" & !is.na(cwe$adcat.5)),
                        "Yes", "No")

#Create subsets for men and women
women <- subset(cwe, sex == "Female" & !is.na(hiv.1))
men <- subset (cwe, sex == "Male" & !is.na(hiv.1))

#Center age on 18
women$age <- women$age-18
men$age <- men$age-18

#Models for mean agedif, women
m1<-glm(hiv.1 ~ admean, data=women, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

m2 <- glm(hiv.1 ~ admean + age + totalpartners, data = women, family = binomial)
summary(m2)
exp(cbind(OR = coef(m2), confint(m2)))

m3 <- glm(hiv.1 ~ admean + age + totalpartners + admean*age, data = women, family = binomial)
summary(m3)
exp(cbind(OR = coef(m3), confint(m3)))

#Model for max agedif, women
m1<-glm(hiv.1 ~ admax, data=women, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

m2 <- glm(hiv.1 ~ admax + age + totalpartners, data = women, family = binomial)
summary(m2)
exp(cbind(OR = coef(m2), confint(m2)))

m3 <- glm(hiv.1 ~ admax + age + totalpartners + admax*age, data = women, family = binomial)
summary(m3)
exp(cbind(OR = coef(m3), confint(m3)))

#Model for AD,  women
m1<-glm(hiv.1 ~ adcatever, data=women, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

m2 <- glm(hiv.1 ~ adcatever + age + totalpartners, data = women, family = binomial)
summary(m2)
exp(cbind(OR = coef(m2), confint(m2)))

m3 <- glm(hiv.1 ~ adcatever + age + totalpartners + adcatever*age, data = women, family = binomial)
summary(m3)
exp(cbind(OR = coef(m3), confint(m3)))

#Model for bridgewidths, women
m1<-glm(hiv.1 ~ bridgewidth, data=women, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

m2 <- glm(hiv.1 ~ bridgewidth + age + totalpartners, data = women, family = binomial)
summary(m2)
exp(cbind(OR = coef(m2), confint(m2)))

m3 <- glm(hiv.1 ~ bridgewidth + age + totalpartners + bridgewidth*age, data = women, family = binomial)
summary(m3)
exp(cbind(OR = coef(m3), confint(m3)))

#Models for mean agedif, men
m1<-glm(hiv.1 ~ admean, data=men, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

m2 <- glm(hiv.1 ~ admean + age + totalpartners, data = men, family = binomial)
summary(m2)
exp(cbind(OR = coef(m2), confint(m2)))

m3 <- glm(hiv.1 ~ admean + age + totalpartners + admean*age, data = men, family = binomial)
summary(m3)
exp(cbind(OR = coef(m3), confint(m3)))

#Models for max agedif, men
m1<-glm(hiv.1 ~ admax, data=men, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

m2 <- glm(hiv.1 ~ admax + age + totalpartners, data = men, family = binomial)
summary(m2)
exp(cbind(OR = coef(m2), confint(m2)))

m3 <- glm(hiv.1 ~ admax + age + totalpartners + admax*age, data = men, family = binomial)
summary(m3)
exp(cbind(OR = coef(m3), confint(m3)))


#Model for AD,  men
m1<-glm(hiv.1 ~ adcatever, data=men, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

m2 <- glm(hiv.1 ~ adcatever + age + totalpartners, data = men, family = binomial)
summary(m2)
exp(cbind(OR = coef(m2), confint(m2)))

m3 <- glm(hiv.1 ~ adcatever + age + totalpartners + adcatever*age, data = men, family = binomial)
summary(m3)
exp(cbind(OR = coef(m3), confint(m3)))

#Model for bridgewidths, men
m1<-glm(hiv.1 ~ bridgewidth, data=men, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

m2 <- glm(hiv.1 ~ bridgewidth + age + totalpartners, data = men, family = binomial)
summary(m2)
exp(cbind(OR = coef(m2), confint(m2)))

m3 <- glm(hiv.1 ~ bridgewidth + age + totalpartners + bridgewidth*age, data = men, family = binomial)
summary(m3)
exp(cbind(OR = coef(m3), confint(m3)))

########
#Table 12: Is participant's HIV status related to bridgewidth and average age difference?
#Stratified by gender
#######

#Create subsets for men and women
women <- subset(cwe, sex == "Female" & !is.na(hiv.1))
men <- subset (cwe, sex == "Male" & !is.na(hiv.1))

#Center age on 18
women$age <- women$age-18
men$age <- men$age-18

#For women
m1 <- glm(hiv.1 ~ admean + bridgewidth + age + totalpartners, data = women, family = binomial)
summary(m1)
exp(cbind(OR = coef(m1), confint(m1)))

#For men
m2 <- glm(hiv.1 ~ admean + bridgewidth + age + totalpartners, data = men, family = binomial)
summary(m2)
exp(cbind(OR = coef(m2), confint(m2)))




rm(list=ls())
detach(package:lme4)
detach(package:MASS)
detach(package:pscl)
detach(package:nlme)
detach(package:gmodels)
detach(package:epicalc)
detach(package:dplyr)
detach(package:mgcv)




