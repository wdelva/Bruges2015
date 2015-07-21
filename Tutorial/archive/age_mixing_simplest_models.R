install.packages("RSimpactCyan", repos="http://research.edm.uhasselt.be/jori")


library(RSimpactCyan)
library(EasyABC)
library(sqldf)
library(nlme)
library(lme4)
library(data.table)
library(ggplot2)



# Adding TOB to r, so that we can plot age mixing pattern
rwithIDmdata <- function(r, p){
  query <-paste0("SELECT *
                 FROM r
                 JOIN p
                 ON r.IDm = p.ID"
  )
  return (sqldf(query))
}

# Setting up the simplest model comparison
# We have two loops: Firstly we loop over scenarios, secondly we loop over n repetitions of each scenario
# For now, let's have 4+1 (+1) scenarios:
# 1. Very narrow, assortative age-mixing pattern (y~x)
# 2. Very narrow, disassortative age-mixing pattern (y~x-c)
# 3. Very narrow, age-changing age-mixing pattern (y~ax)
# 4. Wide, parallel age-mixing pattern with very little within-person variation (y_i~x_i-c_i)

# 5. Wide age-mixing pattern, due to a lot of within-person variation (y~x+e)

#### 6. "The most realistic scenario: age-changing age-mixing pattern, with widening within-person variation
#### (y_i~ax_i+(ex_i)) is not possible with the current hazard function.
#### It requires an formation.hazard.agegap.gap_factor_man and formation.hazard.agegap.gap_factor_woman
#### that are multiplied by a function of age so that the ga_factor changes with (male) age


n <- 10 # number of repetitions per scenario
nscen <-5
scenario <- 1:nscen # number of scenarios that are being compared
repetition <- 1:n
agegap.baseline <- c(6.5, 7.25, 6.55, 8.55, 4.0) #5.4)
# This is the parameter we estimated with EasyABC. 5.3 for scen1, 7.25 for scen2, 6.55 for scen3, 5.8 for scen4, 4#5.4 for scen5
agegap.gap_factor_man <- c(rep(-1,(nscen-1)), -0.1)
agegap.gap_factor_woman <- c(rep(-1,(nscen-1)), -0.1)
agegap.gap_agescale_man <- c(0, 0, 0.1, 0, 0)
agegap.gap_agescale_woman <- c(0, 0, 0.1, 0, 0)
person.agegap.man.dist.type = c("fixed", "fixed", "fixed", "normal", "fixed")
person.agegap.woman.dist.type = c("fixed", "fixed", "fixed", "normal", "fixed")
person.agegap.man.dist.fixed.value <- c(0, 10, 0, NA, 0)
person.agegap.woman.dist.fixed.value <- c(0, 10, 0, NA, 0)
person.agegap.man.dist.normal.mu <- c(NA, NA, NA, 0, NA)
person.agegap.man.dist.normal.sigma <- c(NA, NA, NA, 6.3, NA)
person.agegap.woman.dist.normal.mu <- c(NA, NA, NA, 0, NA)
person.agegap.woman.dist.normal.sigma <- c(NA, NA, NA, 6.3, NA)

for (scen in scenario){
  
  # Before we can run the model multiple times for each scenario, we want to calibrate it
  # So that roughly the same number of relationships are formed each time
  # Write an R function that has a vector of parameter values as input,
  # and a vector of summary statistics as output
  
  #setwd("/Users/wimdelva/Dropbox (Personal)/FWO-VLIR project/AgeDisConcurSurvey/Papers/Age mixing paper/R scripts")
  #source("formation_baseline_estimation.R")
  
  cfg <- list()
  #simpact.getconfig(cfg)
  cfg["population.numwomen"] <- 500
  cfg["population.nummen"] <- 500
  
  cfg["population.simtime"] <- 70#10#45
  cfg["hivseed.time"] <- 10
  cfg["hivseed.fraction"] <- 0.05
  cfg["hivtest.cd4.threshold"] <- 0
  cfg["hivtest.interval.dist.type"] <- "uniform"
  cfg["hivtest.interval.dist.uniform.min"] <- 200  # So that nobody actually gets tested and hence nobody starts ART
  cfg["hivtest.interval.dist.uniform.max"] <- 201  
  
  cfg["transmission.param.a"] <- -0.3
  cfg["transmission.param.b"] <- -8  
  cfg["transmission.param.c"] <- 0.1649
  
  
  cfg["formation.hazard.agegap.numrel_man"] <- -1
  cfg["formation.hazard.agegap.numrel_woman"] <- -1
  
  
  cfg["person.agegap.man.dist.type"] <- person.agegap.man.dist.type[scen]
  cfg["person.agegap.woman.dist.type"] <- person.agegap.woman.dist.type[scen]
  
  if (person.agegap.man.dist.type[scen] == "uniform"){
    cfg["person.agegap.man.dist.uniform.min"] <- person.agegap.man.dist.uniform.min[scen]
    cfg["person.agegap.man.dist.uniform.max"] <- person.agegap.man.dist.uniform.max[scen]    
    cfg["person.agegap.woman.dist.uniform.min"] <- person.agegap.woman.dist.uniform.min[scen]
    cfg["person.agegap.woman.dist.uniform.max"] <- person.agegap.woman.dist.uniform.max[scen]
  } else {
    cfg["person.agegap.man.dist.fixed.value"] <- person.agegap.man.dist.fixed.value[scen]
    cfg["person.agegap.woman.dist.fixed.value"] <- person.agegap.woman.dist.fixed.value[scen]    
  }
  
  cfg["formation.hazard.agegap.gap_factor_man"] <- agegap.gap_factor_man[scen]
  cfg["formation.hazard.agegap.gap_factor_woman"] <- agegap.gap_factor_woman[scen]
  cfg["formation.hazard.agegap.gap_agescale_man"] <- agegap.gap_agescale_man[scen]
  cfg["formation.hazard.agegap.gap_agescale_woman"] <- agegap.gap_agescale_woman[scen]
  
  cfg["formation.hazard.agegap.baseline"] <- agegap.baseline[scen] # This is the parameter we estimated with EasyABC
  #5.3 is good for scen1, 7.25 for scen2, 6.55 is good for scen3, 7.55 is good for scen4
  
  for (repit in repetition){
    cfg["logsystem.filename.events"] <- paste0("events",scen,"_",repit)
    cfg["logsystem.filename.persons"] <- paste0("persons",scen,"_",repit)
    cfg["logsystem.filename.relations"] <- paste0("relations",scen,"_",repit)    
    results <- simpact.run(cfg, "/Users/wimdelva/Dropbox (Personal)/FWO-VLIR project/AgeDisConcurSurvey/Papers/Age mixing paper/Modelling documents/SimOutput")  
  }
}


####
# Running the analysis on the output csv files
####

setwd("/Users/wimdelva/Dropbox (Personal)/FWO-VLIR project/AgeDisConcurSurvey/Papers/Age mixing paper/Modelling documents/SimOutput")

myDir <- "/Users/wimdelva/Dropbox (Personal)/FWO-VLIR project/AgeDisConcurSurvey/Papers/Age mixing paper/Modelling documents/SimOutput"
filenames <- list.files(myDir) 
personfilenames <- filenames[grep("personlog[.]csv", filenames)] 
relationfilenames <- filenames[grep("relationlog[.]csv", filenames)]
configfilenames <- filenames[grep("config[.]txt", filenames)]
outputfilenames <- filenames[grep("output[.]txt", filenames)]



simID <-rep(NA,5)
meanagegap <-rep(NA,5)
varagegap <-rep(NA,5)
BetweenVar <-rep(NA,5)
WithinVar <-rep(NA,5)
cumulincid <-rep(NA,5)
nrels <- rep(NA,5)
allrpdata <- data.frame()


for (i in 1:length(personfilenames)) {
  print(i)
  ptable <- fread(personfilenames[i]) #p <- read.csv(personfilenames[i], header = TRUE)    # To be improved, using fread
  rtable <- fread(relationfilenames[i])  # To be improved, using fread
  
  
  ####
  # To create heatplot, we need to loop through time steps, and calculate age- and gender-specific prevalence and incidence rates.
  alltimestepsdata <- data.table()
  timesteps <- seq(as.numeric(cfg["hivseed.time"]), as.numeric(cfg["population.simtime"]), by=1)
  for (time_i in timesteps) {
    alivetable <- alivetab(time_i)   # First we only take the data of people who were alive at time_i
    alivetable$Age.1 <- time_i -round(alivetable$TOB, 0) # Next we allocate them in discrete age bins with bin size equal as width of timesteps
    alivetable$Infected <- time_i >= alivetable$InfectTime # Now we allocate infection status to all people in our table of living people
    alive_data.table <- data.table(alivetable) # We turn the dataset into a proper data.table
    timestepdata <- alive_data.table[,sum(Infected) / nrow(alive_data.table),by="Gender,Age.1"] # And we calculate HIV prevalence, by gender and age group
    setnames(timestepdata, "V1", "Prevalence")
    timestepdata <- cbind(time_i, timestepdata)
    alltimestepsdata <- rbind(alltimestepsdata, timestepdata)
  }
  alltimestepsdata <- cbind(i, alltimestepsdata) # i is the ID number of the simulation
  
  
  ###
  # A quick try of a heat map for one run
  ###
  alltimestepsdataM <- alltimestepsdata[Gender==0,]
  p <- ggplot(alltimestepsdataM, aes(time_i, Age.1)) +
    geom_tile(aes(fill = Prevalence), colour = "white") 
  p + scale_fill_gradient(low = "white", high = "steelblue")
  
  
  rptable <- rwithIDmdata(rtable,ptable)
  
  
  
  rp <- rwithIDmdata(r,p)
  rp$id <- i
  rp$AgeMale <- -rp$TOB + rp$FormTime
  rp$AgeFemale <- rp$AgeMale - rp$AgeGap
  
  simID[i] <- i
  meanagegap[i] <- mean(r$AgeGap, na.rm=TRUE)
  varagegap[i] <- var(r$AgeGap, na.rm=TRUE)
  MixingModel <- lme(AgeFemale~AgeMale, random=~1|IDm,data=rp)
  #summary(testlme)
  VarianceEstimates <- VarCorr(MixingModel)
  BetweenVar[i] <- as.numeric(VarianceEstimates[1])
  WithinVar[i] <- as.numeric(VarianceEstimates[2])
  
  nrels[i] <- nrow(r)
  cumulincid[i] <- sum(p$InfectTime>=60 & p$InfectTime<Inf) #cfg["hivseed.time"] & p$InfectTime<Inf)
  # How many infections to place in the last 10 years of the simulation? (simulation years 50-60)
  
  allrpdata <- rbind(allrpdata, rp)
} 

summarystats <- (data.frame(simID, meanagegap, varagegap, BetweenVar, WithinVar,
                            nrels, cumulincid))
summarystats

# Let's illustrate the age mixing pattern simulated,
# superimposed on the cumulative number of cases in the last 10 years of the simulation

# Making the plot of cumulative cases, by simulation index number
par(mar=c(5, 5, 4, 2) + 0.1)
matplot(summarystats$cumulincid[c(1:50)], type = "p", pch=16, xlab = "Simulation index", ylab = "New HIV infections 50-60 after onset",
        lwd = 4, lty = 1, bty = "l", col = "black",
        cex.lab = 2,
        cex.axis= 2)
lines(c(10.5,10.5), c(0,135), lty="dashed")
lines(c(20.5,20.5), c(0,135), lty="dashed")
lines(c(30.5,30.5), c(0,135), lty="dashed")
lines(c(40.5,40.5), c(0,135), lty="dashed")
box()

# Making the plot of Within and Between Variance and Age differences, by simulation index number
par(mar=c(5, 5, 4, 2) + 0.1)
matplot(1:50, summarystats[,c(2,4,5)], type = "ppp", pch=c(1:3), cex=c(3,2,1),
        xlab = "Simulation index", ylab = "TBD",
        lwd = 4, lty = 1, bty = "l", col = c("black", "green4", "blue"),
        cex.lab = 2,
        cex.axis= 2)
lines(c(10.5,10.5), c(0,135), lty="dashed")
lines(c(20.5,20.5), c(0,135), lty="dashed")
lines(c(30.5,30.5), c(0,135), lty="dashed")
lines(c(40.5,40.5), c(0,135), lty="dashed")
box()



# Now making the 5 age mixing plots
matplot(allrpdata$AgeMale[allrpdata$id==1], allrpdata$AgeFemale[allrpdata$id==1],
        type = "p", pch=16, xlab = "Male age", ylab = "Female age",
        cex.lab = 2, cex.axis= 2,
        xlim=c(10,100), ylim=c(10,100))
lines(c(15,100), c(15,100), lty="dashed", col="orange2",lwd=8)


matplot(allrpdata$AgeMale[allrpdata$id==11], allrpdata$AgeFemale[allrpdata$id==11],
        type = "p", pch=16, xlab = "Male age", ylab = "Female age",
        cex.lab = 2, cex.axis= 2,
        xlim=c(10,100), ylim=c(10,100))
lines(c(15,100), c(15,100), lty="dashed", col="orange2",lwd=8)


matplot(allrpdata$AgeMale[allrpdata$id==21], allrpdata$AgeFemale[allrpdata$id==21],
        type = "p", pch=16, xlab = "Male age", ylab = "Female age",
        cex.lab = 2, cex.axis= 2,
        xlim=c(10,100), ylim=c(10,100))
lines(c(15,100), c(15,100), lty="dashed", col="orange2",lwd=8)


matplot(allrpdata$AgeMale[allrpdata$id==51], allrpdata$AgeFemale[allrpdata$id==51],
        type = "p", pch=16, xlab = "Male age", ylab = "Female age",
        cex.lab = 2, cex.axis= 2,
        xlim=c(10,100), ylim=c(10,100))
lines(c(15,100), c(15,100), lty="dashed", col="orange2",lwd=8)


matplot(allrpdata$AgeMale[allrpdata$id==41], allrpdata$AgeFemale[allrpdata$id==41],
        type = "p", pch=16, xlab = "Male age", ylab = "Female age",
        cex.lab = 2, cex.axis= 2,
        xlim=c(10,100), ylim=c(10,100))
lines(c(15,100), c(15,100), lty="dashed", col="orange2",lwd=8)


################
# Estimating components of variance: built into the summary stats dataset
################

MixingModel <- lme(AgeFemale~AgeMale, random=~1|IDm,data=rp)
#summary(testlme)
VarianceEstimates <- VarCorr(testlme)
BetweenVar <- as.numeric(VarianceEstimates[1])
WithinVar <- as.numeric(VarianceEstimates[2])

regressionmodel <- lm(summarystats$cumulincid ~ meanagegap*varagegap + nrels)
summary(regressionmodel)

plot(summarystats$cumulincid)

plot(rp$AgeMale, rp$AgeFemale)
lm(rp$AgeFemale ~ rp$AgeMale)


table(p$InfectTime)