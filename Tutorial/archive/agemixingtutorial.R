#########
# Age-mixing tutorial for David and Roxy
#########


# GETTING STARTED

# Loading the RSimpactCyan package and other packages that will use in the
# post-simulation analysis
library(RSimpactCyan)
library(EasyABC)
library(data.table)
library(nlme)
library(ggplot2)
library(RColorBrewer)
library(shape)




# In this tutorial, we will set up a simple age-mixing simulation study, and we
# will calibrate the model to a few target statistics. We will simulate 5
# different age-mixing patterns.



# RUNNING THE DEFAULT MODEL

# But first, let's run the default model that is created in Simpact. This is
# done with the function simpact.run The function needs as arguments a list with
# the configuration of the model parameters. To use the default parameters, you
# can specify an empty list.
cfg <- list()
# You also need to say in which directory the results of the simulation must be
# stored Here we specify a temporary directory that doesn't exist yet. The
# programme will create this directory on the fly. But it will be erase when we
# open a new R session.
DestDir <- "/tmp/testdirectory"
# An example of directory that would not be deleted after you close the R
# session: DestDir <- "/Users/wimdelva/Dropbox (Personal)/FWO Age
# Mixing/Papers/Simplest age mixing simulation/configfiles" Lastly, we need to
# specify the identifier of the output files that we will create and store. Here
# we are only running 1 simulation. So we could call it today's date and attach
# the number 1 to it:
simID <- 1
identifier <- paste0("%T-%y-%m-%d-", simID)

testrun <- simpact.run(cfg,
                       DestDir,
                       identifierFormat = identifier)


# To access the output:

# Make sure you insert today's date
personlogfilename <- paste0(DestDir, "/simpact-cyan-2015-02-25-", simID, "personlog.csv")
relationlogfilename <- paste0(DestDir, "/simpact-cyan-2015-02-25-", simID, "relationlog.csv")
eventlogfilename <- paste0(DestDir, "/simpact-cyan-2015-02-25-", simID, "eventlog.csv")

ptable <- fread(personlogfilename, sep = ",", skip = 0)
rtable <- fread(relationlogfilename, sep = ",", skip = 0)
etable <- fread(eventlogfilename, sep = ",", skip = 0)
# The variables of the eventlog file haven't been named yet, so let's do that now
setnames(etable, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"),
         c("eventtime", "eventname", "p1name", "p1ID", "p1gender", "p1age", "p2name", "p2ID", "p2gender", "p2age"))

# and to view the first 6 rows of these tables:
View(head(ptable)) # Here each row is a person
View(head(rtable)) # Here each row is a relationship
View(head(etable)) # Here each row is an event

nrow(ptable) # In total, 229 people existed at some point in the simulation
nrow(rtable) # These 229 people jointly had 1004 relationships
             # over the course of the simulation (15 years)
tail(etable, 2) # The last event took place at time 15.00298
# Since this last event took place after 15 years, the criterium for stopping the 
# simulation was reached, and no more events took place thereafter


# MODIFYING MODEL PARAMETERS

# To see the parameters that can be modified, what options you have when
# modifying them, and what the default values are:
simpact.showconfig(NULL)
# As you can see, there are many parameters that you can change
# and not all parameters are used.
# For instance, there are two types of hazard function for the event
# of relationship formation that can be specified:
# formation.hazard.type = "simple" OR formation.hazard.type = "agegap"
# And the default hazard function is "agegap".

# Let's say we wanted to run the model for a larger population (500 men and 500
# women) but for a shorter period of time (3 years) We can do this by making
# changes to the cfg list (which was empty thusfar).
cfg["population.numwomen"] <- 500
cfg["population.nummen"] <- 500
cfg["population.simtime"] <- 3

# Let's run the simulation again
simID <- 2
identifier <- paste0("%T-%y-%m-%d-", simID)

testrun <- simpact.run(cfg,
                       DestDir,
                       identifierFormat = identifier)


# To access the output:

# Make sure you insert today's date
personlogfilename <- paste0(DestDir, "/simpact-cyan-2015-02-25-", simID, "personlog.csv")
relationlogfilename <- paste0(DestDir, "/simpact-cyan-2015-02-25-", simID, "relationlog.csv")
eventlogfilename <- paste0(DestDir, "/simpact-cyan-2015-02-25-", simID, "eventlog.csv")

ptable <- fread(personlogfilename, sep = ",", skip = 0)
rtable <- fread(relationlogfilename, sep = ",", skip = 0)
etable <- fread(eventlogfilename, sep = ",", skip = 0)
# The variables of the eventlog file haven't been named yet, so let's do that now
setnames(etable, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"),
         c("eventtime", "eventname", "p1name", "p1ID", "p1gender", "p1age", "p2name", "p2ID", "p2gender", "p2age"))

nrow(ptable) # In total, 1016 people existed at some point in the simulation
nrow(rtable) # These 1016 people jointly had 826 relationships
# over the course of the simulation (3 years)
tail(etable, 2) # The last event took place at time 3.0014

# Let's look at the age mixing pattern now Scatter plot of male and female age
# of couples at the time of relationship formation In order to do this, we need
# to merge the personlog and relationlog tables
DTRandDTP_a <- merge(data.frame(rtable), data.frame(ptable), by.x = "IDm", by.y = "ID")
DTRandDTP_b <- merge(DTRandDTP_a, ptable, by.x = "IDw", by.y = "ID", suffixes = c(".m", ".w"))
names(DTRandDTP_b)

# TOB is the time of birth.
# If you were 20 years old at the start of the simulation (time=0)
# then your TOB is -20

matplot(- DTRandDTP_b$TOB.m + DTRandDTP_b$FormTime,
        - DTRandDTP_b$TOB.w + DTRandDTP_b$FormTime,
        type = "p", pch=16, xlab = "Male age", ylab = "Female age",
        cex.lab = 1, cex.axis= 1,
        xlim=c(10,100), ylim=c(10,100))
lines(c(15,100), c(15,100), lty="dashed", col="orange2",lwd=8)


# FITTING THE MODEL TO DATA (SUMMARY STATISTICS) Let's say we wanted to model an
# age-mixing pattern with a mean age gap of 3 years (i.e. the man is 3 years
# older than the women) and a standard deviation of 2 years, so that 95% of the
# age gaps would fall in the range - 1 to 7 years age gap.

# Furthermore, let's recall that in our second simulation, 1016 people jointly
# had 826 relationships over the course of the simulation (3 years) This means
# that the average person had 0.54 relationships per year 826 / (1016/2) / 3 
# Note that each relationship connects 2 people, hence 1016/2 Suppose we wanted
# to simulate more intensive relationship formation, e.g. an average of 1
# relationship per person per year

# We have 3 target statistics:
# mean age gap = 3 years
# standard deviation in age gaps = 2 years
# average number of relationships per person per year = 1

# Making sure that the mean age gap is about 3 years,is straightforward,
# since this can be achieved by modifying the cfg list:
cfg["person.agegap.man.dist.type"] <- "fixed"
cfg["person.agegap.woman.dist.type"] <- "fixed"
cfg["person.agegap.man.dist.fixed.value"] <- 3
cfg["person.agegap.woman.dist.fixed.value"] <- 3

# Controlling the standard deviation in age gap and the average number of new
# relationships per person per year, is less straightforward. From the
# relationship formation hazard function, we know that we can use the "gap
# factors" to tune how quickly the hazard drops as candidate couples have age
# gaps that deviate from 3 years. But we don't know what the best value for the
# gap factors is to achieve a standard deviation of 2 years. Similarly, we know
# that we can use the "baseline" parameter to tune the rate of relationship
# formation, and we know that we need a value that is higher than 0.1 but we
# don't know how much higher.

# One way is to experiment by trial and error, but this is time-consuming, 
# difficult to reproduce, and it is unclear when you can say your guess is close
# enough.

# A more systematic way to go about this is to use "Approximate Bayesian
# Computation" to find model parameters that lead to simulation output that is 
# very close to the predefined target statistics. In R we can use the EasyABC
# package for this.

# We first create a function that runs the model, based on an "inputvector", and
# returns on "outputvector" The input vector will be the vector of parameter
# values that is sampled from the distribution of possible/likely parameter
# values In the first iteration of the ABC algorithm, this distribution is the 
# prior distribution that we specified, but in later iterations, it is the 
# updated distribution, based on how good or how bad the output from the first
# iteration approximated the target statistics

ABC_DestDir <- "/tmp/ABC/"
ABC_identifier <- "ABC"

simpact4ABC <- function(inputvector){
  cfg <- list()
  cfg["population.numwomen"] <- 500
  cfg["population.nummen"] <- 500
  cfg["population.simtime"] <- 3
  cfg["person.agegap.man.dist.type"] <- "fixed"
  cfg["person.agegap.woman.dist.type"] <- "fixed"
  cfg["person.agegap.man.dist.fixed.value"] <- 3
  cfg["person.agegap.woman.dist.fixed.value"] <- 3
  cfg["formation.hazard.agegap.baseline"] <- inputvector[1]
  cfg["formation.hazard.agegap.gap_factor_man"] <- inputvector[2]
  cfg["formation.hazard.agegap.gap_factor_woman"] <- inputvector[2]
  results <- simpact.run(cfg, "/tmp/testdirectory/") #simpact.run(cfg, ABC_DestDir, identifierFormat = ABC_identifier)  
  # relationlogfilename <- paste0(ABC_DestDir, ABC_identifier, "relationlog.csv")
  # personlogfilename <- paste0(ABC_DestDir, ABC_identifier, "personlog.csv")
  ABC_rtable <- fread(results["logrelations"], sep = ",", skip = 0) # fread(relationlogfilename, sep = ",", skip = 0)
  ABC_ptable <- read.csv(results["logpersons"], sep = ",", skip = 0) # fread(personlogfilename, sep = ",", skip = 0)
  relsperpersonperyear <- nrow(ABC_rtable) / (nrow(ABC_ptable)/2) / 3
  agegapsd <- sd(ABC_rtable$AgeGap)
  outputvector <- c(relsperpersonperyear, agegapsd)
  return(outputvector)
}


# We specify the prior distributions for the input parameters The relationship
# formation rate must roughly double:
# hF' = 2hF hF = exp(a0 + a1X1 + a2X2 + .... + anXn)
# 2hF = exp(a0 + a1X1 + a2X2 + .... + anXn) * 2
# 2hF = exp(a0 + a1X1 + a2X2 + .... + anXn) * exp(log(2))
# 2hF = exp(a0 + log(2) + a1X1 + a2X2 + .... + anXn)
# So we would naively expect that the baseline parameter (0.1) should be
# increased by log(2) ~ 0.7 to 0.8 However, we are also adjusting the "gap
# factors" and making relationships with large age gaps less likely will result
# in an overall decrease in the number of relationships formed per time unit.

simpact_prior <- list(c("unif", 0.8, 5), c("unif", -1, 0))


# Lastly, we specify the target summary statistic
sum_stat_obs <- c(1, 2)

# Here the numbers 1 and 2 are the target statistics, namely that we want to
# have a model that produces an average of 1 new relationship per person per
# year, and a standard deviation of age gaps of 2.sum_stat_obs <- c(1, 2)

# Now we try to run a sequential ABC scheme, according to the method proposed by
# Lenormand et al. 2013 Maxime Lenormand, Franck Jabot and Guillaume Deffuant.
# Adaptive approximate Bayesian computation for complex models. Comput Stat
# (2013) 28:2777–2796 DOI 10.1007/s00180-013-0428-3

# Initial number of simulations
n_init <- 40
alpha <- 0.5 # This is the fraction of n_init that must be retained by the end of the procedure
pacc <- 0.2 #0.02 # This is the criterion for how close the output must approximate the target statistics
# pacc must be between 0 and 1. The smaller, the more strict the criterion.

ABC_LenormandResult <- ABC_sequential(method="Lenormand",
                                      model=simpact4ABC,
                                      prior=simpact_prior,
                                      nb_simul=n_init,
                                      summary_stat_target=sum_stat_obs,
                                      alpha=alpha,
                                      p_acc_min=pacc,
                                      verbose=FALSE)

ABC_LenormandResult

hist(ABC_LenormandResult$param[,1])
hist(ABC_LenormandResult$param[,2])


# Looking at these histograms, it seems like 3.65 and -0.35 are good values for
# the baseline and gap factor parameters, respectively. Let’s update our model
# and have a look at the age-mixing pattern again.

simID <- 3
identifier <- paste0("%T-%y-%m-%d-", simID)

# Just to be complete, we reconstruct the cfg list from scratch
cfg <- list()
cfg["population.numwomen"] <- 500
cfg["population.nummen"] <- 500
cfg["population.simtime"] <- 3
cfg["person.agegap.man.dist.type"] <- "fixed"
cfg["person.agegap.woman.dist.type"] <- "fixed"
cfg["person.agegap.man.dist.fixed.value"] <- 3
cfg["person.agegap.woman.dist.fixed.value"] <- 3
cfg["formation.hazard.agegap.baseline"] <- 3.65
cfg["formation.hazard.agegap.gap_factor_man"] <- -0.35
cfg["formation.hazard.agegap.gap_factor_woman"] <- -0.35



testrun <- simpact.run(cfg,
                       DestDir,
                       identifierFormat = identifier)


# To access the output:

# Make sure you insert today's date
personlogfilename <- paste0(DestDir, "/simpact-cyan-2015-02-26-", simID, "personlog.csv")
relationlogfilename <- paste0(DestDir, "/simpact-cyan-2015-02-26-", simID, "relationlog.csv")
eventlogfilename <- paste0(DestDir, "/simpact-cyan-2015-02-26-", simID, "eventlog.csv")

ptable <- fread(personlogfilename, sep = ",", skip = 0)
rtable <- fread(relationlogfilename, sep = ",", skip = 0)
etable <- fread(eventlogfilename, sep = ",", skip = 0)
# The variables of the eventlog file haven't been named yet, so let's do that now
setnames(etable, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"),
         c("eventtime", "eventname", "p1name", "p1ID", "p1gender", "p1age", "p2name", "p2ID", "p2gender", "p2age"))

nrow(ptable) # In total, 1032 people existed at some point in the simulation
nrow(rtable) # These 1032 people jointly had 1718 relationships
# over the course of the simulation (3 years), that is an average of 1.1
# relationships per person per year
tail(etable, 2) # The last event took place at time 3.0012
mean(rtable$AgeGap) # Average age gap is 2.9
sd(rtable$AgeGap) # And the standard deviation of age gaps is now 1.97

# Let's look at the age mixing pattern now Scatter plot of male and female age
# of couples at the time of relationship formation In order to do this, we need
# to merge the personlog and relationlog tables
DTRandDTP_a <- merge(data.frame(rtable), data.frame(ptable), by.x = "IDm", by.y = "ID")
DTRandDTP_b <- merge(DTRandDTP_a, ptable, by.x = "IDw", by.y = "ID", suffixes = c(".m", ".w"))

# TOB is the time of birth.
# If you were 20 years old at the start of the simulation (time=0)
# then your TOB is -20

matplot(- DTRandDTP_b$TOB.m + DTRandDTP_b$FormTime,
        - DTRandDTP_b$TOB.w + DTRandDTP_b$FormTime,
        type = "p", pch=16, xlab = "Male age", ylab = "Female age",
        cex.lab = 1, cex.axis= 1,
        xlim=c(10,100), ylim=c(10,100))
lines(c(15,100), c(15,100), lty="dashed", col="orange2",lwd=8)





# AGE-MIXING AND HIV TRANSMISSION

# Let's now look at the HIV epidemic that evolves in this population.
# Of particular interest is how HIV prevalence changes over time, for different age and gender groups.
# We will introduce HIV into the population by infecting 1% of the population
# after 10 years of simulation time, so that the population has had time to reach a dynamic equilibrium
# with respect to the sexual network configuration.
# We will further constrain the initial HIV "seed infections" to the age group 20-25 year olds, 
# so that we can easily see how the epidemic spreads from this birth cohort to other birth cohorts.
# We will run the model for 40 years, i.e. 30 years of transmission after HIV is introduced.


simID <- 4
identifier <- paste0("%T-%y-%m-%d-", simID)

# Just to be complete, we reconstruct the cfg list from scratch
cfg <- list()
cfg["population.numwomen"] <- 500
cfg["population.nummen"] <- 500
cfg["population.simtime"] <- 3
cfg["person.agegap.man.dist.type"] <- "fixed"
cfg["person.agegap.woman.dist.type"] <- "fixed"
cfg["person.agegap.man.dist.fixed.value"] <- 3
cfg["person.agegap.woman.dist.fixed.value"] <- 3
cfg["formation.hazard.agegap.baseline"] <- 3.65
cfg["formation.hazard.agegap.gap_factor_man"] <- -0.35
cfg["formation.hazard.agegap.gap_factor_woman"] <- -0.35

cfg["population.simtime"] <- 40 # 40 is 10 years of "burn in", followed by 30 years of HIV transmission
cfg["hivseed.time"] <- 10
# We want to have precise control over the number (not just the fraction) of HIV
# seeds being introduced
cfg["hivseed.type"] <- "amount"
cfg["hivseed.amount"] <- 10     # 10 out of 1000 is 1% of the population
# HIV is introduced at an age that the individual is likely already sexually
# active, in a relationship, and still has lots of time to infect others
cfg["hivseed.age.min"] <- 20
cfg["hivseed.age.max"] <- 25

# Nobody will start ART during the simulations
cfg["monitoring.cd4.threshold"] <- 0.001


# We are inflating the transmission parameters "for effect"
cfg["transmission.param.a"] <- -0.3
cfg["transmission.param.b"] <- -8  
cfg["transmission.param.c"] <- 0.1649

testrun <- simpact.run(cfg,
                       DestDir,
                       identifierFormat = identifier)

# We've written a function (readthedata) to read in the csv files, to make this process leaner and cleaner.
outputID <- testrun["id"]
datalist <- readthedata(DestDir, outputID) # a list that contains ptable, rtable and etable

# We've also written functions to calculate HIV prevalence and HIV incidence by age group, gender and time point,
# and functions to look at the number of people that each infected individual infects (R),
# the age at which people transmit HIV, and the age mixing pattern.
# To load these functions, we can run the R script "postsim.R"
# source("/path/to/postsim.R")
source("/Users/wimdelva/Dropbox (Personal)/FWO Age Mixing/Papers/Simplest age mixing simulation/agemixr/R/postsim.R")

allPrevalencetimestepsdata <- prevalenceheatmapdata(datalist$ptable, cfg)
allIncidencetimestepsdata <- incidenceheatmapdata(datalist$ptable, datalist$etable, cfg)
#Incidencetimestepdata <- allIncidencetimestepsdata[ ,IncidenceTime_i := sum(eventname.incident=="transmission", na.rm = TRUE) / sum(PY),by = "time_i"]
Rtimestepsdata <- allRtimestepsdata(datalist$ptable, datalist$etable, cfg)


ageandtimeattransmissiondata <- transmissionmap(datalist$ptable, datalist$etable, cfg)
agemixingdata <- agemixing(datalist$ptable, datalist$rtable)


# Now let's visualise the HIV prevalence model output
my.pal <- c("#FFFFFF", brewer.pal(9,"YlOrRd"))


# allPrevalencetimestepsdataM <- allPrevalencetimestepsdata[Gender==0, ]
# p <- ggplot(allPrevalencetimestepsdataM, aes(time_i, Age)) +
#   geom_tile(aes(fill = Prevalence), colour = "white") 
# p + scale_fill_gradient(low = "white", high = "steelblue")
# 
# # And now we look at HIV prevalence over time and by age, in women:
# allPrevalencetimestepsdataW <- allPrevalencetimestepsdata[Gender==1, ]
# p <- ggplot(allPrevalencetimestepsdataW, aes(time_i, Age)) +
#   geom_tile(aes(fill = Prevalence), colour = "white") 
# p + scale_fill_gradient(low = "white", high = "darkred")

# We can create a 2-panel plot to look at HIV prevalence in men and women
allPrevalencetimestepsdata$Gender <- as.factor(allPrevalencetimestepsdata$Gender)
levels(allPrevalencetimestepsdata$Gender) <- c("Men", "Women")

p <- ggplot(allPrevalencetimestepsdata, aes(time_i, Age, Prevalence)) +
  geom_tile(aes(fill = Prevalence))
p + scale_fill_gradientn(colours = my.pal) +
  facet_wrap(~Gender, nrow=1) +
  theme(panel.grid = element_blank()) +
  xlab("Simulation time")

# Now let's visualise the HIV incidence model output
allIncidencetimestepsdata$Gender <- as.factor(allIncidencetimestepsdata$Gender)
levels(allIncidencetimestepsdata$Gender) <- c("Men", "Women")
# We don't want weirdly high incidence values for strata where the number of PY
# exposure time was very small

p <- ggplot(allIncidencetimestepsdata[allIncidencetimestepsdata$PY_G_A>=5, ], aes(time_i, Age, Incidence_G_A)) +
  geom_tile(aes(fill = Incidence_G_A))
p + scale_fill_gradientn(colours = my.pal) +
  facet_wrap(~Gender, nrow=1) +
  theme(panel.grid = element_blank()) +
  xlab("Simulation time")


# But we are also interested in the population-average HIV incidence over time
names(Incidencetimestepdata)
p <- ggplot(allIncidencetimestepsdata, aes(x=time_i, y=Incidence_G, group=1)) +
  geom_line() +
  xlab("Simulation time") +
  ylab(" HIV incidence") +
  facet_wrap(~Gender, nrow=1)
print(p)
  
  #aes(fill = Incidence))
p + scale_fill_gradientn(colours = my.pal) +
  facet_wrap(~Gender, nrow=1) +
  theme(panel.grid = element_blank()) +
  xlab("Simulation time")


# VISUALISING THE HIV TRANSMISSION NETWORK
# Since we keep track of all the events, we can visualise when exactly HIV was
# transmitted and how old the infectors and newly infecteds were at the time of
# transmission.

agetrans <- ageandtimeattransmissiondata$transmissiondata

plot(c(0, 45), c(0, 75), type="n",
     xlab="Simulation time",
     ylab="Ages of transmission pairs")

Arrows(agetrans$timeattransm,
       agetrans$ageoftransmitter,
       agetrans$timeattransm,
       agetrans$ageofreceptor,
       arr.length = 0.2, code = 2,
       arr.type = "triangle", col = - agetrans$genderoftransmitter + 3)
legend("topleft",
       legend = c("male-to-female", "female-to-male"),
       fill = c(3,2),
       border = c(3,2))


# CHANGING THE AGE_MIXING PATTERN "MID-COURSE"
# We will now introduce an "interventions" in the middle of our simulation.
# Specifically, we will make the gap factor more negative, such that the
# age-mixing scatter becomes narrower, and fewer age-disparate relationships are
# formed.



simID <- 6
identifier <- paste0("%T-%y-%m-%d-", simID)

# Just to be complete, we reconstruct the cfg list from scratch
cfg <- list()
cfg["population.numwomen"] <- 1000
cfg["population.nummen"] <- 1000
cfg["person.agegap.man.dist.type"] <- "fixed"
cfg["person.agegap.woman.dist.type"] <- "fixed"
cfg["person.agegap.man.dist.fixed.value"] <- 3
cfg["person.agegap.woman.dist.fixed.value"] <- 3
cfg["formation.hazard.agegap.baseline"] <- 3.65
cfg["formation.hazard.agegap.gap_factor_man"] <- -0.35
cfg["formation.hazard.agegap.gap_factor_woman"] <- -0.35

cfg["population.simtime"] <- 70 # 40 is 10 years of "burn in", followed by 30 years of HIV transmission
cfg["hivseed.time"] <- 10
# We want to have precise control over the number (not just the fraction) of HIV
# seeds being introduced
cfg["hivseed.type"] <- "amount"
cfg["hivseed.amount"] <- 10     # 10 out of 1000 is 1% of the population
# HIV is introduced at an age that the individual is likely already sexually
# active, in a relationship, and still has lots of time to infect others
cfg["hivseed.age.min"] <- 20
cfg["hivseed.age.max"] <- 30

# Nobody will start ART during the simulations
cfg["monitoring.cd4.threshold"] <- 0.001


# We are inflating the transmission parameters "for effect"
cfg["transmission.param.a"] <- -0.3
cfg["transmission.param.b"] <- -8  
cfg["transmission.param.c"] <- 0.1649

gapfactorchange <- list()
gapfactorchange["time"] <- 40 # This is after 30 years of transmission
gapfactorchange["formation.hazard.agegap.gap_factor_man"] <- -1 #-0.8
gapfactorchange["formation.hazard.agegap.gap_factor_woman"] <- -1 #-0.8

interventions <- list(gapfactorchange)


testrunwithintervention <- simpact.run(cfg,
                                       DestDir,
                                       intervention = interventions,
                                       identifierFormat = identifier)

# Reading in the data
outputID <- testrunwithintervention["id"]
datalist <- readthedata(DestDir, outputID) # a list that contains ptable, rtable and etable

# Viewing the simulation output
allPrevalencetimestepsdata <- prevalenceheatmapdata(datalist$ptable, cfg)
allIncidencetimestepsdata <- incidenceheatmapdata(datalist$ptable, datalist$etable, cfg)
#Rtimestepsdata <- allRtimestepsdata(datalist$ptable, datalist$etable, cfg)
ageandtimeattransmissiondata <- transmissionmap(datalist$ptable, datalist$etable, cfg)
agemixingdata <- agemixing(datalist$ptable, datalist$rtable)

# age mixing pattern for period 10-40 years
# of couples at the time of relationship formation

agescatter1 <- agemixingdata$agescatterdata[agemixingdata$agescatterdata$FormTime >= 10 & agemixingdata$agescatterdata$FormTime < 40, ]
matplot(agescatter1$AgeMaleatForm,
        agescatter1$AgeFemaleatForm,
        type = "p", pch=16, xlab = "Male age", ylab = "Female age",
        cex.lab = 1, cex.axis= 1,
        xlim=c(10,100), ylim=c(10,100))
lines(c(15,100), c(15,100), lty="dashed", col="orange2",lwd=8)

# age mixing pattern for period 40-70 years
agescatter2 <- agemixingdata$agescatterdata[agemixingdata$agescatterdata$FormTime >= 25, ]
matplot(agescatter2$AgeMaleatForm,
        agescatter2$AgeFemaleatForm,
        type = "p", pch=16, xlab = "Male age", ylab = "Female age",
        cex.lab = 1, cex.axis= 1,
        xlim=c(10,100), ylim=c(10,100))
lines(c(15,100), c(15,100), lty="dashed", col="orange2",lwd=8)


# Creating a variable that indicates if the relationship was formed before or after the intervention
agemixingdata$agescatterdata$intervention <- as.factor(agemixingdata$agescatterdata$FormTime >= 40)
levels(agemixingdata$agescatterdata$intervention) <- c("Before", "After")
# Or a 2d density plot, which can handle the many superimposed points better.

ggplot(data = agemixingdata$agescatterdata, aes(x = AgeMaleatForm,
                                                y = AgeFemaleatForm)) +
  geom_bin2d(binwidth = c(1, 1)) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed",
              size = 2,
              colour = "orange2") +
  facet_wrap(~ intervention, nrow = 1) +
  xlab("Male age") +
  ylab("Female age")



# Now let's visualise the HIV prevalence model output
my.pal <- c("#FFFFFF", brewer.pal(9,"YlOrRd"))


# allPrevalencetimestepsdataM <- allPrevalencetimestepsdata[Gender==0, ]
# p <- ggplot(allPrevalencetimestepsdataM, aes(time_i, Age)) +
#   geom_tile(aes(fill = Prevalence), colour = "white") 
# p + scale_fill_gradient(low = "white", high = "steelblue")
# 
# # And now we look at HIV prevalence over time and by age, in women:
# allPrevalencetimestepsdataW <- allPrevalencetimestepsdata[Gender==1, ]
# p <- ggplot(allPrevalencetimestepsdataW, aes(time_i, Age)) +
#   geom_tile(aes(fill = Prevalence), colour = "white") 
# p + scale_fill_gradient(low = "white", high = "darkred")

# We can create a 2-panel plot to look at HIV prevalence in men and women
allPrevalencetimestepsdata$Gender <- as.factor(allPrevalencetimestepsdata$Gender)
levels(allPrevalencetimestepsdata$Gender) <- c("Men", "Women")

p <- ggplot(allPrevalencetimestepsdata, aes(time_i, Age, Prevalence)) +
  geom_tile(aes(fill = Prevalence))
p + scale_fill_gradientn(colours = my.pal) +
  facet_wrap(~Gender, nrow=1) +
  theme(panel.grid = element_blank()) +
  xlab("Simulation time")

# Now let's visualise the HIV incidence model output
allIncidencetimestepsdata$Gender <- as.factor(allIncidencetimestepsdata$Gender)
levels(allIncidencetimestepsdata$Gender) <- c("Men", "Women")
# We don't want weirdly high incidence values for strata where the number of PY
# exposure time was very small

p <- ggplot(allIncidencetimestepsdata[allIncidencetimestepsdata$PY_G_A>=5, ], aes(time_i, Age, Incidence_G_A)) +
  geom_tile(aes(fill = Incidence_G_A))
p + scale_fill_gradientn(colours = my.pal) +
  facet_wrap(~Gender, nrow=1) +
  theme(panel.grid = element_blank()) +
  xlab("Simulation time")


# But we are also interested in the population-average HIV incidence over time
names(Incidencetimestepdata)
p <- ggplot(allIncidencetimestepsdata, aes(x=time_i, y=Incidence_G, group=1)) +
  geom_line() +
  xlab("Simulation time") +
  ylab(" HIV incidence") +
  facet_wrap(~Gender, nrow=1)
print(p)



# VISUALISING THE HIV TRANSMISSION NETWORK
# Since we keep track of all the events, we can visualise when exactly HIV was
# transmitted and how old the infectors and newly infecteds were at the time of
# transmission.

agetrans <- ageandtimeattransmissiondata$transmissiondata

plot(c(0, 75), c(0, 95), type="n",
     xlab="Simulation time",
     ylab="Ages of transmission pairs")

Arrows(agetrans$timeattransm,
       agetrans$ageoftransmitter,
       agetrans$timeattransm,
       agetrans$ageofreceptor,
       arr.length = 0.2, code = 2,
       arr.type = "triangle", col = - agetrans$genderoftransmitter + 3)
legend("topleft",
       legend = c("male-to-female", "female-to-male"),
       fill = c(3,2),
       border = c(3,2))









# COMPARING CONTRASTING SCENARIOS

# In this final illustration, we will compare two senarios. In the first
# scenario, HIV is introduced in a population with a very wide age-mixing
# pattern, and the evolving HIV epidemic is simulated for 60 years. The second
# scenario is identical to the first one, for the first 30 years of the
# epidemic. At that point (i.e. 40 years into the simulation because HIV is only
# introduced after 10 years of simulation time), we introduce an intervention
# that narrows down the age-mixing scatter. However, the overall rate of
# relationship formation does not change.

# To ensure that the simulation output for the 2 scenarios is identical for the
# first 30 years, we will set an explicit "seed" value of the random number
# generator. By setting the same seed in both scenarios, we ensure that the
# series of random numbers generated are the same in both scenarios.

# First, we run ABC to find out what the baseline factor needs to change into at
# the time of the intervention.

# Our target statistics are:
# 1. Number of partners per year of 1 in the period 10-40 years into the simulation
# 2. Number of partners per year of 1 in the period 40-70 years into the simulation

# 3. Standard deviation of age gaps of 3 years in the period 10-20 years into the simulation
# 4. Standard deviation of age gaps of 1.5 years in the period 20-30 years into the simulation


ABC_DestDir <- "/tmp/ABC/"

simpact4ABCintervention <- function(inputvector){
  cfg <- list()
  # Let's start with a small population first, to keep computation time low
  cfg["population.numwomen"] <- 200
  cfg["population.nummen"] <- 200
  cfg["population.simtime"] <- 7 # 70
  cfg["person.agegap.man.dist.type"] <- "fixed"
  cfg["person.agegap.woman.dist.type"] <- "fixed"
  cfg["person.agegap.man.dist.fixed.value"] <- 3
  cfg["person.agegap.woman.dist.fixed.value"] <- 3
  cfg["formation.hazard.agegap.baseline"] <- inputvector[1]
  cfg["formation.hazard.agegap.gap_factor_man"] <- inputvector[2]
  cfg["formation.hazard.agegap.gap_factor_woman"] <- inputvector[2]
  cfg["hivseed.time"] <- 1 # 10
  # We want to have precise control over the number (not just the fraction) of HIV
  # seeds being introduced
  cfg["hivseed.type"] <- "amount"
  cfg["hivseed.amount"] <- 4     # 4 out of 400 is 1% of the population
  # HIV is introduced at an age that the individual is likely already sexually
  # active, in a relationship, and still has lots of time to infect others
  cfg["hivseed.age.min"] <- 20
  cfg["hivseed.age.max"] <- 30
  
  # Nobody will start ART during the simulations
  cfg["monitoring.cd4.threshold"] <- 0.001
  
  
  # We are inflating the transmission parameters "for effect"
  cfg["transmission.param.a"] <- -0.3
  cfg["transmission.param.b"] <- -8  
  cfg["transmission.param.c"] <- 0.1649
  
  gapfactorchange <- list()
  gapfactorchange["time"] <- 4 # 40 # This is after 30 years of transmission
  gapfactorchange["formation.hazard.agegap.baseline"] <- inputvector[3]
  gapfactorchange["formation.hazard.agegap.gap_factor_man"] <- inputvector[4]
  gapfactorchange["formation.hazard.agegap.gap_factor_woman"] <- inputvector[4]
  
  interventions <- list(gapfactorchange)
  
  
  
  
  resultsABCintervention <- simpact.run(cfg,
                                        "/tmp/testdirectory/",
                                        seed = 1,
                                        parallel = TRUE,
                                        intervention = interventions)
  
  ABC_rtable <- fread(resultsABCintervention["logrelations"], sep = ",", skip = 0) # fread(relationlogfilename, sep = ",", skip = 0)
  ABC_ptable <- fread(resultsABCintervention["logpersons"], sep = ",", skip = 0) # fread(personlogfilename, sep = ",", skip = 0)
  relsperpersonperyearbefore <- nrow(ABC_rtable[FormTime < 4]) / (nrow(ABC_ptable[TOB < 4])/2) / 4
  agegapsdbefore <- sd(ABC_rtable[FormTime < 4]$AgeGap)
  relsperpersonperyearafter <- nrow(ABC_rtable[FormTime >= 4]) / (nrow(ABC_ptable[TOB < 4])/2) / 4
  agegapsdafter <- sd(ABC_rtable[FormTime >= 4]$AgeGap)  
  outputvector <- c(relsperpersonperyearbefore, agegapsdbefore, relsperpersonperyearafter, agegapsdafter)
  return(outputvector)
}




simpact_prior <- list(c("unif", 3.2, 4.2),
                      c("unif", -0.27, -0.21),
                      c("unif", 3.9, 4.7),
                      c("unif", -0.52, -0.42))


# Lastly, we specify the target summary statistic
sum_stat_obs <- c(1, 3, 1, 1.5)


# Initial number of simulations
n_init <- 50
alpha <- 0.5 # This is the fraction of n_init that must be retained by the end of the procedure
pacc <- 0.2 #0.02 # This is the criterion for how close the output must approximate the target statistics
# pacc must be between 0 and 1. The smaller, the more strict the criterion.

ABC_LenormandResultintervention <- ABC_sequential(method="Lenormand",
                                      model=simpact4ABCintervention,
                                      prior=simpact_prior,
                                      nb_simul=n_init,
                                      summary_stat_target=sum_stat_obs,
                                      alpha=alpha,
                                      p_acc_min=pacc,
                                      verbose=FALSE)

ABC_LenormandResultintervention

hist(ABC_LenormandResultintervention$param[,1])
hist(ABC_LenormandResultintervention$param[,2])
hist(ABC_LenormandResultintervention$param[,3])
hist(ABC_LenormandResultintervention$param[,4])

# Highest univariate frequencies:
# 3.65, -0.235, 4.45, -0.5



###################
# NOW WE ARE READY FOR THE ACTUAL COMPARISON
###################


simID <- "7A"
identifierA <- paste0("%T-%y-%m-%d-", simID)

# Just to be complete, we reconstruct the cfg list from scratch
cfg <- list()
cfg["population.numwomen"] <- 2500
cfg["population.nummen"] <- 2500
cfg["population.eyecap.fraction"] <- 0.5
cfg["person.agegap.man.dist.type"] <- "fixed"
cfg["person.agegap.woman.dist.type"] <- "fixed"
cfg["person.agegap.man.dist.fixed.value"] <- 3
cfg["person.agegap.woman.dist.fixed.value"] <- 3
cfg["formation.hazard.agegap.baseline"] <- 3.65
cfg["formation.hazard.agegap.gap_factor_man"] <- -0.235
cfg["formation.hazard.agegap.gap_factor_woman"] <- -0.235

cfg["population.simtime"] <- 70 # 40 is 10 years of "burn in", followed by 30 years of HIV transmission
cfg["hivseed.time"] <- 10
# We want to have precise control over the number (not just the fraction) of HIV
# seeds being introduced
cfg["hivseed.type"] <- "amount"
cfg["hivseed.amount"] <- 10     # 10 out of 1000 is 1% of the population
# HIV is introduced at an age that the individual is likely already sexually
# active, in a relationship, and still has lots of time to infect others
cfg["hivseed.age.min"] <- 20
cfg["hivseed.age.max"] <- 30

# Nobody will start ART during the simulations
cfg["monitoring.cd4.threshold"] <- 0.001


# We are inflating the transmission parameters "for effect"
cfg["transmission.param.a"] <- -0.35 # -1.3997 # -0.3
cfg["transmission.param.b"] <- -9 # -12.022 # -8  
cfg["transmission.param.c"] <- 0.1649

gapfactorchange <- list()
gapfactorchange["time"] <- 70 # This is after 30 years of transmission
gapfactorchange["formation.hazard.agegap.baseline"] <- 4.45
gapfactorchange["formation.hazard.agegap.gap_factor_man"] <- -0.5 #-0.8
gapfactorchange["formation.hazard.agegap.gap_factor_woman"] <- -0.5 #-0.8

interventionsA <- list(gapfactorchange)


testrunA <- simpact.run(cfg,
                        DestDir,
                        seed = 1,
                        parallel = TRUE,
                        intervention = interventionsA,
                        identifierFormat = identifierA)


# Reading in the data
outputID <- testrunA["id"]
datalist <- readthedata(DestDir, outputID) # a list that contains ptable, rtable and etable

# Viewing the simulation output
allPrevalencetimestepsdata <- prevalenceheatmapdata(datalist$ptable, cfg)
allIncidencetimestepsdata <- incidenceheatmapdata(datalist$ptable, datalist$etable, cfg)
#Rtimestepsdata <- allRtimestepsdata(datalist$ptable, datalist$etable, cfg)
ageandtimeattransmissiondata <- transmissionmap(datalist$ptable, datalist$etable, cfg)
agemixingdata <- agemixing(datalist$ptable, datalist$rtable)

# age mixing pattern for period 10-40 years
# of couples at the time of relationship formation

agescatter1 <- agemixingdata$agescatterdata[agemixingdata$agescatterdata$FormTime >= 10 & agemixingdata$agescatterdata$FormTime < 40, ]
matplot(agescatter1$AgeMaleatForm,
        agescatter1$AgeFemaleatForm,
        type = "p", pch=16, xlab = "Male age", ylab = "Female age",
        cex.lab = 1, cex.axis= 1,
        xlim=c(10,100), ylim=c(10,100))
lines(c(15,100), c(15,100), lty="dashed", col="orange2",lwd=8)

# age mixing pattern for period 40-70 years
agescatter2 <- agemixingdata$agescatterdata[agemixingdata$agescatterdata$FormTime >= 25, ]
matplot(agescatter2$AgeMaleatForm,
        agescatter2$AgeFemaleatForm,
        type = "p", pch=16, xlab = "Male age", ylab = "Female age",
        cex.lab = 1, cex.axis= 1,
        xlim=c(10,100), ylim=c(10,100))
lines(c(15,100), c(15,100), lty="dashed", col="orange2",lwd=8)


# Creating a variable that indicates if the relationship was formed before or after the intervention
agemixingdata$agescatterdata$intervention <- as.factor(agemixingdata$agescatterdata$FormTime >= 40)
levels(agemixingdata$agescatterdata$intervention) <- c("Before", "After")
# Or a 2d density plot, which can handle the many superimposed points better.

ggplot(data = agemixingdata$agescatterdata, aes(x = AgeMaleatForm,
                                                y = AgeFemaleatForm)) +
  geom_bin2d(binwidth = c(1, 1)) +
  scale_fill_gradientn(limits=c(0,5000), breaks=seq(0, 5000, by=500), colours=rainbow(10)) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed",
              size = 2,
              colour = "orange2") +
  facet_wrap(~ intervention, nrow = 1) +
  xlab("Male age") +
  ylab("Female age")



# Now let's visualise the HIV prevalence model output
my.pal <- c("#FFFFFF", brewer.pal(9,"YlOrRd"))


# allPrevalencetimestepsdataM <- allPrevalencetimestepsdata[Gender==0, ]
# p <- ggplot(allPrevalencetimestepsdataM, aes(time_i, Age)) +
#   geom_tile(aes(fill = Prevalence), colour = "white") 
# p + scale_fill_gradient(low = "white", high = "steelblue")
# 
# # And now we look at HIV prevalence over time and by age, in women:
# allPrevalencetimestepsdataW <- allPrevalencetimestepsdata[Gender==1, ]
# p <- ggplot(allPrevalencetimestepsdataW, aes(time_i, Age)) +
#   geom_tile(aes(fill = Prevalence), colour = "white") 
# p + scale_fill_gradient(low = "white", high = "darkred")

# We can create a 2-panel plot to look at HIV prevalence in men and women
allPrevalencetimestepsdata$Gender <- as.factor(allPrevalencetimestepsdata$Gender)
levels(allPrevalencetimestepsdata$Gender) <- c("Men", "Women")

p <- ggplot(allPrevalencetimestepsdata, aes(time_i, Age, Prevalence)) +
  geom_tile(aes(fill = Prevalence))
p + scale_fill_gradientn(colours = my.pal) +
  facet_wrap(~Gender, nrow=1) +
  theme(panel.grid = element_blank()) +
  xlab("Simulation time")

# Now let's visualise the HIV incidence model output
allIncidencetimestepsdata$Gender <- as.factor(allIncidencetimestepsdata$Gender)
levels(allIncidencetimestepsdata$Gender) <- c("Men", "Women")
# We don't want weirdly high incidence values for strata where the number of PY
# exposure time was very small

p <- ggplot(allIncidencetimestepsdata[allIncidencetimestepsdata$PY_G_A>=5, ], aes(time_i, Age, Incidence_G_A)) +
  geom_tile(aes(fill = Incidence_G_A))
p + scale_fill_gradientn(colours = my.pal) +
  facet_wrap(~Gender, nrow=1) +
  theme(panel.grid = element_blank()) +
  xlab("Simulation time")


# But we are also interested in the population-average HIV incidence over time
names(Incidencetimestepdata)
p <- ggplot(allIncidencetimestepsdata, aes(x=time_i, y=Incidence_G, group=1)) +
  geom_line() +
  xlab("Simulation time") +
  ylab(" HIV incidence") +
  facet_wrap(~Gender, nrow=1)
print(p)



# VISUALISING THE HIV TRANSMISSION NETWORK
# Since we keep track of all the events, we can visualise when exactly HIV was
# transmitted and how old the infectors and newly infecteds were at the time of
# transmission.

agetrans <- ageandtimeattransmissiondata$transmissiondata

plot(c(0, 75), c(0, 95), type="n",
     xlab="Simulation time",
     ylab="Ages of transmission pairs")

Arrows(agetrans$timeattransm,
       agetrans$ageoftransmitter,
       agetrans$timeattransm,
       agetrans$ageofreceptor,
       arr.length = 0.2, code = 2,
       arr.type = "triangle", col = - agetrans$genderoftransmitter + 3)
legend("topleft",
       legend = c("male-to-female", "female-to-male"),
       fill = c(3,2),
       border = c(3,2))







simID <- "7B"
identifierB <- paste0("%T-%y-%m-%d-", simID)

# Just to be complete, we reconstruct the cfg list from scratch
cfgB <- list()
cfgB["population.numwomen"] <- 1000
cfgB["population.nummen"] <- 1000
cfgB["person.agegap.man.dist.type"] <- "fixed"
cfgB["person.agegap.woman.dist.type"] <- "fixed"
cfgB["person.agegap.man.dist.fixed.value"] <- 3
cfgB["person.agegap.woman.dist.fixed.value"] <- 3
cfgB["formation.hazard.agegap.baseline"] <- 3.65
cfgB["formation.hazard.agegap.gap_factor_man"] <- -0.235
cfgB["formation.hazard.agegap.gap_factor_woman"] <- -0.235

cfgB["population.simtime"] <- 70 # 40 is 10 years of "burn in", followed by 30 years of HIV transmission
cfgB["hivseed.time"] <- 10
# We want to have precise control over the number (not just the fraction) of HIV
# seeds being introduced
cfgB["hivseed.type"] <- "amount"
cfgB["hivseed.amount"] <- 10     # 10 out of 1000 is 1% of the population
# HIV is introduced at an age that the individual is likely already sexually
# active, in a relationship, and still has lots of time to infect others
cfgB["hivseed.age.min"] <- 20
cfgB["hivseed.age.max"] <- 30

# Nobody will start ART during the simulations
cfgB["monitoring.cd4.threshold"] <- 0.001


# We are inflating the transmission parameters "for effect"
cfgB["transmission.param.a"] <- -1.3997
cfgB["transmission.param.b"] <- -12.022  
cfgB["transmission.param.c"] <- 0.1649

# The B arm of the simulation does not have an intervention

ptm <- proc.time()
testrunB <- simpact.run(cfgB,
                        DestDir,
                        seed = 1,
                        parallel = TRUE,
                        identifierFormat = identifierB)

proc.time() - ptm


# Reading in the data
outputID <- testrunB["id"]
datalist <- readthedata(DestDir, outputID) # a list that contains ptable, rtable and etable

# Viewing the simulation output
allPrevalencetimestepsdata <- prevalenceheatmapdata(datalist$ptable, cfgB)
allIncidencetimestepsdata <- incidenceheatmapdata(datalist$ptable, datalist$etable, cfgB)
#Rtimestepsdata <- allRtimestepsdata(datalist$ptable, datalist$etable, cfgB)
ageandtimeattransmissiondata <- transmissionmap(datalist$ptable, datalist$etable, cfgB)
agemixingdata <- agemixing(datalist$ptable, datalist$rtable)

# age mixing pattern for period 10-40 years
# of couples at the time of relationship formation

agescatter1 <- agemixingdata$agescatterdata[agemixingdata$agescatterdata$FormTime >= 10 & agemixingdata$agescatterdata$FormTime < 40, ]
matplot(agescatter1$AgeMaleatForm,
        agescatter1$AgeFemaleatForm,
        type = "p", pch=16, xlab = "Male age", ylab = "Female age",
        cex.lab = 1, cex.axis= 1,
        xlim=c(10,100), ylim=c(10,100))
lines(c(15,100), c(15,100), lty="dashed", col="orange2",lwd=8)

# age mixing pattern for period 40-70 years
agescatter2 <- agemixingdata$agescatterdata[agemixingdata$agescatterdata$FormTime >= 25, ]
matplot(agescatter2$AgeMaleatForm,
        agescatter2$AgeFemaleatForm,
        type = "p", pch=16, xlab = "Male age", ylab = "Female age",
        cex.lab = 1, cex.axis= 1,
        xlim=c(10,100), ylim=c(10,100))
lines(c(15,100), c(15,100), lty="dashed", col="orange2",lwd=8)


# Creating a variable that indicates if the relationship was formed before or after the intervention
agemixingdata$agescatterdata$intervention <- as.factor(agemixingdata$agescatterdata$FormTime >= 40)
levels(agemixingdata$agescatterdata$intervention) <- c("Before", "After")
# Or a 2d density plot, which can handle the many superimposed points better.

ggplot(data = agemixingdata$agescatterdata, aes(x = AgeMaleatForm,
                                                y = AgeFemaleatForm)) +
  geom_bin2d(binwidth = c(1, 1)) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed",
              size = 2,
              colour = "orange2") +
  facet_wrap(~ intervention, nrow = 1) +
  xlab("Male age") +
  ylab("Female age")



# Now let's visualise the HIV prevalence model output
my.pal <- c("#FFFFFF", brewer.pal(9,"YlOrRd"))


# allPrevalencetimestepsdataM <- allPrevalencetimestepsdata[Gender==0, ]
# p <- ggplot(allPrevalencetimestepsdataM, aes(time_i, Age)) +
#   geom_tile(aes(fill = Prevalence), colour = "white") 
# p + scale_fill_gradient(low = "white", high = "steelblue")
# 
# # And now we look at HIV prevalence over time and by age, in women:
# allPrevalencetimestepsdataW <- allPrevalencetimestepsdata[Gender==1, ]
# p <- ggplot(allPrevalencetimestepsdataW, aes(time_i, Age)) +
#   geom_tile(aes(fill = Prevalence), colour = "white") 
# p + scale_fill_gradient(low = "white", high = "darkred")

# We can create a 2-panel plot to look at HIV prevalence in men and women
allPrevalencetimestepsdata$Gender <- as.factor(allPrevalencetimestepsdata$Gender)
levels(allPrevalencetimestepsdata$Gender) <- c("Men", "Women")

p <- ggplot(allPrevalencetimestepsdata, aes(time_i, Age, Prevalence)) +
  geom_tile(aes(fill = Prevalence))
p + scale_fill_gradientn(colours = my.pal) +
  facet_wrap(~Gender, nrow=1) +
  theme(panel.grid = element_blank()) +
  xlab("Simulation time")

# Now let's visualise the HIV incidence model output
allIncidencetimestepsdata$Gender <- as.factor(allIncidencetimestepsdata$Gender)
levels(allIncidencetimestepsdata$Gender) <- c("Men", "Women")
# We don't want weirdly high incidence values for strata where the number of PY
# exposure time was very small

p <- ggplot(allIncidencetimestepsdata[allIncidencetimestepsdata$PY_G_A>=5, ], aes(time_i, Age, Incidence_G_A)) +
  geom_tile(aes(fill = Incidence_G_A))
p + scale_fill_gradientn(colours = my.pal) +
  facet_wrap(~Gender, nrow=1) +
  theme(panel.grid = element_blank()) +
  xlab("Simulation time")


# But we are also interested in the population-average HIV incidence over time
names(Incidencetimestepdata)
p <- ggplot(allIncidencetimestepsdata, aes(x=time_i, y=Incidence_G, group=1)) +
  geom_line() +
  xlab("Simulation time") +
  ylab(" HIV incidence") +
  facet_wrap(~Gender, nrow=1)
print(p)



# VISUALISING THE HIV TRANSMISSION NETWORK
# Since we keep track of all the events, we can visualise when exactly HIV was
# transmitted and how old the infectors and newly infecteds were at the time of
# transmission.

agetrans <- ageandtimeattransmissiondata$transmissiondata

plot(c(0, 75), c(0, 95), type="n",
     xlab="Simulation time",
     ylab="Ages of transmission pairs")

Arrows(agetrans$timeattransm,
       agetrans$ageoftransmitter,
       agetrans$timeattransm,
       agetrans$ageofreceptor,
       arr.length = 0.2, code = 2,
       arr.type = "triangle", col = - agetrans$genderoftransmitter + 3)
legend("topleft",
       legend = c("male-to-female", "female-to-male"),
       fill = c(3,2),
       border = c(3,2))








# From here on is old code. Please ignore







# Now making the 5 age mixing plots

agemixing <- function(DT = datalist$ptable, DTR = datalist$rtable){
  DTRandDT_a <- merge(data.frame(DTR), data.frame(DT), by.x = "IDm", by.y = "ID")
  DTRandDT_b <- merge(DTRandDT_a, DT, by.x = "IDw", by.y = "ID", suffixes = c(".m", ".w"))
  DTRandDT_b$relID <- factor(1:nrow(DTRandDT_b))
  setnames(DTRandDT_b, c("IDw", "IDm"), c("ID.w", "ID.m"))
  DTRlong <- reshape(data = DTRandDT_b,
                     idvar = "relID",
                     varying = names(DTRandDT_b)[c(1:2, 6:29)],
                     timevar = "g",
                     direction = "long")
  DTRlong$g.ID <- as.character(paste(DTRlong$g, DTRlong$ID, sep="."))
  
  # Analysis of age mixing pattern, using the DTRlong dataset
  # Mean age gap, between-, and within-subject variance of age gaps
  meanagegap <- mean(DTRlong$AgeGap)
  varagegap <- var(DTRlong$AgeGap[DTRlong$g == "w"]) # We only want to count each relationship once
  MixingModel <- lme(AgeGap~1, random=~1|g.ID, data = DTRlong)
  VarianceEstimates <- VarCorr(MixingModel)
  BetweenVar <- as.numeric(VarianceEstimates[1])
  WithinVar <- as.numeric(VarianceEstimates[2])
  agemixingdata = data.frame(meanagegap = meanagegap,
                             varagegap = varagegap,
                             BetweenVar = BetweenVar,
                             WithinVar = WithinVar)
  return(agemixingdata)
}







relationlogfilename <- paste0("simpact-cyan-2015-02-05-", simID, "relationlog.csv")
eventlogfilename <- paste0("simpact-cyan-2015-02-05-", simID, "eventlog.csv")

ptable <- fread(personlogfilename, sep = ",", skip = 0)
rtable <- fread(relationlogfilename, sep = ",", skip = 0)



# 1. We want to replace the default age distribution with one that is consistent with the non-HIV-related mortality hazard.
# mortality.normal.weibull.shape                               = 4
# mortality.normal.weibull.scale                               = 70
# mortality.normal.weibull.genderdiff                          = 0 (default is 5 and will need to be changed to 0 further down when the cfg object is created)

agedist.create <- function(shape = 4, scale = 70){
  agebins <- seq(0.5, 99.5, 1)
  probofstillalive <- 1 - pweibull(agebins, shape = shape, scale = scale)
  #plot(agebins,probofstillalive, type="l")
  #meanweib <- 70 * gamma(1+ (1/4))
  fractionsinagebins <- 100 * probofstillalive/sum(probofstillalive)
  simple.age.data.frame <- data.frame(Age = c(agebins, 100.5), Percent.Male = c(fractionsinagebins, 0), Percent.Female = c(fractionsinagebins, 0))
  
  return(simple.age.data.frame)
}

# 2. We write a function to create the configuration files that contain the specific parameter values for the 5 age-mixing patterns
cfg.create <- function(scen = 1:5){
  
  cfg <- list()
  #simpact.getconfig(cfg)
  cfg["population.numwomen"] <- 500
  cfg["population.nummen"] <- 500
  
  cfg["population.simtime"] <- 120#70#10#45   120 is 20 years of burn in, followed by 100 years of HIV transmission
  cfg["hivseed.time"] <- 20
  cfg["hivseed.type"] <- "amount" # We want to have precise control over the number (not just the fraction) of HIV seeds being introduced
  cfg["hivseed.amount"] <- 10     # 10 out of 10'000 is 0.1% of the population
  cfg["hivseed.age.min"] <- 20    # This should be at a time that the individual is likely already sexually active, in a relationship, and still has lots of time to infect others
  cfg["hivseed.age.max"] <- 25
  
  # Nobody will start ART during the simulations
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
  return(cfg)
}




# We first need to find the parameter value for the formation.hazard.agegap.baseline
# that gives us the desired number of relationships, for each of the age-mixing patterns

n <- 10 # number of repetitions per scenario
nscen <- 5
#scenario <- 1:nscen # number of scenarios that are being compared

agegap.baseline <- rep(5, 5) # This is a naive starting value, which will be overwritten # c(6.5, 7.25, 6.55, 8.55, 4.0) #5.4)
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


simpact4ABC <- function(inputvector){
  cfg <- cfg.create(scen = 1:nscen)
  # Now we overwrite the naive estimate of the baseline parameter with a guess from the prior distribution
  cfg["formation.hazard.agegap.baseline"] <- inputvector[1] # This is the parameter we estimate with EasyABC
  results <- simpact.run(cfg, "/tmp/testdirectory")
  datalist <- readthedata(simID) # a list that contains ptable, rtable and etable
  
  
  rABC <- read.csv(results["logrelations"])
  outputvector <- nrow(rABC)
}


# We specify the prior distributions for the input parameter
simpact_prior <- list(c("unif", 8.4, 8.75)) #6.5 is good for scen1, 7.25 for scen2, 6.55 is good for scen3, 7.55 is good for scen4

# Lastly, we specify the target summary statistic
# Target is to have 10000 relationships at the end of the simulation
sum_stat_obs <- 5000

# Now we try to run a sequential ABC scheme, according to the method proposed by Lenormand et al. 2012
# (Adaptive ABC)

# Initial number of simulations
n_init <- 30
alpha <- 0.5
pacc <- 0.2 #0.02

ABC_LenormandResult <- ABC_sequential(method="Lenormand",
                                      model=simpact4ABC,
                                      prior=simpact_prior,
                                      nb_simul=n_init,
                                      summary_stat_target=sum_stat_obs,
                                      alpha=alpha,
                                      p_acc_min=pacc,
                                      verbose=FALSE)

ABC_Lenormand


