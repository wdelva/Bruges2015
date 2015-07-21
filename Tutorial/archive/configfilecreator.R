library(RSimpactCyan)
cfg <- list()
DestDir <- "/Users/wimdelva/Dropbox (Personal)/FWO Age Mixing/Papers/Simplest age mixing simulation/configfiles"

# We want to replace the default age distribution with one that is consistent with the non-HIV-related mortality hazard.
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




agedist.data.frame <- agedist.create(shape = 4, scale = 70)

# Designing the experiment
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
nscen <- 5
#scenario <- 1:nscen # number of scenarios that are being compared

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




for (sim in seq(100 * 1:nscen)){
  cfg <- cfg.create(scen = sim/100)
  
  for (rep in 1:n){
    simID <- sim + rep - 1
    identifier <- paste0("%T-%y-%m-%d-", simID)
    testruns <- simpact.run(cfg,
                            DestDir,
                            dryrun = TRUE,
                            identifierFormat = identifier,
                            agedist = agedist.data.frame)
  }
}



str(results0)


help(package="RSimpactCyan")
