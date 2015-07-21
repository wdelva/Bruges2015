simpact4ABC <- function(inputvector){
  
  
  cfg <- list()
  #simpact.getconfig(cfg)
  cfg["population.numwomen"] <- 500
  cfg["population.nummen"] <- 500
  
  cfg["population.simtime"] <- 10
  cfg["hivseed.time"] <- 5
  cfg["hivseed.fraction"] <- 0.05
  cfg["hivtest.cd4.threshold"] <- 0
  cfg["hivtest.interval.dist.uniform.min"] <- 200  # So that nobody actually gets tested and hence nobody starts ART
  cfg["hivtest.interval.dist.uniform.max"] <- 201 
  
    
  cfg["transmission.param.a"] <- -0.3
  cfg["transmission.param.b"] <- -8  
  cfg["transmission.param.c"] <- 0.1649
  
  
  cfg["formation.hazard.agegap.numrel_man"] <- -1
  cfg["formation.hazard.agegap.numrel_woman"] <- -1
  
  
  cfg["person.agegap.man.dist.type"] <- person.agegap.man.dist.type[scen]
  cfg["person.agegap.woman.dist.type"] <- person.agegap.woman.dist.type[scen]
  
  if (person.agegap.man.dist.type[scen] == "normal"){
    cfg["person.agegap.man.dist.normal.mu"] <- person.agegap.man.dist.normal.mu[scen]
    cfg["person.agegap.man.dist.normal.sigma"] <- person.agegap.man.dist.normal.sigma[scen]    
    cfg["person.agegap.woman.dist.normal.mu"] <- person.agegap.woman.dist.normal.mu[scen]
    cfg["person.agegap.woman.dist.normal.sigma"] <- person.agegap.woman.dist.normal.sigma[scen]
  } else {
    cfg["person.agegap.man.dist.fixed.value"] <- person.agegap.man.dist.fixed.value[scen]
    cfg["person.agegap.woman.dist.fixed.value"] <- person.agegap.woman.dist.fixed.value[scen]    
  }
  
  cfg["formation.hazard.agegap.gap_factor_man"] <- agegap.gap_factor_man[scen]
  cfg["formation.hazard.agegap.gap_factor_woman"] <- agegap.gap_factor_woman[scen]
  cfg["formation.hazard.agegap.gap_agescale_man"] <- agegap.gap_agescale_man[scen]
  cfg["formation.hazard.agegap.gap_agescale_woman"] <- agegap.gap_agescale_woman[scen]
  
  cfg["formation.hazard.agegap.baseline"] <- inputvector[1] # This is the parameter we estimate with EasyABC
  results <- simpact.run(cfg, "/tmp/testdirectory")
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

ABC_Lenormand <- ABC_sequential(method="Lenormand",
                                model=simpact4ABC,
                                prior=simpact_prior,
                                nb_simul=n_init,
                                summary_stat_target=sum_stat_obs,
                                alpha=alpha,
                                p_acc_min=pacc,
                                verbose=FALSE)

ABC_Lenormand
hist(ABC_Lenormand$param[,1])
##hist(ABC_Lenormand$param[,2])

