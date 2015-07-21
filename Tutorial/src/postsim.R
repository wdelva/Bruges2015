library(data.table)
library(shape)
library(reshape2)
library(nlme)
library(phylobase)
library(sna)


#library(sqldf)

readthedata <- function(modeloutput){
  path <- as.character(modeloutput["outputfile"])
  outputID <- as.character(modeloutput["id"])
  DestDir <- sub(pattern = paste0(outputID, "output.txt"), replacement = "", x = path, fixed=T)
  personlogfilename <- paste0(DestDir, outputID, "personlog.csv")
  relationlogfilename <- paste0(DestDir, outputID, "relationlog.csv")
  eventlogfilename <- paste0(DestDir, outputID, "eventlog.csv")
  treatmentlogfilename <- paste0(DestDir, outputID, "treatmentlog.csv")
    
  ptable <- fread(personlogfilename, sep = ",", skip = 0)
  rtable <- fread(relationlogfilename, sep = ",", skip = 0)
  etable <- fread(eventlogfilename, sep = ",", skip = 0)
  ttable <- fread(treatmentlogfilename, sep  = ",", skip = 0)
  setnames(etable, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"),
           c("eventtime", "eventname", "p1name", "p1ID", "p1gender", "p1age", "p2name", "p2ID", "p2gender", "p2age"))
  outputtables <- return(list(ptable = ptable, rtable = rtable, etable = etable, ttable = ttable))
}

####
# To create heatplot, we need to loop through time steps, and calculate age- and gender-specific prevalence and incidence rates.
prevalenceheatmapdata <- function(DT = datalist$ptable, cfg = configfile){
  output <- data.table()
  timestep.width <- 1
  decplaces <- decimalplaces(timestep.width)
  timesteps <- seq(as.numeric(cfg["hivseed.time"]), as.numeric(cfg["population.simtime"]), by = timestep.width)
  for (time_i in timesteps){
    DTalive.infected <- alive.infected(DT = datalist$ptable, time = time_i, dec.places = decplaces) # First we only take the data of people who were alive at time_i
    Prevalencetimestepdata <- DTalive.infected[ ,Prevalence := sum(Infected) / sum(!is.na(Infected)),by = "Gender,Age"] # And we calculate HIV prevalence, by gender and age group
    #setnames(Prevalencetimestepdata, "V1", "Prevalence")
    Prevalencetimestepdata <- cbind(time_i, Prevalencetimestepdata)
    output <- rbind(output, Prevalencetimestepdata)
  }
  return(output)
}

incidenceheatmapdata <- function(DTP = datalist$ptable, DTE = datalist$etable, cfg = configfile){
  output <- data.table()
  timestep.width <- 1
  decplaces <- decimalplaces(timestep.width)
  timesteps <- seq(as.numeric(cfg["hivseed.time"]), as.numeric(cfg["population.simtime"]), by = timestep.width)
  for (time_i in head(timesteps, -1)){
    DTcases <- incidentcase(DT = DTE, interval = c(time_i, (time_i + timestep.width)))  # New infections that happened at time x with: time_i < x <= time_(i+1)
    DTnonHIVdeaths <- nonHIVdeath(DT = DTE, interval = c(time_i, (time_i + timestep.width))) # New non-HIV deaths that happened at time x with: time_i < x <= time_(i+1)
    DTcases$ID <- DTcases$p2ID
    DTcases$Gender <- DTcases$p2gender
    DTnonHIVdeaths$ID <- DTnonHIVdeaths$p1ID
    DTnonHIVdeaths$Gender <- DTnonHIVdeaths$p1gender
    
    # Personyears at risk during the interval
    DTalive.infected <- alive.infected(DT = DTP, time = time_i, dec.places = decplaces) # First we only take the data of people who were alive at time_i
    DTalive.uninfected <- DTalive.infected[Infected == "FALSE"] # We only keep those that were not infected at the start of the interval
    # A merge with etable will find those that died and/or got infected during the interval
    # To do this, we need to add a key to DTalive.uninfected
    DTPkeycols <- c("ID", "Gender")
    DTcaseskeycols <- c("ID", "Gender") # The ID and gender of the newly infected person  
    DTnonHIVdeathskeycols <- c("ID", "Gender") # The ID and gender of the newly died person (non-HIV death)
    setkeyv(DTalive.uninfected, DTPkeycols)
    setkeyv(DTcases, DTcaseskeycols)
    setkeyv(DTnonHIVdeaths, DTnonHIVdeathskeycols)
    DTforincidencecalc_a <- merge(DTalive.uninfected, DTcases, all.x=TRUE)
    DTforincidencecalc <- merge(DTforincidencecalc_a, DTnonHIVdeaths, all.x=TRUE, suffixes = c(".incident", ".dead"))
    DTforincidencecalc$eventtime.incident[is.na(DTforincidencecalc$eventtime.incident)] <- time_i + timestep.width
    DTforincidencecalc$eventtime.dead[is.na(DTforincidencecalc$eventtime.dead)] <- time_i + timestep.width
    
    DTforincidencecalc$PYend <- pmin(DTforincidencecalc$eventtime.incident, DTforincidencecalc$eventtime.dead)
    DTforincidencecalc$PY <- DTforincidencecalc$PYend - time_i
    
    Incidencetimestepdata <- DTforincidencecalc[ ,Incidence_G_A := sum(eventname.incident=="transmission", na.rm = TRUE) / sum(PY),by = "Gender,Age"] # And we calculate HIV Incidence, by gender and age group
    Incidencetimestepdata <- DTforincidencecalc[ ,PY_G_A := sum(PY),by = "Gender,Age"] # And we calculate PY, by gender and age group
    Incidencetimestepdata <- DTforincidencecalc[ ,Incidence_G := sum(eventname.incident=="transmission", na.rm = TRUE) / sum(PY),by = "Gender"] # And we calculate HIV Incidence, by gender and age group
    Incidencetimestepdata <- DTforincidencecalc[ ,PY_G := sum(PY),by = "Gender"] # And we calculate PY, by gender and age group
    Incidencetimestepdata <- DTforincidencecalc[ ,Incidence_All := sum(eventname.incident=="transmission", na.rm = TRUE) / sum(PY)] # And we calculate HIV Incidence, by gender and age group
    Incidencetimestepdata <- DTforincidencecalc[ ,PY_All := sum(PY)] # And we calculate PY, by gender and age group
    # (so we can be selective when plotting, and not plot if PY per age-gender group is too low)
    # setnames(Incidencetimestepdata, "V1", "Incidence") # Cases per personyear at risk
    Incidencetimestepdata <- cbind(time_i, Incidencetimestepdata)
    output <- rbind(output, Incidencetimestepdata)
  }
  return(output)
}

allRtimestepsdata <- function(DT = datalist$ptable, DTE = datalist$etable, cfg = configfile){
  output <- data.table()
  timestep.width <- 1
  decplaces <- decimalplaces(timestep.width)
  timesteps <- seq(as.numeric(cfg["hivseed.time"]), as.numeric(cfg["population.simtime"]), by = timestep.width)
  for (time_i in head(timesteps, -1)){
    DTnewlyinfected <- DT[InfectTime > time_i & InfectTime <= (time_i + timestep.width) & TOD <= as.numeric(cfg["population.simtime"])] # People who got infected during the time slice and had a completed infectious period 
    # small hack: DT[ID == 151]$TOD <- as.numeric(cfg["population.simtime"]) - 1
    transmissions <- DTE[eventname == "transmission" & eventtime > time_i] # transmissions that could possibly be caused by the pople in DTnewlyinfected
    transmissions$ID <- transmissions$p1ID
    transmissions$Gender <- transmissions$p1gender

    DTnewlyinfectedkeycols <- c("ID", "Gender")
    transmissionskeycols <- c("ID", "Gender") # The ID and gender of the receptors, acquiring infection from people in DTnewlyinfected  
    setkeyv(DTnewlyinfected, DTnewlyinfectedkeycols)
    setkeyv(transmissions, transmissionskeycols)
    DTforRcalc <- merge(DTnewlyinfected, transmissions, all.x=TRUE)
    DTforRcalc$ID.Gender <- paste(DTforRcalc$ID, DTforRcalc$Gender, sep=".")
    DTwithR <- DTforRcalc[ ,infcaused :=  sum(!is.na(eventtime)), by = "ID.Gender"] # number of infections caused for each person in DTnewlyinfected 
    DTwithR$R <- sum(DTwithR$infcaused) / nrow(DTwithR)
    Rtimestepdata <- data.frame(time_i = time_i, suminfcaused = sum(DTwithR$infcaused), R = sum(DTwithR$infcaused) / nrow(DTwithR))
    output <- rbind(output, Rtimestepdata)
  }
  return(output)
}

transmissionmap <- function(DT = datalist$ptable, DTE = datalist$etable, cfg = configfile){
  # requires the shape package
  # xlim <- c(0, 15)
  # ylim <- c(0, 108)
  # plot(0, type = "n", xlim = xlim, ylim = ylim)
  # First we insert diagonal "arrows" (age and time living with HIV)
  Infecteds <- DT[InfectTime < Inf]
  timeatinf <- Infecteds$InfectTime
  ageatinf <- - Infecteds$TOB + Infecteds$InfectTime
  timeatdeath <- pmin(as.numeric(cfg["population.simtime"]), Infecteds$TOD)
  ageatdeath <- ageatinf + (timeatdeath - timeatinf)
  # Arrows(timeatinf, ageatinf, timeatdeath, ageatdeath, arr.length = 0.02, code = 2,
  #        arr.type = "T", arr.col = Infecteds$Gender + 2, col = Infecteds$Gender + 2)
  # Next we insert vertical dashed arrows (age and time at transmission)
  DTE <- datalist$etable
  DTincident <- DTE[eventname == "transmission"]
  timeattransm <- DTincident$eventtime
  ageoftransmitter <- DTincident$p1age
  genderoftransmitter <- DTincident$p1gender
  ageofreceptor <- DTincident$p2age + ((DTincident$p1age - DTincident$p2age) > 0) - ((DTincident$p1age - DTincident$p2age) < 0)
  # Arrows(timeattransm, ageoftransmitter, timeattransm, ageofreceptor, arr.length = 0.1, code = 2,
  #        arr.type = "triangle", col = DTincident$p1gender + 2)
  result <- list(livingwithHIVdata = data.frame(timeatinf = timeatinf,
                                                ageatinf = ageatinf,
                                                timeatdeath = timeatdeath,
                                                ageatdeath = ageatdeath),
                 transmissiondata = data.frame(timeattransm = timeattransm,
                                               ageoftransmitter = ageoftransmitter,
                                               genderoftransmitter = genderoftransmitter,
                                               ageofreceptor = ageofreceptor))
  return(result)
}


agemixing <- function(DT = datalist$ptable, DTR = datalist$rtable){
  DTRandDT_a <- merge(data.frame(DTR), data.frame(DT), by.x = "IDm", by.y = "ID")
  DTRandDT_b <- merge(DTRandDT_a, data.frame(DT), by.x = "IDw", by.y = "ID", suffixes = c(".m", ".w"))
  DTRandDT_b$relID <- factor(1:nrow(DTRandDT_b))
  setnames(DTRandDT_b, c("IDw", "IDm"), c("ID.w", "ID.m"))
  DTRlong <- reshape(data = DTRandDT_b,
                     idvar = "relID",
                     varying = names(DTRandDT_b)[c(1:2, 6:(ncol(DTRandDT_b)-1))], #, 6:ncol(DTRandDT_b))],
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
  
  DTRandDT_m <- merge(data.frame(DTR), data.frame(DT[Gender==0]), by.x = "IDm", by.y = "ID")
  DTRandDT_m$AgeMaleatForm <- - DTRandDT_m$TOB + DTRandDT_m$FormTime
  DTRandDT_m$AgeFemaleatForm <- DTRandDT_m$AgeMaleatForm - DTRandDT_m$AgeGap  
  
  result <- list(agemixingdata = data.frame(meanagegap = meanagegap,
                             varagegap = varagegap,
                             BetweenVar = BetweenVar,
                             WithinVar = WithinVar),
                 agescatterdata = data.frame(DTRandDT_m))
  return(result)
}


# Calculate point prevalence of concurrency,
# Lifetime number of sex partners
# Cross-sectional degree distribution

degreecalc <- function(DT = datalist$ptable, DTR = datalist$rtable){
  DT$ID.Gender <- paste(DT$ID, DT$Gender, sep=".")
  setkey(DT, ID.Gender)
  
  men_degree <- data.table(data.frame(table(DTR$IDm)))
  men_degree$ID.Gender <- paste(men_degree$Var1, 0, sep=".")
  setnames(men_degree, c("IDm", "degree", "ID.Gender"))
  setkey(men_degree, ID.Gender)
  
  women_degree <- data.table(data.frame(table(DTR$IDw)))
  women_degree$ID.Gender <- paste(women_degree$Var1, 1, sep=".")
  setnames(women_degree, c("IDw", "degree", "ID.Gender"))
  setkey(women_degree, ID.Gender)
  
  
  DT_degree_MW <- merge(DT, men_degree, all.x=TRUE)
  DT_degree_W <- merge(DT, women_degree, all.x=TRUE)
  DT_degree_MW$degree[DT_degree_MW$Gender==1] <-  DT_degree_W$degree[DT_degree_W$Gender==1]
  DT_degree_MW$degree[is.na(DT_degree_MW$degree)] <- 0

  
  result <- DT_degree_MW
  return(result)
}
  
# allRtimestepsdata <- function(DT = datalist$ptable, DTE = datalist$etable, cfg = configfile){
#   DTE$ID <- DTE$
#   output <- data.table()
#   
#     transmissions <- DTE[eventname == "transmission" & eventtime > time_i] # transmissions that could possibly be caused by the pople in DTnewlyinfected
#     transmissions$ID <- transmissions$p1ID
#     transmissions$Gender <- transmissions$p1gender
#     
#     DTnewlyinfectedkeycols <- c("ID", "Gender")
#     transmissionskeycols <- c("ID", "Gender") # The ID and gender of the receptors, acquiring infection from people in DTnewlyinfected  
#     setkeyv(DTnewlyinfected, DTnewlyinfectedkeycols)
#     setkeyv(transmissions, transmissionskeycols)
#     DTforRcalc <- merge(DTnewlyinfected, transmissions, all.x=TRUE)
#     DTforRcalc$ID.Gender <- paste(DTforRcalc$ID, DTforRcalc$Gender, sep=".")
#     DTwithR <- DTforRcalc[ ,infcaused :=  sum(!is.na(eventtime)), by = "ID.Gender"] # number of infections caused for each person in DTnewlyinfected 
#     DTwithR$R <- sum(DTwithR$infcaused) / nrow(DTwithR)
#     Rtimestepdata <- data.frame(time_i = time_i, suminfcaused = sum(DTwithR$infcaused), R = sum(DTwithR$infcaused) / nrow(DTwithR))
#     output <- rbind(output, Rtimestepdata)
#   }
#   return(output)
# }

################
# Create cumulative network dataset
################
cumulnetwork <- function(DTR = datalist$rtable, cfg = configfile, windowwidth = 20){
  endwindow <- 30 # as.numeric(cfg["hivseed.time"]) + windowwidth
  startwindow <- 0 # endwindow - windowwidth
  Rels <- DTR[FormTime >=startwindow & FormTime < endwindow]
  rdata <- data.frame(
    tails = paste0(Rels$IDm, "m"),
    heads = paste0(Rels$IDw, "w"),
    form = Rels$FormTime,
    diss = Rels$DisTime,
    durat = pmin(endwindow, Rels$DisTime) - Rels$FormTime,
    stringsAsFactors=FALSE
  )
  cn <- network(rdata[ , 1:2], directed = FALSE, matrix.type = "edgelist") # cumulative network
  cnvertexnames <- network.vertex.names(cn)
  maleindices <- grep("m", cnvertexnames)
  Sex <- rep(0, network.size(cn))
  Sex[maleindices] <- 1
  cnseedindex <- grep(SeedID, cnvertexnames)
  Sex[cnseedindex] <- 2
  VertexSize <- rep(1, network.size(cn))
  VertexSize[cnseedindex] <- 2
  
  set.vertex.attribute(cn, attrname="Sex", value=Sex, v=seq_len(network.size(cn)))
  #cnplot<- plot(cn, label = "", vertex.col = Sex)
  # Saving the plot
  png(file="CN.png", height = 1800, width = 1800, res=300)
  cnplot <- plot(cn, label = "", vertex.col = Sex, vertex.cex = VertexSize)
  dev.off()
  
  
}



################
# Create cross-sectional network dataset
################
csnetwork <- function(DTR = datalist$rtable, cfg = configfile, time = as.numeric(cfg["hivseed.time"])){
  Rels <- DTR[FormTime <= time & DisTime > time]
  rdata <- data.frame(
    tails = paste0(Rels$IDm, "m"),
    heads = paste0(Rels$IDw, "w"),
    form = Rels$FormTime,
    diss = Rels$DisTime,
    durat = pmin(endwindow, Rels$DisTime) - Rels$FormTime,
    stringsAsFactors=FALSE
  )
  csn <- network(rdata[ , 1:2], directed = FALSE, matrix.type = "edgelist") # cumulative network
  csnvertexnames <- network.vertex.names(csn)
  maleindices <- grep("m", csnvertexnames)
  SexCSN <- rep(0, network.size(csn))
  SexCSN[maleindices] <- 1
  csnseedindex <- grep(SeedID, csnvertexnames)
  SexCSN[csnseedindex] <- 2
  VertexSize <- rep(1, network.size(csn))
  VertexSize[csnseedindex] <- 2
  set.vertex.attribute(csn, attrname="Sex", value=SexCSN, v=seq_len(network.size(csn)))
  #csnplot <- plot(csn, label = "", vertex.col = SexCSN)
  # Saving the plot
  png(file="CSN.png", height = 1800, width = 1800, res=300)
  plot(csn, label = "", vertex.col = SexCSN, vertex.cex = VertexSize)
  dev.off()
  
}

# Overlaying cn with csn
# matchingindices <- pmatch(csnvertexnames, cnvertexnames)
# csncoord <- cnplot[matchingindices, ]
# png(file="CSN_CS.png", height = 1800, width = 1800, res=300)
# plot(csn, label = "", vertex.col = SexCSN, coord = csncoord, vertex.cex = VertexSize)
# dev.off()


################
# Create potential transmission network (directed network) (=PTN)
################
# The potential transmission routes, starting from the seed infection
ptnetwork <- function(DT = datalist$ptable, DTR = datalist$rtable, cfg = configfile){
  
  DTR[, IDm := paste0(DTR[, IDm], "m")]
  DTR[, IDw := paste0(DTR[, IDw], "w")]
  Seed <- datalist$ptable[InfectTime == as.numeric(cfg["hivseed.time"])]
  gender_suffix <- c("m", "w")
  SeedID <- paste0(Seed$ID, gender_suffix[Seed$Gender + 1])
  
  PTN <- data.table()
  # Now we run through the partners of the seed and their partners and so on, while NEW relationships are being formed
  DTRfr <- data.frame(DTR)
  # Rels of seed
  SeedRels <- DTRfr[, (Seed$Gender + 1)] == SeedID
  AfterSeedTimeRels <- DTR$DisTime > as.numeric(cfg["hivseed.time"])
  RelsList <- DTR[SeedRels & AfterSeedTimeRels, ] # This is the list of ongoing Rels of the seed, after HIV was introduced
  
  iteration <- 0
  while (nrow(RelsList) > 0){
    iteration <- iteration + 1
    RelsListfr <- data.frame(RelsList)
    RelsListCopy <- data.table(RelsListfr)
    # Now we list the partners of the seed
    if (iteration == 1){
      PartnerColumnIndex <- abs(Seed$Gender - 1) + 1
      
      if (Seed$Gender == 1){
        setcolorder(RelsListCopy, c("IDw", "IDm", "FormTime", "DisTime", "AgeGap"))
      }
    } else {
      PartnerColumnIndex <- - PartnerColumnIndex + 3 # 1 becomes 2 and 2 becomes 1
      if (PartnerColumnIndex == 1){
        setcolorder(RelsListCopy, c("IDw", "IDm", "FormTime", "DisTime", "AgeGap"))
      }
    }
    RelsListCopyfr <- data.frame(RelsListCopy)
    if (iteration == 1){
    PTN <- as.matrix(RelsListCopyfr)#[, 1:2])
    } else {
      PTN <- rbind(as.matrix(PTN), as.matrix(RelsListCopyfr))#[, 1:2]))
    }
    
    PartnerIDs <- unique(RelsListfr[, PartnerColumnIndex])
    FirstRelIndex <- match(PartnerIDs, RelsListfr[, PartnerColumnIndex])
    StartPotentialTransm <- pmax(as.numeric(cfg["hivseed.time"]), RelsListfr[FirstRelIndex, "FormTime"])
    
    # Now list the relationships that these partners formed AFTER they formed the relationship with the seed that could have led to transmission
    RelsList <- data.table()
    for (i in (1:length(PartnerIDs))){
      PartnerRels <- DTRfr[, PartnerColumnIndex] == PartnerIDs[i]
      AfterRels <- DTRfr[, "FormTime"] > StartPotentialTransm[i]
      RelsListi <- DTR[PartnerRels & AfterRels, ]
      RelsList <- rbind(RelsList, RelsListi)
    }
  }
  rPTNdata <- data.frame(
    tails = PTN[, 1],
    heads = PTN[, 2],
    form = as.numeric(PTN[, 3]),
    diss = as.numeric(PTN[, 4]),
    durat = pmin(as.numeric(cfg["population.simtime"]), as.numeric(PTN[, 4])) - as.numeric(PTN[, 3]),
    stringsAsFactors=FALSE
  )
  ptn <- network(rPTNdata[ , 1:2], directed = TRUE, matrix.type = "edgelist") # potential tranmission network
  ptnvertexnames <- network.vertex.names(ptn)
  ptnmaleindices <- grep("m", ptnvertexnames)
  ptnseedindex <- grep(SeedID, ptnvertexnames)
  SexPTN <- rep(0, network.size(ptn))
  SexPTN[ptnmaleindices] <- 1
  SexPTN[ptnseedindex] <- 2
  set.vertex.attribute(ptn, attrname="Sex", value=list(SexPTN), v=seq_len(network.size(ptn)))
  VertexSize <- rep(1, network.size(ptn))
  VertexSize[ptnseedindex] <- 2
  # Plotting the network
  ptnplot <- plot(ptn,
                  label = "",
                  vertex.col = SexPTN,
                  vertex.cex = VertexSize,
                  jitter=T)
  # Saving the plot
  png(file="PTN.png", height = 1800, width = 1800, res=300)
  plot(ptn,
       label = "",
       vertex.col = SexPTN,
       vertex.cex = VertexSize,
       jitter=T)

  dev.off()
  
}


# Overlaying cn with ptn
# matchingindicesPTN <- pmatch(ptnvertexnames, cnvertexnames)
# ptncoord <- cnplot[matchingindicesPTN, ]
# png(file="PTN_CN.png", height = 1800, width = 1800, res=300)
# plot(ptn, label = "", vertex.col = SexPTN, coord = ptncoord, vertex.cex = VertexSize)
# dev.off()







################
# Create transmission tree dataset
################
transmissiontree <- function(DTE = datalist$etable, cfg = configfile){
  TE = DTE[eventname=="transmission"] # TE = transmission events
  edata <- data.frame(
    tails = paste0(TE$p1ID, gender_suffix[TE$p1gender + 1]), #TE$p1name,
    heads = paste0(TE$p2ID, gender_suffix[TE$p2gender + 1]),
    time = TE$eventtime,
    tailname = TE$p1name,
    tailgender = TE$p1gender,
    tailage = TE$p1age,
    headname = TE$p2name,
    headgender = TE$p2gender,
    headage = TE$p2age,
    stringsAsFactors=FALSE
  )
  
# Before we turn this edge list into a network object, we need to add the seeding transmission event
  edataplusroot <- rbind(c("root", edata$tails[1], as.numeric(cfg["hivseed.time"]),
                   "root", NA, NA, edata$tailname[1], edata$tailgender[1],
                   edata$tailage[1] - edata$time[1] + as.numeric(cfg["hivseed.time"])),
                 edata)
#   
#   ntips <- length(unique(c(TE$p1name, TE$p2name)))
#   tip.labels <- unique(c(TE$p1name, TE$p2name))
#                      
#   ninternalnodes <- nrow(edataplusroot)

tn <- network(edata[ , 1:2], directed = TRUE, matrix.type = "edgelist")
tnvertexnames <- network.vertex.names(tn)
tnmaleindices <- grep("m", tnvertexnames)
tnseedindex <- grep(SeedID, tnvertexnames) # Temporary solution for SACEMA Research days
SexTN <- rep(0, network.size(tn))
SexTN[tnmaleindices] <- 1
SexTN[tnseedindex] <- 2
set.vertex.attribute(tn, attrname="Sex", value=list(SexTN), v=seq_len(network.size(tn)))
VertexSizeTN <- rep(1, network.size(tn))
VertexSizeTN[tnseedindex] <- 2
# Plotting the network
tnplot <- plot(tn,
                label = "",
                vertex.col = SexTN,
                vertex.cex = VertexSizeTN,
                jitter=T)
# Saving the plot
png(file="TN.png", height = 1800, width = 1800, res=300)
plot(tn,
     label = "",
     vertex.col = SexTN,
     vertex.cex = VertexSizeTN,
     jitter=T)

dev.off()

}


# Overlaying cn with tn
# matchingindicesTN <- pmatch(tnvertexnames, cnvertexnames)
# tncoord <- cnplot[matchingindicesTN, ]
# png(file="TN_CN.png", height = 1800, width = 1800, res=300)
# plot(tn, label = "", vertex.col = SexTN, coord = tncoord, vertex.cex = VertexSizeTN)
# dev.off()

######
# Last step: plot phylogeny of transmission tree
######

# pg <- phylo4(nj(geodist(tn)$gdist))
# apepg <- as(pg, "phylo")
# apepgrooted <- multi2di(apepg, random = FALSE)
# apepgrooted$edge.length <- c(30-10.18568, # 7-6
#                              15.05207-10.18568, # 7-8
#                              30-15.05207, # 8-5
#                              16.13587-15.05207, # 8-9
#                              22.99193-16.13587, # 9-10
#                              24.29339-22.99193, # 10-11
#                              30-24.29339, # 11-2
#                              30-24.29339, # 11-1
#                              30-22.99193, # 10-3
#                              30-16.13587  #  9-4
# )
# 
# apepgrooted$root.edge <- 10.18568-10
# apepgrooted$tip.label <- c("Woman 46",
#                            "Woman 37",
#                            "Man 23",
#                            "Man 4",
#                            "Man 2",
#                            "Woman 28")
# png(file="PG.png", height = 1800, width = 1800, res=300)
# plot(apepgrooted, root.edge = TRUE)
# dev.off()

# 
# 
#   g <- network(edata[ , 1:2], directed = TRUE, matrix.type = "edgelist")
# vertexnames <- network.vertex.names(g)
# plot(network(edata, directed=T, matrix.type = "edgelist"), label = vertexnames)
# 
# plot(g, label = vertexnames)
# 
#   g2 <- network(edataplusroot[ , 1:2], directed = TRUE, matrix.type = "edgelist")
# 
# 
# test <- phylo4(nj(geodist(g)$gdist))
# apetest <- as(test, "phylo")
# apetestrooted <- multi2di(apetest, random = FALSE)
# plot(apetestrooted)
# 
# phylo4(nj(geodist(g2)$gdist))


  # Now we add node attributes
  # Now we add edge attributes

  # Now we extract essential information from edata to create phylo tree.
  # Constraints to input for phylo4 tree constructor:
  
  #   if the tree has edge lengths defined, the number of edge lengths must match
  #   the number of edges; the number of tip labels must match the number of tips;
  #   in a tree with ntips tips and nnodes (total) nodes, nodes 1 to ntips must be
  #   tips if the tree is rooted, the root must be node number ntips+1 and the
  #   root node must be the first row of the edge matrix; tip labels, node labels,
  #   edge labels, edge lengths must have proper internal names (i.e. internal
  #   names that match the node numbers they document); tip and node labels must be
  #   unique.

  # First essential element is x, the matrix of edges.
#   nodes <- 1:(ntips+ninternalnodes)
#   rootnodenumber <- ntips+1
#   internalnodenumbers <- (rootnodenumber+1):max(nodes)
#   tipnumbers <- 1:ntips
#   
#   firstrow <- c(rootnodenumber, 1)
  


################
# Create phylogeny dataset
################


#   nrels[i] <- nrow(r)
#   cumulincid[i] <- sum(p$InfectTime>=60 & p$InfectTime<Inf) #cfg["hivseed.time"] & p$InfectTime<Inf)
#   # How many infections to place in the last 10 years of the simulation? (simulation years 50-60)
#   
#   allrpdata <- rbind(allrpdata, rp)
# } 
# 
# summarystats <- (data.frame(simID, meanagegap, varagegap, BetweenVar, WithinVar,
#                             nrels, cumulincid))
# summarystats
# }
  
  



# hivincdata <- data.frame(timesteps,incidence)
# ggplot(hivincdata,aes(x=timesteps, y=1000*incidence)) +
#   geom_line(colour = "darkred", size=2) +
#   guides(colour = FALSE) +
#   xlab("Simulation time") +
#   ylab("HIV incidence (per 100 PY)")



incidentcase <- function(DT, interval){ # arguments are the eventlog data.table and a vector of start and end times of the interval
  DTincident <- DT[eventname == "transmission" & eventtime > interval[1] & eventtime <= interval[2]]
}

nonHIVdeath <- function(DT, interval){ # arguments are the eventlog data.table and a vector of start and end times of the interval
  DTnonHIVdeath <- DT[eventname == "normalmortality" & eventtime > interval[1] & eventtime <= interval[2]]
}


alive.infected <- function(DT, time, dec.places){ # arguments are the personlog data.table and a point in time
  DTalive <- DT[TOB <= time & TOD > time]
  DTalive$Age <- time -round(ceiling(DTalive$TOB), dec.places) # Next we allocate them in discrete age bins with bin size as wide as timestep
  DTalive$Infected <- time >= DTalive$InfectTime # Now we allocate infection status to all people in our table of living people
  return(DTalive)
}

decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

# ####
# # To create heatplot, we need to loop through time steps, and calculate age- and gender-specific prevalence and incidence rates.
# alltimestepsdata <- data.table()
# timesteps <- seq(as.numeric(cfg["hivseed.time"]), as.numeric(cfg["population.simtime"]), by=1)
# for (time_i in timesteps) {
#   alivetable <- alivetab(time_i)   # First we only take the data of people who were alive at time_i
#   alivetable$Age.1 <- time_i -round(alivetable$TOB, 0) # Next we allocate them in discrete age bins with bin size equal as width of timesteps
#   alivetable$Infected <- time_i >= alivetable$InfectTime # Now we allocate infection status to all people in our table of living people
#   alive_data.table <- data.table(alivetable) # We turn the dataset into a proper data.table
#   timestepdata <- alive_data.table[,sum(Infected) / nrow(alive_data.table),by="Gender,Age.1"] # And we calculate HIV prevalence, by gender and age group
#   setnames(timestepdata, "V1", "Prevalence")
#   timestepdata <- cbind(time_i, timestepdata)
#   alltimestepsdata <- rbind(alltimestepsdata, timestepdata)
# }
# alltimestepsdata <- cbind(i, alltimestepsdata) # i is the ID number of the simulation
# 
# 
# 
# # Test sqldf query
# # x is evaltime (e.g. at end of simulation)
# relsandconcurM <- function(x){
#   query <- paste0("SELECT IDm, rels,
#                   CASE WHEN rels > 1 THEN 1 ELSE 0 END concur
#                   FROM (
#                   SELECT IDm, COUNT(*) rels
#                   FROM r
#                   WHERE ", x, "> FormTime AND ", x, "< DisTime
#                   GROUP BY IDm
#                   )"
#                   )
#   return (sqldf(query))
# }
# relsandconcurW <- function(x){
#   query <- paste0("SELECT IDw, rels,
#                   CASE WHEN rels > 1 THEN 1 ELSE 0 END concur
#                   FROM (
#                   SELECT IDw, COUNT(*) rels
#                   FROM r
#                   WHERE ", x, "> FormTime AND ", x, "< DisTime
#                   GROUP BY IDw
#                   )"
#                   )
#   return (sqldf(query))
# }
# # selecting people who are alive at time x
# alive <- function(x){
#   query <- paste0("SELECT ID, TOB
#                   FROM p
#                   WHERE ", x, ">= TOB AND ", x, "< TOD
#                   "#                  GROUP BY ID"
#   )
#   return (sqldf(query))
# }
# # number of infections caused
# infcaused <- function(){
#   query <- paste0("SELECT InfectOrigID, COUNT(*) infcaused
#                   FROM p
#                   GROUP BY InfectOrigID"
#   )
#   return (sqldf(query))
# }
# # infected at evaltime x
# infected <- function(x){
#   query <- paste0("SELECT ID, log10SPVL
#                   FROM p
#                   WHERE ", x, ">= InfectTime"
#   )
#   return (sqldf(query))
# }
# 
# # chronically infected at evaltime x (past the acute phase already)
# chroninfected <- function(x){
#   query <- paste0("SELECT ID, log10SPVL
#                   FROM p
#                   WHERE ", x, ">= InfectTime + 0.25"
#   )
#   return (sqldf(query))
# }
# 
# # dead and infected at evaltime x
# deadandinfected <- function(x){
#   query <- paste0("SELECT *
#                   FROM p
#                   WHERE ", x, ">= InfectTime AND ", x, "> TOD"
#   )
#   return (sqldf(query))
# }
# 
# 
# # dead and infected with number of infections at evaltime x
# secondaryinfections <- function(completedinfperiod, infectionscaused){
#   query <- paste0("SELECT *
#                   FROM completedinfperiod
#                   LEFT OUTER JOIN infectionscaused
#                   ON completedinfperiod.ID = infectionscaused.InfectOrigID"
#   )
#   return (sqldf(query))
# }
# 
# # alive and infected at evaltime x
# aliveandinfected <- function(infected, alive){
#   query <- paste0("SELECT ID, log10SPVL
#                   FROM infected
#                   JOIN alive
#                   USING (ID)"
#   )
#   return (sqldf(query))
# }
# 
# # alive and chronically infected at evaltime x
# aliveandchroninfected <- function(chroninfected, alive){
#   query <- paste0("SELECT ID, log10SPVL
#                   FROM chroninfected
#                   JOIN alive
#                   USING (ID)"
#   )
#   return (sqldf(query))
# }
# 
# # number of relationships per person in pop
# relsandpop <- function(relsMandW, pop){
#   query <-paste0("SELECT *
#                  FROM pop
#                  LEFT OUTER JOIN relsMandW
#                  USING (ID)"
#   )
#   return (sqldf(query))
# }
# 
# # Adding TOB to r, so that we can plot age mixing pattern
# rwithIDmdata <- function(r, p){
#   query <-paste0("SELECT *
#                  FROM r
#                  JOIN p
#                  ON r.IDm = p.ID"
#   )
#   return (sqldf(query))
# }
# 
# # Adding TOB to r, so that we can plot age mixing pattern
# pGender <- function(p, vIDs){
#   query <-paste0("SELECT ID, Gender
#                  FROM p
#                  JOIN e
#                  ON p.ID = vIDs.ID"
#   )
#   return (sqldf(query))
# }
# 
# alive <- function(x){
#   query <- paste0("SELECT ID, TOB
#                   FROM p
#                   WHERE ", x, ">= TOB AND ", x, "< TOD
#                   "#                  GROUP BY ID"
#   )
#   return (sqldf(query))
# }
# 
# alivetab <- function(x){
#   query <- paste0("SELECT *
#                   FROM DT
#                   WHERE ", x, ">= TOB AND ", x, "< TOD
# "#                  GROUP BY ID"
#   )
#   return (sqldf(query))
# }
