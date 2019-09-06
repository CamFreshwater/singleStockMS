#******************************************************************************
# varyAllocation_singleStock_sens.R
# Date revised: June 10, 2019
# Explainer: Runs closed loop simulation model to examine changes in performance 
# across different sensitivity parameterizations. 
#******************************************************************************

# Parameterization notes:

# -Single stock HCR is retrospective (i.e. generational median abundance used to 
# open/close fishery)
# -Includes generic precautionary approach mixed-stock HCR
# -En route mortality included
# -Focus on low productivity scenarios primarily since that's the focus of
# primary analysis but adjust to match (may need to include both)

listOfPackages <- c("here", "parallel", "doParallel", "foreach", "tidyverse",
                    "tictoc", "samSim")
lapply(listOfPackages, require, character.only = TRUE)

simPar <- read.csv(here("data", "manProcScenarios",
                        "fraserMPInputs_varyAllocationSens.csv"), 
                   stringsAsFactors = F)
srDat <- read.csv(here("data/fraserRecDatTrim.csv"), stringsAsFactors = F)
catchDat <- read.csv(here("data/fraserCatchDatTrim.csv"), stringsAsFactors = F)
cuPar <- read.csv(here("data/summOnlyCUPars.csv"), stringsAsFactors = F)
ricPars <- read.csv(here("data/trimRecursiveRickerMCMCPars.csv"),
                    stringsAsFactors = F)
larkPars <- read.csv(here("data/trimRecursiveLarkinMCMCPars.csv"),
                     stringsAsFactors = F)
tamFRP <- read.csv(here("data/tamRefPts.csv"), stringsAsFactors = F)

## Define simulations to be run
nTrials <- 1300

## Make unique MP vector 
simPar$nameMP <- paste(simPar$propMixHigh, simPar$singleHCR, "_", 
                       simPar$harvContRule, sep = "")

simParTrim <- simPar 
# %>%
#   filter(scenario %in% c("ageTau"))

scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species), 
                                                sep = "_"))

# Focal CUs and CU-specific PMs
summCUs <- cuPar %>% 
  filter(manUnit == "Summ") %>% 
  mutate(abbStkName = abbreviate(stkName, minlength = 4)) %>% 
  select(abbStkName) %>% 
  unlist()

#------------------------------------------------------------------------------

## Run simulation
# recoverySim(simParTrim[1,], cuPar, catchDat = catchDat, srDat = srDat,
#             variableCU = FALSE, ricPars, larkPars = larkPars, tamFRP = tamFRP,
#             cuCustomCorrMat = cuCustomCorrMat, dirName = "Test", nTrials = 10,
#             makeSubDirs = FALSE, random = FALSE)

for (i in seq_along(dirNames)) {
  dirName <- dirNames[i]
  d <- subset(simParTrim, scenario == scenNames[i])
  simsToRun <- split(d, seq(nrow(d)))
  Ncores <- detectCores()
  cl <- makeCluster(Ncores - 1) #save two cores
  registerDoParallel(cl)
  clusterEvalQ(cl, c(library(MASS),
                     library(here),
                     library(sensitivity),
                     library(mvtnorm),
                     library(scales), #shaded colors for figs
                     library(here),
                     library(synchrony),
                     library(zoo),
                     library(viridis), #color blind gradient palette
                     library(ggplot2),
                     library(gsl), #to calculate exact estimate of MS
                     library(dplyr),
                     library(Rcpp),
                     library(RcppArmadillo),
                     library(sn),
                     library(samSim)))
  clusterExport(cl, c("simsToRun","cuPar","dirName","nTrials",
                      "catchDat","srDat","ricPars","larkPars",
                      "tamFRP"), envir=environment())
  tic("run in parallel")
  parLapply(cl, simsToRun, function(x) {
    recoverySim(x, cuPar, catchDat=catchDat, srDat=srDat, variableCU=FALSE,
                ricPars, larkPars=larkPars, tamFRP=tamFRP, dirName=dirName,
                nTrials=nTrials, makeSubDirs = FALSE, random = FALSE)
  })
  stopCluster(cl) #end cluster
  toc()
}

## Clean data
vars <- c("medRecRY", "medSpawners", "ppnCULower", "medCatch", "stabilityCatch",
          "medER")

pdat = NULL
for (h in seq_along(dirNames)) {
  agList <- genOutputList(dirNames[h], agg = TRUE)
  singleScen = NULL
  for (i in seq_along(vars)) {
    dum <- data.frame(scen = as.factor(rep(scenNames[h], 
                                           length.out =  length(agList))),
                      var = rep(vars[i], length.out = length(agList)),
                      om = as.factor(sapply(agList, function(x)
                        unique(x$opMod))),
                      mp = as.factor(sapply(agList, function(x)
                        unique(x$manProc))),
                      avg = sapply(agList, function(x) median(x[,vars[i]])),
                      lowQ = sapply(agList, function(x) qLow(x[,vars[i]])),
                      highQ = sapply(agList, function(x) qHigh(x[,vars[i]])),
                      row.names = NULL
    )
    singleScen <- rbind(singleScen, dum)
  }
  #merge multiple scenarios into one dataframe
  pdat <- rbind(pdat, singleScen) 
}

plotDat <- pdat %>% 
  mutate(labelledVar = case_when(
    var == "medRecRY" ~ "Recruits",
    var == "medSpawners" ~ "Spawners",
    var == "ppnCULower" ~ "Ppn. Stocks Above BM",
    var == "medCatch" ~ "Catch",
    var == "stabilityCatch" ~ "Catch Stability",
    var == "medER" ~ "Exp. Rate"),
    pmClass = case_when(
      var %in% c("medRecRY", "medSpawners", "ppnCULower") ~ "cons",
      TRUE ~ "catch"
    ))

# write.csv(plotDat, here("outputs", "generatedData",
#                         "sensAnalysisSummaryData.rds"),
#           row.names = FALSE)

sensPlot <- function(Var) {
  dum <- plotDat %>% 
    filter(var == Var,
           # mp = MP
           !scen == "genericRetro") 
  dumRef <- plotDat %>% 
    filter(var == Var,
           # mp = MP
           scen == "genericRetro")
  
  p <- ggplot(dum, aes(x = scen, y = avg, ymin = lowQ, ymax = highQ,
                       fill = om)) +
    labs(x = "", y = unique(dum$labelledVar)) +
    geom_pointrange(shape = 21, fatten = 4, size = 1,
                    position = position_dodge(width = 0.65)) +
    geom_hline(dumRef, mapping = aes(yintercept = avg), linetype = 1) +
    geom_hline(dumRef, mapping = aes(yintercept = lowQ), linetype = 2, 
               size = 0.5) +
    geom_hline(dumRef, mapping = aes(yintercept = highQ), linetype = 2, 
               size = 0.5) +
    scale_x_discrete(labels = c("ageTau" = "Mat.\nAge", "mixOU" = 
                                  "Mixed\nOU", 
                                "singOU" = "Single\nOU",
                                "rho" = "Temp.\nAuto.", 
                                "enRouteSig" = "Mig.\nMort.")) + 
    scale_fill_discrete(name = "Operating\nModel",
                        labels = c("High\nVariability", "Low\nVariability")) +
    facet_wrap(~mp)
  return(p)
}

outList <- lapply(seq_along(vars), function(x) sensPlot(Var = vars[x]))
outList

## allocatiosn don't dramatically change interpretation of sensitivity analysis
## present 50/50 with justification in text