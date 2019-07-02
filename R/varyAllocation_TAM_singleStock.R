#******************************************************************************
# varyAllocation_TAM_singleStock.R
# Date revised: Mar. 27, 2019
# Explainer: Runs closed loop simulation model to examine changes in performance
# with TAM and different allocations to single vs. mixed stock fisheries. 
# Includes alternative productivity and en route mortality operating models.
#******************************************************************************

# Parameterization notes:

### NOTE RERUN AFTER REPLACING SR PARAMETERS WITH CONSTRAINED BETAS ###
### ALSO RECALCULATE MEDIAN DECLINES IN ALPHAS WITH CHANGED BETAS ###

# -Single stock HCR is retrospective (i.e. generational median abundance used to 
# open/close fishery)
# -Spawner observation error (SD = 0.2) impacts the accuracy of these estimates
# -Catch observation error present, but no impact unless estimated SR BMs used
# -Simulation includes moderate covariance (0.4), autocorrelation (0.2), and OU
# (0.07)
# -En route mortality included
# -Focus on low productivity scenarios primarily since reference alpha results 
# in overly rapid recovery; discuss with co-authors what best parameterization
# should be

listOfPackages <- c("here", "parallel", "doParallel", "foreach", 
                    "tidyverse", "tictoc", "samSim")
newPackages <- listOfPackages[!(listOfPackages %in% 
                                  installed.packages()[ , "Package"])]
if (length(newPackages)) {
  install.packages(newPackages)
}
lapply(listOfPackages, require, character.only = TRUE)

simPar <- read.csv(here("data", "manProcScenarios",
                        "fraserMPInputs_varyAllocationTAM.csv"),
                   stringsAsFactors = F)
cuPar <- read.csv(here("data/fraserCUpars.csv"), stringsAsFactors = F)
srDat <- read.csv(here("data/fraserRecDatTrim.csv"), stringsAsFactors = F)
catchDat <- read.csv(here("data/fraserCatchDatTrim.csv"), stringsAsFactors = F)
ricPars <- read.csv(here("data/pooledRickerMCMCPars.csv"), stringsAsFactors = F)
larkPars <- read.csv(here("data/pooledLarkinMCMCPars.csv"), 
                     stringsAsFactors = F)
tamFRP <- read.csv(here("data/tamRefPts.csv"), stringsAsFactors = F)

## Define simulations to be run
nTrials <- 250

## Make unique MP vector 
simPar$nameMP <- paste(simPar$propMixHigh, simPar$singleHCR, "_", 
                       simPar$harvContRule, sep = "")

simParTrim <- simPar %>%
  filter(!nameOM %in% c("noObsErr")) %>% 
  mutate(benchmark = "stockRecruit")
scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species), 
                                                sep = "_"))

# Focal CUs and CU-specific PMs
summCUs <- cuPar %>% 
  filter(manUnit == "Summ") %>% 
  mutate(abbStkName = abbreviate(stkName, minlength = 4)) %>% 
  select(abbStkName) %>% 
  unlist()

lateCUs <- cuPar %>% 
  filter(manUnit == "Lat") %>% 
  mutate(abbStkName = abbreviate(stkName, minlength = 4)) %>% 
  select(abbStkName) %>% 
  unlist()


#------------------------------------------------------------------------------

## Run simulation
# recoverySim(simParTrim[6,], cuPar, catchDat = catchDat, srDat = srDat,
#             variableCU = FALSE, ricPars, larkPars = larkPars, tamFRP = tamFRP,
#             cuCustomCorrMat = cuCustomCorrMat, dirName = "Test", nTrials = 4,
#             makeSubDirs = FALSE, random = FALSE)

# for (i in seq_along(dirNames)) {
#   dirName <- dirNames[i]
#   d <- subset(simParTrim, scenario == scenNames[i])
#   simsToRun <- split(d, seq(nrow(d)))
#   Ncores <- detectCores()
#   cl <- makeCluster(Ncores - 1) #save two cores
#   registerDoParallel(cl)
#   clusterEvalQ(cl, c(library(MASS),
#                      library(here),
#                      library(sensitivity),
#                      library(mvtnorm),
#                      library(scales), #shaded colors for figs
#                      library(here),
#                      library(synchrony),
#                      library(zoo),
#                      library(viridis), #color blind gradient palette
#                      library(ggplot2),
#                      library(gsl), #to calculate exact estimate of MS
#                      library(dplyr),
#                      library(Rcpp),
#                      library(RcppArmadillo),
#                      library(sn),
#                      library(samSim)))
#   clusterExport(cl, c("simsToRun","cuPar","dirName","nTrials",
#                       "catchDat","srDat","ricPars","larkPars",
#                       "tamFRP"), envir=environment())
#   tic("run in parallel")
#   parLapply(cl, simsToRun, function(x) {
#     recoverySim(x, cuPar, catchDat=catchDat, srDat=srDat, variableCU=FALSE,
#                 ricPars, larkPars=larkPars, tamFRP=tamFRP, dirName=dirName,
#                 nTrials=nTrials, makeSubDirs=TRUE, random = FALSE)
#   })
#   stopCluster(cl) #end cluster
#   toc()
# }

#------------------------------------------------------------------------------
# Aggregate PMs to plot
agVarsToPlot <-  c("medRecRY", "medSpawners", "ppnCULower", "ppnCUStable",
                   "ppnCUExtant",
                   "medCatch", "ppnYrsHighCatch", "stabilityCatch", 
                   "ppnMixedOpen", "medER")

simParTrim <- simPar %>% 
  filter(nameOM %in% c("lowProd", "ref", "fullEnRoute_lowProd"))
scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species), 
                                                sep = "_"))

multOMTamDat <- buildDataAgg(dirNames, agVars =  agVarsToPlot, 
                            keyVarName = "ppnMixed") %>%
  filter(om %in% c("ref", "lowProd", "fullEnRoute_lowProd")) %>% 
  mutate(om = fct_relevel(om, "ref", after = Inf)) %>% 
  mutate(mp = as.factor(mp),
         om = fct_recode(om, "Low Productivity" = "lowProd",
                         "En Route Mortality" = "fullEnRoute_lowProd",
                         "High Productivity" = "ref")) 

consVec <- rep(c("medSpawners", "ppnCULower", "medRecRY", "ppnCUStable"), 
               each = 1)
# catchVec <- rep(c("medCatch"), times = 2)

toAgPlotList <- lapply(seq_along(consVec), function(x) {
  p <- plotAgTradeoff(multOMTamDat, consVar = consVec[x], 
                      catchVar = "medCatch",
                      facet = "om", showUncertainty = TRUE,
                      yLab = consVec[x], xLab = "medCatch", 
                      legendLab = "Proportion\nTAC Mix\nFishery",
                      axisSize = 11, dotSize = 4, legendSize = 13,
                      lineSize = 0.5, scaleAxis = "fixed")
  return(p)
})
pdf(file = paste(here(),"/figs/tamTrends/aggTO_multOM.pdf", sep = ""),
    height = 4, width = 7)
sapply(toAgPlotList, function(x) print(x))
dev.off()
