#******************************************************************************
# varyAllocation_mainSimRun.R
# Date revised: May 10, 2019
# Explainer: Runs closed loop simulation model to examine impacts of divergent
# productivity.
#******************************************************************************

# This script is approximately equivalent varyAllocation_fullProdRange.R but
# does not include multiple productivity or synchrony treatments, only alternative
# en route mortality.

# Parameterization notes:

listOfPackages <- c("here", "parallel", "doParallel", "foreach", 
                    "tidyverse", "tictoc", "samSim")
newPackages <- listOfPackages[!(listOfPackages %in% 
                                  installed.packages()[ , "Package"])]
if (length(newPackages)) {
  install.packages(newPackages)
}
lapply(listOfPackages, require, character.only = TRUE)

simPar <- read.csv(here("data", "manProcScenarios",
                        "fraserMPInputs_varyAllocationVaryMixHCR.csv"),
                   stringsAsFactors = F)
cuPar <- read.csv(here("data/summOnlyCUPars.csv"), stringsAsFactors = F)
srDat <- read.csv(here("data/fraserRecDatTrim.csv"), stringsAsFactors = F)
catchDat <- read.csv(here("data/fraserCatchDatTrim.csv"), stringsAsFactors = F)
ricPars <- read.csv(here("data/trimRecursiveRickerMCMCPars.csv"),
                    stringsAsFactors = F)
larkPars <- read.csv(here("data/trimRecursiveLarkinMCMCPars.csv"),
                     stringsAsFactors = F)
tamFRP <- read.csv(here("data/tamRefPts.csv"), stringsAsFactors = F)

## Define simulations to be run
nTrials <- 250

## Make unique MP vector 
simPar$nameMP <- paste(simPar$propMixHigh, simPar$singleHCR, "_", 
                       simPar$harvContRule, sep = "")

simParTrim <- simPar %>%
  filter(nameOM %in% c("ref", "divProd", "lowProd"),
         benchmark == "stockRecruit",
         harvContRule == "genPA")
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
recoverySim(simParTrim[6, ], cuPar, catchDat = catchDat, srDat = srDat,
            variableCU = FALSE, ricPars, larkPars = larkPars, tamFRP = tamFRP,
            cuCustomCorrMat = cuCustomCorrMat, dirName = "Test", nTrials = 4,
            makeSubDirs = FALSE, random = FALSE)
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
                nTrials=nTrials, makeSubDirs=TRUE, random = FALSE)
  })
  stopCluster(cl) #end cluster
  toc()
}

#------------------------------------------------------------------------------

# Aggregate PMs to plot
agVarsToPlot <-  c("medRecRY", "medSpawners", "ppnCULower", "ppnCUExtant", 
                   # "ppnCUStable", "ppnCUExtant", "ppnYrsHighCatch", 
                   # "ppnMixedOpen", 
                   "medCatch", "stabilityCatch", "medER")

agDat <- buildDataAgg(dirNames, agVars =  agVarsToPlot, 
                      keyVarName = "ppnMixed") %>%
  mutate(mp = case_when(
    str_detect(mp, "genPA") ~ "Generic",
    str_detect(mp, "TAM") ~ "TAM"
  ))

consVec <- c("medRecRY", "ppnCULower", "ppnCUExtant")
conLabels <- c("Median Aggregate Return Abundance", 
               "Proportion of CUs Above Benchmark",
               "Proportion of CUs Extant")

toAgList <- lapply(seq_along(consVec), function(x) {
  p <- plotAgTradeoff(agDat, consVar = consVec[x], 
                      catchVar = "medCatch",
                      facet = NULL, shape = "om", showUncertainty = TRUE,
                      xLab = "Median Aggregate Catch", yLab = conLabels[x], 
                      legendLab = "Proportion\nTAC Mix\nFishery",
                      axisSize = 11, dotSize = 4, legendSize = 13,
                      lineSize = 0.5, scaleAxis = "fixed")
  return(p)
})
ggpubr::ggarrange(toAgList[[1]], toAgList[[2]], toAgList[[3]],
                  ncol = 1, nrow = 3, common.legend = TRUE, legend = "right")

