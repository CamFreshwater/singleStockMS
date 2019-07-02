#******************************************************************************
# runModelVaryMP_fixedComparison.R
# Date revised: Sep. 28, 2018
# Inputs: recoverySim.R
# Outputs: pdf plots, matrices for plotting
# Explainer: Runs closed loop simulation model with different OMs and MPs based
# on input csv files; structured so that MPs vary within a scenario, not OMs. 
# This script runs a "vanilla" comaprison between a range of fixed exploitation
# rates and the pseudo TAM rule.
#******************************************************************************


# Check if required packages are installed and run
listOfPackages <- c("plyr", "here", "parallel", "doParallel", "foreach", 
                    "reshape2", "tidyr", "gsl", "tictoc", "stringr", "dplyr",
                    "synchrony", "zoo", "Rcpp", "RcppArmadillo", "sn", 
                    "sensitivity", "mvtnorm", "devtools")
newPackages <- listOfPackages[!(listOfPackages %in% 
                                  installed.packages()[ , "Package"])]
if (length(newPackages)) {
  install.packages(newPackages)
}
lapply(listOfPackages, require, character.only = TRUE)

here <- here::here()

source(here("scripts/func/postProcessing.R"))
source(here("scripts/func/simUtilityFunc.R"))
source(here("scripts/recoverySim.R"))


simPar <- read.csv(here("data/manProcScenarios/fraserMPInputs_fixedComp.csv"), 
                   stringsAsFactors = F)
cuPar <- read.csv(here("data/fraserCUpars.csv"), stringsAsFactors=F)
srDat <- read.csv(here("data/fraserRecDatTrim.csv"), stringsAsFactors=F)
catchDat <- read.csv(here("data/fraserCatchDatTrim.csv"), stringsAsFactors=F)
ricPars <- read.csv(here("data/fraserDat/rickerMCMCPars.csv"), stringsAsFactors=F)
larkPars <- read.csv(here("data/fraserDat/larkinMCMCPars.csv"), stringsAsFactors=F)
tamFRP <- read.csv(here("data/fraserDat/tamRefPts.csv"), stringsAsFactors=F)


### SET UP MODEL RUN ----------------------------------------------------------

## Define simulations to be run
nTrials <- 100

simParTrim <- subset(simPar, scenario == "tamComp")
scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species), 
                                                sep = "_"))
# recoverySim(simParTrim[1,], cuPar, catchDat = catchDat, srDat = srDat,
#             variableCU = FALSE, ricPars, larkPars = larkPars, tamFRP = tamFRP,
#             cuCustomCorrMat = cuCustomCorrMat, dirName = "Test", nTrials = 2,
#             multipleMPs = FALSE)
# simsToRun <- split(simParTrim, seq(nrow(simParTrim)))
# for(i in seq_along(simsToRun)) {
#   recoverySim(simsToRun[[i]], cuPar, catchDat=catchDat, srDat=srDat, variableCU=FALSE,
#               ricPars, larkPars=larkPars, tamFRP=tamFRP, dirName=dirNames[1],
#               nTrials = nTrials, multipleMPs=TRUE)
# }

for (i in seq_along(dirNames)) {
  dirName <- dirNames[i]
  d <- subset(simParTrim, scenario == scenNames[i])
  simsToRun <- split(d, seq(nrow(d)))
  Ncores <- detectCores()
  cl <- makeCluster(Ncores - 2) #save two cores
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
                     library(sn)))
  if(simsToRun[[1]]$species == "sockeye"){ #export custom function and objects
    clusterExport(cl, c("simsToRun","recoverySim","cuPar","dirName","nTrials",
                        "catchDat","srDat","ricPars","dirName","larkPars",
                        "tamFRP"), envir=environment()) 
    tic("run in parallel")
    parLapply(cl, simsToRun, function(x) {
      recoverySim(x, cuPar, catchDat=catchDat, srDat=srDat, variableCU=FALSE, 
                  ricPars, larkPars=larkPars, tamFRP=tamFRP, dirName=dirName, 
                  nTrials=nTrials, multipleMPs=TRUE)
    })
    stopCluster(cl) #end cluster
    toc()
  }
}


### GENERATE OUTPUT FIGS ------------------------------------------------------
## Read in data from directories and convert CU-lists into dataframes then plot
cuDat <- buildDataCU(dirNames = dirNames, 
                      cuVars = c("medSpawners", "medCatch", "ppnYrsUpper", 
                                 "ppnYrsLower"),
                      keyVarName = "expRate", 
                      selectedCUs = c("E.St", "Bwrn", "Chlk", "L.Sh", "Clts", 
                                      "Hrrs")) %>% 
  mutate(muName = recode(muName, EStu = "E. Stuart", ESumm = "E. Summer",
                         Summ = "Summer", Lat = "Late"))
cuDat$muName <- factor(cuDat$muName, levels(cuDat$muName)[c(1, 2, 4, 3)])
cuDat$om <- factor(cuDat$om, levels(cuDat$om)[c(2, 1, 3)])

# trade off plots
consVec <- c("medSpawners", "ppnYrsLower")
yLabVec <- c("Median Spawner Abundance", "Proportion Yrs Above Lower BM")

toPlotList <- lapply(seq_along(consVec), function (h) {
  p <- plotCUTradeoff(cuDat, consVar = consVec[h], catchVar = "medCatch", 
                      facet = "mu", panel = "om", showUncertainty = TRUE,
                      xLab = "Median Catch", yLab = yLabVec[h], 
                      legendLab = "Proportion Mixed\nStock TAC",
                      lineSize = 0.3, main = FALSE)
  p2 <- p[c(2, 1, 3)]
  return(p2)
})

pdf(here("outputs/varyHCR/Fig1_CUTradeoffs.pdf"), 
    height = 6.5, width = 8.5)
sapply(toPlotList, function (h) 
  sapply(h, function(i) print(i))
)
dev.off()

# dot plots
varVec <- c("medSpawners", "medCatch", "ppnYrsLower")
yLabVec <- c("Median Spawner Abundance", "Median Catch Abundance", 
             "Proportion Above Lower BM")
cuDatRef <- cuDat %>% 
  filter(om == "ref")
dotPlotList <- lapply(seq_along(varVec), function(x) {
  p <- plotCUDot(cuDatRef, plotVar = varVec[x], group = "om", 
                 legendLab = "Operating\nModel", 
                 xLab = "Exploitation Rate", 
                 yLab = yLabVec[x], axisSize = 14, dotSize = 2.5, 
                 lineSize = 0.9, legendSize = 14)
  return(p)
})

pdf(here("outputs/varyHCR/Fig2_CUDotplots.pdf"), 
    height = 6.5,
    width = 9.5)
sapply(dotPlotList, function(x) print(x))
dev.off()


## Read in data from directories and convert aggregate lists into df then plot
agDat <- buildDataAgg(dirNames, agVars =  c("medSpawners", "medCatch", 
                                             "ppnCULower", "ppnCUStable",
                                             "ppnYrsHighCatch", 
                                             "ppnFisheriesOpen"), 
                       keyVarName = "ppnMixed") %>% 
  mutate(om = recode(om, ref = "Reference", lowProd = "Low Productivity"))
agDat$om <- factor(agDat$om, levels(agDat$om)[c(2, 1)])
# agDat <- agDatF %>% 
#   filter(mp == "TAM", !plotOrder == "6")

# trade off plots
consVec <- c("medSpawners", "ppnCULower")
yLabVec <- c("Median Spawner Abundance", "Proportion CUs Above Upper BM")
toAgPlotList <- lapply(seq_along(consVec), function(x) {
  p <- plotAgTradeoff(agDat, consVar = consVec[x], catchVar = "medCatch",
                      facet = "om", showUncertainty = TRUE,
                      xLab = "Median Catch", yLab = yLabVec[x], 
                      legendLab = "Propotion\nTAC Mix\nFishery",
                      axisSize = 14, dotSize = 4, legendSize = 14,
                      lineSize = 0.5, freeY = FALSE)
  return(p)
})

pdf(here("outputs/varyHCR/Fig3_AgTradeoffs_2om.pdf"), 
    height = 4.5, width = 6.5)
sapply(toAgPlotList, function(x) print(x))
dev.off()

agDatRef <- agDat %>%
  filter(om == "ref")
pdf(here("outputs/varyHCR/Fig4_AgDotPlot.pdf"), 
    height = 4.5, width = 8.5)
plotAgDot(agDatRef, group = "om", legendLab = "Operating\nModel",
          xLab = "Proportion of TAC in Mixed Fishery", yLab = "", axisSize = 14,
          dotSize = 4, lineSize = 1, legendSize = 14)
dev.off()


### Radar plots
dat <- buildDataAgg(dirNames, agVars =  c("ppnCUExtant", "ppnYrsHighCatch", 
                                          "ppnCULower", "ppnFisheriesOpen",
                                          "ppnCUStable"), 
                    keyVarName = "expRate") %>% 
  filter(om == "ref", !expRate == "0.2" & !expRate == "0.6")
pdf(here("outputs/varyHCR/Fig5_AgRadarPlot.pdf"), 
    height = 6, width = 7.5)
plotRadar(dum, xLab = c("Extant", "Years\nHigh Catch", "CUs\nLower BM", 
                        "Fisheries\n Open", "CUs\nStable"),
          plotVars = NULL, groupingVar = NULL, cu = FALSE, 
          legendLab = "Exploitation Rate", axisSize = 11)
dev.off()

