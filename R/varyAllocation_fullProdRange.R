#******************************************************************************
# varyAllocation_VaryFixedER_fullProdRange.R
# Date revised: May 4, 2019
# Explainer: Runs closed loop simulation model to examine changes in performance
# across different  allocations to single vs. mixed stock fisheries and a wide 
# range of productivity levels.
#******************************************************************************

# Parameterization notes:

# -Multistock HCR is generic precautionary approach
# -Single stock HCR is retrospective (i.e. generational median abundance used to 
# open/close fishery)
# -Spawner observation error (SD = 0.2) impacts the accuracy of these estimates
# -Catch observation error present, but no impact unless estimated SR BMs used
# -Simulation includes moderate covariance (0.4), autocorrelation (0.2), and OU
# (0.07)
# -En route mortality included
# -Reference productivtiy is based on most recent estimates from Kalman filter 
# models then look at +/- 35% by 5% intervals

listOfPackages <- c("here", "parallel", "doParallel", "foreach", 
                    "tidyverse", "tictoc", "samSim")
newPackages <- listOfPackages[!(listOfPackages %in% 
                                  installed.packages()[ , "Package"])]
if (length(newPackages)) {
  install.packages(newPackages)
}
lapply(listOfPackages, require, character.only = TRUE)

simPar <- read.csv(here("data", "manProcScenarios",
                        "fraserMPInputs_varyAllocationProdRange.csv"),
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
nTrials <- 100

## Expand dataset to include full range of productivity regimes 
prodScalars <- seq(0.5, 1.1, by = 0.05)
synchScalars <- seq(0.1, 0.7, by = 0.05)
addScalar <- function(prodScalarIn, synchScalarIn = NULL, dat) {
  dum <- dat %>%
    filter(scenario == "varyProd") %>% 
    mutate(prodScalar = prodScalarIn)
  if (!is.null(synchScalarIn)) {
    dum2 <- dat %>% 
      filter(scenario == "varySynch") %>% 
      mutate(correlCU = synchScalarIn)
    dum <- rbind(dum, dum2)
  }
  return(dum)
}
simParNew <- lapply(seq_along(prodScalars), function(x) 
  addScalar(prodScalarIn = prodScalars[x], synchScalarIn = synchScalars[x],
            dat = simPar)) %>%
  do.call(rbind, .) %>% 
  mutate(nameMP = paste(singleHCR, propMixHigh, "Mix", "_", harvContRule, 
                        sep = ""),
         nameOM = case_when(
           scenario == "varyProd" ~ paste(prodScalar, "Prod", sep = ""),
           scenario == "varySynch" ~ paste(correlCU, "Synch", sep = ""))) 

dirNames <- unique(simParNew$nameOM)
# dirNames <- simParNew %>% 
#   filter(scenario == "varySynch") %>% 
#   select(nameOM) %>% 
#   unique() %>%
#   sapply(as.character) %>% 
#   as.vector()

# recoverySim(simParNew[12, ], cuPar, catchDat = catchDat, srDat = srDat,
#             variableCU = FALSE, ricPars, larkPars = larkPars, tamFRP = tamFRP,
#             dirName = "test", nTrials = 5, makeSubDirs = FALSE, random = FALSE)

for (i in seq_along(dirNames)) {
  dirName <- dirNames[i]
  d <- subset(simParNew, nameOM == dirNames[i])
  simsToRun <- split(d, seq(nrow(d)))
  Ncores <- detectCores()
  cl <- makeCluster(Ncores - 1) #save two cores
  registerDoParallel(cl)
  clusterEvalQ(cl, c(library(MASS),
                     library(here),
                     library(sensitivity),
                     library(mvtnorm),
                     library(scales), #shaded colors for figs
                     library(viridis), #color blind gradient palette
                     library(gsl),
                     library(dplyr),
                     library(Rcpp),
                     library(RcppArmadillo),
                     library(sn),
                     library(samSim)))
  #export custom function and objects
  clusterExport(cl, c("simsToRun", "recoverySim", "cuPar", "dirName", "nTrials",
                      "catchDat", "srDat", "ricPars", "dirName", "larkPars",
                      "tamFRP"), envir = environment())
  tic("run in parallel")
  parLapply(cl, simsToRun, function(x) {
    recoverySim(x, cuPar, catchDat = catchDat, srDat = srDat,
                variableCU = FALSE, ricPars, larkPars = larkPars,
                tamFRP = tamFRP, cuCustomCorrMat = NULL, dirName = dirName,
                nTrials = nTrials, makeSubDirs = FALSE, random = FALSE)
  })
  stopCluster(cl) #end cluster
  toc()
}


## Prep data
agVarsToPlot <-  c("medRecRY", "medSpawners", "ppnCULower", "medCatch", 
                   "stabilityCatch", "ppnCUExtant", "ppnCUUpper")

agDat <- buildDataAgg(dirNames, agVars =  agVarsToPlot, 
                      keyVarName = "ppnMixed") %>% 
  mutate(mp =  "genPA",
         scenario = str_extract(om, "[A-Z]+[a-z]+"),
         ppnMixedNew = as.numeric(as.character(ppnMixed))) %>% 
  mutate(prod = case_when(
           scenario == "Prod" ~ str_extract(om, "\\d+\\.*\\d*"),
           scenario == "Synch" ~ "1"), 
         synch = case_when(
           scenario == "Synch" ~ str_extract(om, "\\d+\\.*\\d*"),
           scenario == "Prod" ~ "0.4")) %>% 
  mutate(prod = as.numeric(prod),
         synch = as.numeric(synch))

         
## Make heat map function
plotHeat <- function(plotVar = "medRecRY", legName = "rec", bins = 20,
                     dat = agDat) {
  dum <- dat %>%
    filter(var == plotVar) %>% 
    mutate(xAxis = case_when(
      scenario == "Synch" ~ synch,
      scenario == "Prod" ~ prod,
    ))
  xLab <- ifelse(dum$scenario == "Synch", "Pairwise Correlation", 
                 "Productivity Scalar")
  dumInterp <- akima::interp(x = dum$xAxis, y = dum$ppnMixedNew, 
                             z = dum$avg, nx = bins, ny = bins) %>% 
    akima::interp2xyz() %>% 
    as.data.frame()
  names(dumInterp) <- c("opMod", "ppnMixedNew", "value")
  ggplot(dumInterp, aes(x = opMod, y = ppnMixedNew, fill = value)) + 
    geom_tile() + 
    viridis::scale_fill_viridis(name = legName, option = "viridis") +
    # scale_fill_gradient(low = "red", high = "blue", name = legName) +
    labs(x = xLab, y = "Proportion Mixed Stock") +
    theme_sleekX()
}
synchDat <- agDat %>% 
  filter(scenario == "Synch")
prodDat <- agDat %>% 
  filter(scenario == "Prod")

pdf(here("figs", "heatMaps", "recHeat.pdf"), height = 4, width = 5)
plotHeat(plotVar = "medRecRY", legName = "Mean\nRecruits", dat = prodDat)
plotHeat(plotVar = "medRecRY", legName = "Mean\nRecruits", dat = synchDat)
dev.off()

pdf(here("figs", "heatMaps", "ppnHeat.pdf"), height = 4, width = 5)
plotHeat(plotVar = "ppnCULower", legName = "Ppn CU\n Above\nLower BM", dat = prodDat)
plotHeat(plotVar = "ppnCULower", legName = "Ppn CU\n Above\nLower BM", dat = synchDat)
dev.off()

pdf(here("figs", "heatMaps", "catchStblHeat.pdf"), height = 4, width = 5)
plotHeat(plotVar = "stabilityCatch", legName = "Catch\nStability", 
         dat = prodDat)
plotHeat(plotVar = "stabilityCatch", legName = "Catch\nStability", 
         dat = synchDat)
dev.off()

#Presentation figs
png(here("figs", "heatMaps", "recHeat.png"), height = 4, width = 5, 
    units = "in", res = 400)
plotHeat(plotVar = "medRecRY", legName = "Mean\nRecruits", dat = prodDat)
dev.off()

png(here("figs", "heatMaps", "catchStblHeat.png"), height = 4, width = 5, 
    units = "in", res = 400)
plotHeat(plotVar = "stabilityCatch", legName = "Catch\nStability", 
         dat = synchDat)
dev.off()



## Trade off plot
consVec <- rep(c("medRecRY", "ppnCULower"), each = 1)
catchVec <- rep(c("medCatch"), times = 2)

synchDat <- agDat %>% 
  filter(scenario == "Synch")

plotAgTradeoff(synchDat, consVar = "medRecRY", catchVar = "medCatch", 
               facet = "om", shape = "mp", showUncertainty = TRUE, 
               xLab = "medCatch", yLab = "medRecRY", 
               legendLab = "Proportion\nTAC Mix\nFishery",
               axisSize = 11, dotSize = 4, legendSize = 13, lineSize = 0.5, 
               scaleAxis = "fixed")

prodDat <- agDat %>% 
  filter(scenario == "Prod")

plotAgTradeoff(prodDat, consVar = "medRecRY", catchVar = "medCatch", 
               facet = "om", shape = "mp", showUncertainty = TRUE, 
               xLab = "medCatch", yLab = "medRecRY", 
               legendLab = "Proportion\nTAC Mix\nFishery",
               axisSize = 11, dotSize = 4, legendSize = 13, lineSize = 0.5, 
               scaleAxis = "fixed")
