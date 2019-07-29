#******************************************************************************
# varyAllocation_VaryFixedER_singleStock.R
# Date revised: Mar. 23, 2019
# Explainer: Runs closed loop simulation model to examine changes in performance
# across different fixed exploitation rates and allocations to single vs. mixed
# stock fisheries. Includes alternative productivity and en route mortality 
# operating models.
#******************************************************************************

# Parameterization notes:

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
                        "fraserMPInputs_varyAllocationVaryFixedER.csv"),
                   stringsAsFactors = F)
cuPar <- read.csv(here("data/summOnlyCUPars.csv"), stringsAsFactors = F)
srDat <- read.csv(here("data/fraserRecDatTrim.csv"), stringsAsFactors = F)
catchDat <- read.csv(here("data/fraserCatchDatTrim.csv"), stringsAsFactors = F)
# ricPars <- read.csv(here("data/constrainedRickerMCMCPars.csv"), stringsAsFactors = F)
# larkPars <- read.csv(here("data/constrainedLarkinMCMCPars.csv"), 
#                      stringsAsFactors = F)
ricPars <- read.csv(here("data/trimRecursiveRickerMCMCPars.csv"),
                    stringsAsFactors = F)
larkPars <- read.csv(here("data/trimRecursiveLarkinMCMCPars.csv"),
                     stringsAsFactors = F)
tamFRP <- read.csv(here("data/tamRefPts.csv"), stringsAsFactors = F)

## Define simulations to be run
nTrials <- 500

## Make unique MP vector 
simPar$nameMP <- paste(simPar$propMixHigh, simPar$singleHCR, "_", simPar$canER, 
                       simPar$harvContRule, sep = "")

simParTrim <- simPar %>%
  filter(nameOM %in% c("lowProd", "ref", "fullEnRoute_lowProd"))
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
# recoverySim(simParTrim[68,], cuPar, catchDat = catchDat, srDat = srDat,
#             variableCU = FALSE, ricPars, larkPars = larkPars, tamFRP = tamFRP,
#             cuCustomCorrMat = cuCustomCorrMat, dirName = "Test", nTrials = 4,
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
                nTrials=nTrials, makeSubDirs=TRUE, random = FALSE)
    })
  stopCluster(cl) #end cluster
  toc()
}

#------------------------------------------------------------------------------
# Aggregate PMs to plot
agVarsToPlot <-  c("medRecRY", "medSpawners", "ppnCULower", "ppnCUStable",
                   "ppnCUExtant",
                   "medCatch", "ppnYrsHighCatch", "stabilityCatch", 
                   "ppnMixedOpen", "medER")

## Make line plots demonstrating multiple fixed ERs and allocations
lpAgDat <- buildDataAgg(dirNames, agVars =  agVarsToPlot, 
                        keyVarName = "ppnMixed") %>%
  #pool MPs by fixedER
  mutate(mp = as.numeric(case_when(
    str_detect(mp, "0.1fixedER") ~ "0.1",
    str_detect(mp, "0.3fixedER") ~ "0.3",
    str_detect(mp, "0.5fixedER") ~ "0.5",
    str_detect(mp, "0.7fixedER") ~ "0.7",
    str_detect(mp, "0.9fixedER") ~ "0.9"
  ))) %>% 
  filter(om == "lowProd",
         mp %in% c("0.1", "0.3", "0.5", "0.7", "0.9")) %>%
  mutate(labelledVar = var) %>%
  mutate(labelledVar = fct_recode(labelledVar, 
                                  "Aggregate\nReturn Abundance" = "medRecRY",
                                  "Aggregate\nSpawner Abundance" = 
                                    "medSpawners",
                                  "Ppn. MUs Above\nEsc. Goal" = "ppnMixedOpen",
                                  "Ppn. CUs Not\nRed (Throughout)" = 
                                    "ppnCULower",
                                  "Ppn. CUs Not\nRed (Sim. End)" = 
                                    "ppnCUStable",
                                  "Ppn. CUs Extant" = 
                                    "ppnCUExtant",
                                  "Aggregate Catch" = "medCatch",
                                  "Ppn. Yrs. Above\nCatch Target" = 
                                    "ppnYrsHighCatch",
                                  "Catch Stability" = "stabilityCatch"))

colPal <- viridis::viridis(length(levels(lpAgDat$ppnMixed)), begin = 0, 
                           end = 1)
names(colPal) <- levels(lpAgDat$ppnMixed)

ggplot(lpAgDat %>% filter(!var == "medER"), 
       aes(x = mp, y = avg, col = ppnMixed)) +
  geom_line() +
  scale_color_manual(name = "Mixed-stock Fishery\nAllocation",
                     values = colPal) +
  labs(x = "Fixed Exploitation Rate", y = "Performance Metric Median") +
  theme_sleekX() + 
  facet_wrap(~labelledVar, scales = "free_y")


## Changes in realized exploitation rate 
ggplot(lpAgDat %>% filter(var == "medER"), 
       aes(x = mp, y = avg, col = ppnMixed)) +
  geom_line() +
  geom_ribbon(aes(ymin = lowQ, ymax = highQ, fill = ppnMixed), alpha = 0.2) +
  scale_fill_manual(values = colPal, guide = FALSE) +
  scale_color_manual(name = "Mixed-stock Fishery\nAllocation",
                     values = colPal) +
  labs(x = "Target Mixed-Stock Exploitation Rate", 
       y = "Realized Exploitation Rate") +
  theme_sleekX()

## Trade-off plots across different fixed exploitation rates
simParTrim <- simPar %>% 
  filter(nameOM %in% c("lowProd", "ref", "fullEnRoute_lowProd"))
scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species), 
                                                sep = "_"))

multOMAgDat <- buildDataAgg(dirNames, agVars =  agVarsToPlot, 
                         keyVarName = "ppnMixed") %>%
  #pool MPs by fixedER
  mutate(mp = as.numeric(case_when(
    str_detect(mp, "0.1fixedER") ~ "0.1",
    str_detect(mp, "0.3fixedER") ~ "0.3",
    str_detect(mp, "0.5fixedER") ~ "0.5",
    str_detect(mp, "0.7fixedER") ~ "0.7",
    str_detect(mp, "0.9fixedER") ~ "0.9"
  ))) %>% 
  filter(mp %in% c("0.3", "0.7"),
         om %in% c("ref", "lowProd", "fullEnRoute_lowProd")) %>% 
  mutate(om = fct_relevel(om, "ref", after = Inf)) %>% 
  mutate(mp = as.factor(mp),
         om = fct_recode(om, "Low Productivity" = "lowProd",
                         "En Route Mortality" = "fullEnRoute_lowProd",
                         "High Productivity" = "ref")) 

consVec <- rep(c("medSpawners", "ppnCULower"), each = 1)
catchVec <- rep(c("medCatch"), times = 2)
toAgPlotList <- lapply(seq_along(consVec), function(x) {
  p <- plotAgTradeoff(multOMAgDat, consVar = consVec[x], 
                      catchVar = catchVec[x],
                      facet = "om", shape = "mp", showUncertainty = TRUE,
                      xLab = catchVec[x], yLab = consVec[x], 
                      legendLab = "Proportion\nTAC Mix\nFishery",
                      axisSize = 11, dotSize = 4, legendSize = 13,
                      lineSize = 0.5, scaleAxis = "fixed")
  return(p)
})
sapply(toAgPlotList, function(x) print(x))


## CU-specific plots
cuVarsToPlot <- c("medSpawners", 
                  # "medRecRY", 
                  "medCatch", "ppnYrsUpper")

cuDat <- buildDataCU(dirNames = dirNames, 
                     cuVars = c(cuVarsToPlot, "medTotalER"),
                     keyVarName = "ppnMixed", 
                     selectedCUs = summCUs) %>% 
  mutate(mp = as.numeric(case_when(
    str_detect(mp, "0.1fixedER") ~ "0.1",
    str_detect(mp, "0.3fixedER") ~ "0.3",
    str_detect(mp, "0.5fixedER") ~ "0.5",
    str_detect(mp, "0.7fixedER") ~ "0.7",
    str_detect(mp, "0.9fixedER") ~ "0.9")),
    ppnMixed = as.factor(ppnMixed),
    cuName = fct_recode(cuName, "L. Stuart" = "L.St", "Stellako" = "Stll",
                        "Quesnel" = "Qsnl", "Chilko" = "Chlk", "Harrison" = 
                          "Hrrs")) %>% 
  filter(om == "lowProd",
         mp %in% c("0.1", "0.3", "0.5", "0.7", "0.9")) %>% 
  select(-plotOrder, -muName)

#Line plots
colPal <- viridis::viridis(length(unique(cuDat$ppnMixed)), begin = 0, 
                           end = 1, option = "plasma")
names(colPal) <- levels(cuDat$ppnMixed)

for(i in seq_along(cuVarsToPlot)) {
  dum <- cuDat %>% 
    filter(var == cuVarsToPlot[i])
  p <- ggplot(dum, aes(x = mp, y = avg, col = ppnMixed)) +
    geom_line() +
    scale_color_manual(name = "Mixed-stock Fishery\nAllocation",
                       values = colPal) +
    labs(x = "Fixed Exploitation Rate", y = cuVarsToPlot[i]) +
    theme_sleekX() + 
    facet_wrap(~cuName, scales = "free_y")
  print(p)
}

png(file = paste(here(),"/figs/suppMS/FigS1_cuRealizedER.png", sep = ""),
    height = 3.75, width = 6, units = "in", res = 300)
ggplot(cuDat %>% filter(var == "medTotalER"), 
       aes(x = mp, y = avg, col = ppnMixed)) +
  geom_ribbon(aes(ymin = lowQ, ymax = highQ, fill = ppnMixed), alpha = 0.2,
              colour = NA) +
  geom_line() +
  scale_fill_manual(values = colPal, guide = FALSE) +
  scale_color_manual(name = "Mixed-stock Fishery\nAllocation",
                     values = colPal) +
  labs(x = "Target Mixed-Stock Exploitation Rate", 
       y = "Realized Exploitation Rate") +
  theme_sleekX(legendSize = 0.8, axisSize = 9) +
  facet_wrap(~cuName, scales = "free_y")
dev.off()

png(file = paste(here(),"/figs/suppMS/FigS2_cuSpwnrs.png", sep = ""),
    height = 3.75, width = 6, units = "in", res = 300)
ggplot(cuDat %>% filter(var == "medSpawners"), 
       aes(x = mp, y = avg, col = ppnMixed)) +
  # geom_ribbon(aes(ymin = lowQ, ymax = highQ, fill = ppnMixed), alpha = 0.2,
  #             colour = NA) +
  geom_line() +
  scale_fill_manual(values = colPal, guide = FALSE) +
  scale_color_manual(name = "Mixed-stock Fishery\nAllocation",
                     values = colPal) +
  labs(x = "Target Mixed-Stock Exploitation Rate", 
       y = "Median Spawner Abundance (millions)") +
  theme_sleekX(legendSize = 0.8, axisSize = 9) +
  facet_wrap(~cuName, scales = "free_y")
dev.off()

png(file = paste(here(),"/figs/suppMS/FigS3_cuCatch.png", sep = ""),
    height = 3.75, width = 6, units = "in", res = 300)
ggplot(cuDat %>% filter(var == "medCatch"), 
       aes(x = mp, y = avg, col = ppnMixed)) +
  # geom_ribbon(aes(ymin = lowQ, ymax = highQ, fill = ppnMixed), alpha = 0.2,
  #             colour = NA) +
  geom_line() +
  scale_fill_manual(values = colPal, guide = FALSE) +
  scale_color_manual(name = "Mixed-stock Fishery\nAllocation",
                     values = colPal) +
  labs(x = "Target Mixed-Stock Exploitation Rate", 
       y = "Median Catch (millions)") +
  theme_sleekX(legendSize = 0.8, axisSize = 9) +
  facet_wrap(~cuName, scales = "free_y")
dev.off()

#Trade off plots
aggCatchList <- lapply(seq_along(summCUs), function (i) {
  refAgDat %>% 
    filter(var == "medCatch") %>% 
    mutate(cuName = summCUs[i],
           var = fct_recode(var, "medAggCatch" = "medCatch")) %>% 
    select(ppnMixed, mp, om, hcr, cuName, var, avg, lowQ, highQ)
}) 
cuDatFull <- do.call(rbind, aggCatchList) %>% 
  rbind(., cuRefDat) %>% 
  mutate(cuName = fct_recode(cuName, "L. Stuart" = "L.St", "Stellako" = "Stll",
                             "Quesnel" = "Qsnl", "Chilko" = "Chlk", "Harrison" = 
                               "Hrrs"))
cuDatLowERLow <- cuDatFull %>%
  filter(mp == "0.3")

# cuDatModER <- cuDatFull %>% 
#   filter(mp == "0.5")


consVec <- cuVarsToPlot#c("medSpawners", "medRecRY", "ppnYrsUpper")
yLabVec <- c("Median Spawner Abundance", "Median Return Abundance", "Proportion
             of Years Above Upper Benchmark")

toPlotList <- lapply(seq_along(consVec), function (h) {
  p <- plotCUTradeoff(cuDatLowERLow, consVar = consVec[h], catchVar = "medAggCatch",
                      facet = "cu", panel = "mp", showUncertainty = TRUE,
                      xLab = "Median Aggregate Catch", yLab = yLabVec[h],
                      legendLab = "Proportion Mixed\nStock TAC",
                      lineSize = 0.3)
  return(p)
})

sapply(toPlotList, function (h)
  sapply(h, function(i) print(i))
)


refAgDat <- buildDataAgg(dirNames, agVars =  agVarsToPlot, 
                        keyVarName = "ppnMixed") %>%
  #pool MPs by fixedER
  mutate(mp = as.numeric(case_when(
    str_detect(mp, "0.1fixedER") ~ "0.1",
    str_detect(mp, "0.3fixedER") ~ "0.3",
    str_detect(mp, "0.5fixedER") ~ "0.5",
    str_detect(mp, "0.7fixedER") ~ "0.7",
    str_detect(mp, "0.9fixedER") ~ "0.9"
  ))) %>% 
  filter(om == "ref",
         mp %in% c("0.1", "0.3", "0.5", "0.7", "0.9")) %>%
  mutate(labelledVar = var) %>%
  mutate(labelledVar = fct_recode(labelledVar, 
                                  "Aggregate\nReturn Abundance" = "medRecRY",
                                  "Aggregate\nSpawner Abundance" = 
                                    "medSpawners",
                                  "Ppn. MUs Above\nEsc. Goal" = "ppnMixedOpen",
                                  "Ppn. CUs Not\nRed (Throughout)" = 
                                    "ppnCULower",
                                  "Ppn. CUs Not\nRed (Sim. End)" = 
                                    "ppnCUStable",
                                  "Ppn. CUs Extant" = 
                                    "ppnCUExtant",
                                  "Aggregate Catch" = "medCatch",
                                  "Ppn. Yrs. Above\nCatch Target" = 
                                    "ppnYrsHighCatch",
                                  "Catch Stability" = "stabilityCatch"))


cuRefDat <- buildDataCU(dirNames = dirNames, 
                     cuVars = cuVarsToPlot,
                     keyVarName = "ppnMixed", 
                     selectedCUs = summCUs) %>% 
  mutate(mp = as.numeric(case_when(
    str_detect(mp, "0.1fixedER") ~ "0.1",
    str_detect(mp, "0.3fixedER") ~ "0.3",
    str_detect(mp, "0.5fixedER") ~ "0.5",
    str_detect(mp, "0.7fixedER") ~ "0.7",
    str_detect(mp, "0.9fixedER") ~ "0.9")),
    ppnMixed = as.factor(ppnMixed)) %>% 
  filter(om == "ref",
         mp %in% c("0.1", "0.3", "0.5", "0.7", "0.9")) %>% 
  select(-plotOrder, -muName)
