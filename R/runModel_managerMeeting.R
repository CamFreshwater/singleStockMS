#******************************************************************************
# managerMeeting.R
# Date revised: Mar. 4, 2019
# Explainer: Contains simulation runs and figures generated for April meeting
# with managers
#******************************************************************************


listOfPackages <- c("here", "parallel", "doParallel", "foreach", "viridis",
                    "tidyverse", "tictoc", "samSim")
newPackages <- listOfPackages[!(listOfPackages %in% 
                                  installed.packages()[ , "Package"])]
if (length(newPackages)) {
  install.packages(newPackages)
}
lapply(listOfPackages, require, character.only = TRUE)

simPar <- read.csv(here("data", "manProcScenarios",
                        "fraserMPInputs_managerMeeting.csv"),
                   stringsAsFactors = F)
cuPar <- read.csv(here("data/fraserCUparsNEW.csv"), stringsAsFactors = F)
srDat <- read.csv(here("data/fraserRecDatTrim.csv"), stringsAsFactors = F)
catchDat <- read.csv(here("data/fraserCatchDatTrim.csv"), stringsAsFactors = F)
ricPars <- read.csv(here("data/constrainedRickerMCMCPars.csv"), stringsAsFactors = F)
larkPars <- read.csv(here("data/constrainedLarkinMCMCPars.csv"), 
                     stringsAsFactors = F)
tamFRP <- read.csv(here("data/tamRefPts.csv"), stringsAsFactors = F)

## Define simulations to be run
nTrials <- 150

for(k in 1:nrow(simPar)) {
  if(is.na(simPar$nameMP[k])) {
    simPar$nameMP[k] <- paste(simPar$propMixHigh[k], simPar$singleHCR[k], "_", 
                              simPar$canER[k], simPar$harvContRule[k], sep = "")
  } 
}

simParTrim <- simPar  
# %>%
#     filter(keyVar == "ppnMix",
#            scenario == "fixed0.7Retro_enRoute")
scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species), 
                                                sep = "_"))
agVarsToPlot <-  c("medRecRY", "medSpawners", "ppnCUUpper", "ppnCUExtinct",
                   "medCatch")

#------------------------------------------------------------------------------

## Run model
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

### Reference productivity plots
refAgDat <- buildDataAgg(dirNames, agVars =  agVarsToPlot, 
                         keyVarName = "expRate") %>% 
  filter(om == "ref",
         !str_detect(mp, "retro"))

# Double y-axis plot
png(file = paste(here(),"/figs/april2019Meeting/singAxisPlot_noCatch.png", 
                 sep = ""), height = 4, width = 5, units = "in", res = 300)
plotContTradeOffs(refAgDat %>% filter(!var %in% ("medCatch")), 
                  keyVar = "expRate", double = FALSE)
dev.off()

png(file = paste(here(),"/figs/april2019Meeting/singAxisPlot.png", sep = ""),
    height = 4, width = 5, units = "in", res = 300)
plotContTradeOffs(refAgDat, keyVar = "expRate", double = FALSE)
dev.off()

png(file = paste(here(),"/figs/april2019Meeting/dblAxisPlot.png", sep = ""),
    height = 4, width = 5, units = "in", res = 300)
plotContTradeOffs(refAgDat, keyVar = "expRate", double = TRUE)
dev.off()

png(file = paste(here(),"/figs/april2019Meeting/aggTO_noTAM_spwn_ref.png", 
                 sep = ""), height = 4, width = 5, units = "in", res = 300)
plotAgTradeoff(refAgDat %>% filter(!hcr == "TAM"), consVar = "medSpawners", 
               catchVar = "medCatch",
               facet = NULL, shape = NULL, showUncertainty = TRUE, 
               mainLab = "", legendLab = "Exploitation Rate",
               xLab = "Median Catch", yLab = "Median Spawners", 
               axisSize = 13, legendSize = 12)
dev.off()

png(file = paste(here(),"/figs/april2019Meeting/aggTO_spwn_ref.png", 
                 sep = ""), height = 4, width = 5, units = "in", res = 300)
plotAgTradeoff(refAgDat, consVar = "medSpawners", 
               catchVar = "medCatch",
               facet = NULL, shape = NULL, showUncertainty = TRUE, 
               mainLab = "", legendLab = "Exploitation Rate",
               xLab = "Median Catch", yLab = "Median Spawners", 
               axisSize = 13, legendSize = 12)
dev.off()

png(file = paste(here(),"/figs/april2019Meeting/aggTO_ppnCU_ref.png", sep = ""),
    height = 4, width = 5, units = "in", res = 300)
plotAgTradeoff(refAgDat, consVar = "ppnCUUpper", catchVar = "medCatch",
               facet = NULL, shape = NULL, showUncertainty = TRUE, 
               mainLab = "", legendLab = "Exploitation Rate",
               xLab = "Median Catch", yLab = "Ppn CUs Above Upper BM", 
               axisSize = 13, legendSize = 12)
dev.off()


# CU-specific trade-off plots
summCUs <- cuPar %>% 
  filter(manUnit == "Summ") %>% 
  mutate(abbStkName = abbreviate(stkName, minlength = 4)) %>%
  select(abbStkName) %>%
  unlist()

cuVarsToPlot <- c("medSpawners", "ppnYrsLower", "medCatch")
cuDat <- buildDataCU(dirNames = dirNames, 
                     cuVars = cuVarsToPlot,
                     keyVarName = "expRate", 
                     selectedCUs = summCUs) %>% 
  filter(om == "ref",
         !str_detect(mp, "retro")) %>%
  select(-plotOrder, -muName)

# Merge aggregate catch data w/ cuDat
aggCatchList <- lapply(seq_along(summCUs), function (i) {
  refAgDat %>% 
    filter(var == "medCatch") %>% 
    mutate(cuName = summCUs[i],
           var = fct_recode(var, "medAggCatch" = "medCatch")) %>% 
    select(expRate, mp, om, hcr, cuName, var, avg, lowQ, highQ)
}) 
cuDatFull <- do.call(rbind, aggCatchList) %>% 
  rbind(., cuDat) %>% 
  mutate(cuName = fct_recode(cuName, "L. Stuart" = "L.St", "Stellako" = "Stll",
                             "Quesnel" = "Qsnl", "Chilko" = "Chlk", "Harrison" = 
                               "Hrrs"))

png(file = paste(here(),"/figs/april2019Meeting/cuTO_ppnBM_ref.png", sep = ""),
    height = 4, width = 7, units = "in", res = 300)
plotCUTradeoff(cuDatFull, consVar = "ppnYrsLower", catchVar = "medAggCatch", 
               facet = "cu", panel = "om", showUncertainty = TRUE,
               xLab = "Median Aggregate Catch", 
               yLab = "Ppn of Years Above Benchmark", 
               legendLab = "Exploitation Rate", main = FALSE, legendSize = 10,
               lineSize = 0.3, dotSize = 3, axisSize = 11, freeY = FALSE)
dev.off()

#------------------------------------------------------------------------------

### Changes in productivity plots
simParTrim <- simPar %>%
  filter(scenario %in% c("tam"))
scenNames <- unique(simParTrim$scenario)
prodDir <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species), 
                                                sep = "_"))
subProdDir <- unique(simParTrim$nameOM)[1:3]
stkNames <- genOutputList(prodDir, subDirName = subProdDir[1], 
                          agg = FALSE)[[1]][["stkName"]]
nCUs <- length(stkNames)
arrayNames <- sapply(subProdDir, function(x) { 
  list.files(paste(here("outputs/simData"), prodDir, x, sep="/"), 
             pattern = "\\Arrays.RData$")
})

# Function to pull and clean array data
outList <-  lapply(seq_along(arrayNames), function(h) {
  datList <- readRDS(paste(here("outputs/simData"), prodDir, subProdDir[h], 
                             arrayNames[h], sep = "/"))
  alpha <- datList[["alpha"]] %>% 
    reshape2::melt() %>% 
    dplyr::rename("yr" = "Var1", "cu" =  "Var2", "trial" = "Var3") %>% 
    mutate(om = as.factor(datList[["nameOM"]]), 
           var = "alpha")
  spwn <- datList[["S"]] %>% 
    reshape2::melt() %>% 
    dplyr::rename("yr" = "Var1", "cu" =  "Var2", "trial" = "Var3") %>% 
    mutate(om = as.factor(datList[["nameOM"]]), 
           var = "spwn")
  t <- rbind(alpha, spwn) %>% 
    mutate(stkName = as.factor(plyr::mapvalues(cu, from = unique(cu), 
                                               to = stkNames)))
})
outDat <- do.call(rbind, outList) %>% 
  mutate(om = fct_recode(om, "Median" = "ref", "Low" = "lowProd", 
                         "Divergent" = "divProd"))
  
# sample individual trials
# subCU <- sample.int(nCUs, size = 4)
trialID <- sample.int(unique(outDat$trial), size = 1)
standDat <- outDat %>% 
  filter(stkName %in% c("L.Sh", "Chlk", "Hrrs"),
         trial == trialID,
         yr > 50) 

colPal <- viridis(3, begin = 0, end = 1)
names(colPal) <- unique(standDat$stkName)

png(file = paste(here(),"/figs/april2019Meeting/prod_lines.png", sep = ""),
    height = 2, width = 6, units = "in", res = 300)
ggplot(standDat %>% filter(var == "alpha"), 
       aes(x = yr, y = value, col = stkName)) +
  geom_line() +
  scale_color_manual(values = colPal, guide = FALSE) +
  geom_vline(xintercept = 60, linetype = "dashed", size = 1) +
  labs(x = "Year", y = "Productivity") +
  theme_sleekX() + 
  facet_wrap(~om)
dev.off()

png(file = paste(here(),"/figs/april2019Meeting/spwn_lines.png", sep = ""),
    height = 2, width = 6, units = "in", res = 300)
ggplot(standDat %>% filter(var == "spwn"), 
       aes(x = yr, y = value, col = stkName)) +
  geom_line() +
  scale_color_manual(values = colPal, guide = FALSE) +
  geom_vline(xintercept = 60, linetype = "dashed", size = 1) +
  labs(x = "Year", y = "Spawner Abundance") +
  theme_sleekX() + 
  facet_wrap(~om)
dev.off()


### Synchrony impact plots
## Load in data from synchrony analyses to visualize impacts of greater 
#variance/covariance
synchDirNames <- c("lowSig_sockeye", "highSig_sockeye")
stkNames <- genOutputList(synchDirNames[1], 
                          agg = FALSE)[["medSynch_TAM"]][["stkName"]]
nCUs <- length(stkNames)
#matrix of array names to be passed
arrayNames <- sapply(synchDirNames[1], function(x) { 
  list.files(paste(here("outputs/simData"), x, sep="/"), 
             pattern = "\\Arrays.RData$")
})
synchNames <- c("high", "low", "med")
sigNames <- c("low", "high")

# Function to pull and clean array data
outList2 <- lapply(seq_along(synchDirNames), function(i) {
  outList1 <-  lapply(seq_along(arrayNames), function(h) {
    datList <- readRDS(paste(here("outputs/simData"), synchDirNames[i], 
                             arrayNames[h], 
                             sep = "/"))
    datList[["recDev"]] %>% 
      reshape2::melt() %>% 
      dplyr::rename("yr" = "Var1", "cu" =  "Var2", "trial" = "Var3", 
                    "recDev" = "value") %>% 
      mutate(sigma = sigNames[i], synch = synchNames[h])
  })
  do.call(rbind, outList1)
})
outDat <- do.call(rbind, outList2) 
outDat <- outDat %>% 
  mutate(scen = paste(sigma, "Sig", synch, "Synch", sep = ""),
         stkName = as.factor(plyr::mapvalues(outDat$cu, 
                                             from = unique(outDat$cu),
                                             to = stkNames)))
# standardize abundance within a CU and trim
subCU <- sample.int(nCUs, size = 3)
trialID <- sample.int(unique(outDat$trial), size = 1)
standDat <- outDat %>% 
  filter(yr > 60,
         cu %in% subCU,
         trial == trialID,
         !scen %in% c("lowSigmedSynch", "highSigmedSynch")) %>% 
  mutate(stkName = factor(stkName), 
         scen = recode(factor(scen), "lowSiglowSynch" = "Low Var. - Low Synch.",
                       "highSiglowSynch" = "High Var. - Low Synch.",
                       "highSighighSynch" = "High Var. - High Synch.",
                       "lowSighighSynch" = "Low Var. - High Synch.")) %>% 
  mutate(scen = factor(scen, levels(scen)[c(4, 2, 3, 1)]))

colPal <- viridis(length(subCU), begin = 0, end = 1)
names(colPal) <- unique(standDat$stkName)

png(file = paste(here(),"/figs/april2019Meeting/recDev_lines.png", sep = ""),
    height = 4, width = 6, units = "in", res = 300)
ggplot(standDat, aes(x = yr, y = recDev, col = stkName)) +
  geom_line() +
  scale_color_manual(values = colPal, name = "Conservation\nUnit") +
  labs(x = "Year", y = "Recruitment Deviations") +
  theme_sleekX() + 
  facet_wrap(~scen)
dev.off()

### Aggregate tradeoff plots
agVarsToPlot <-  c("medRecRY", "medSpawners", "ppnCUUpper", "ppnCUExtinct",
                   "medCatch")
prodAgDat <- buildDataAgg(dirNames, agVars =  agVarsToPlot, 
                         keyVarName = "expRate") %>% 
  filter(!mp %in% c("fixed0.2", "fixed0.4", "fixed0.6", "fixed0.8"),
         !str_detect(mp, "retro"),
         !om %in% c("highOU", "divProd")) %>% 
  mutate(om = fct_recode(om, "Reference Productivity" = "ref",
                         "Low Productivity" = "lowProd"))

png(file = paste(here(),"/figs/april2019Meeting/aggTO_spwn_low.png", sep = ""),
    height = 4, width = 7, units = "in", res = 300)
plotAgTradeoff(prodAgDat, consVar = "medSpawners", catchVar = "medCatch",
               facet = "om", shape = NULL, showUncertainty = TRUE, 
               mainLab = "", legendLab = "Exploitation Rate",
               xLab = "Median Catch", yLab = "Median Spawners", 
               axisSize = 13, legendSize = 11, freeY = FALSE)
dev.off()

# png(file = paste(here(),"/figs/april2019Meeting/aggTO_ppnCU_low.png", sep = ""),
#     height = 4, width = 7, units = "in", res = 300)
# plotAgTradeoff(prodAgDat, consVar = "ppnCUUpper", catchVar = "medCatch",
#                facet = "om", shape = NULL, showUncertainty = TRUE, 
#                mainLab = "", legendLab = "Exploitation Rate",
#                xLab = "Median Catch", yLab = "Ppn CUs Above Upper BM", 
#                axisSize = 13, legendSize = 11, freeY = FALSE)
# dev.off()

#------------------------------------------------------------------------------

### Changes in single-stock harvest rates
# simParSing
# simParTrim <- simPar %>%
#    filter(keyVar == "ppnMix")
# scenNames <- unique(simParTrim$scenario)
# dirNames <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species), 
#                                                 sep = "_"))
singleDirNames <- c("fixed0.3Retro_lowProd_sockeye", 
                    "fixed0.7Retro_lowProd_sockeye",
                    "fixed0.3Retro_enRouteLow_sockeye",
                    "fixed0.7Retro_enRouteLow_sockeye")
allocationDat <- buildDataAgg(singleDirNames, agVars =  agVarsToPlot, 
                          keyVarName = "ppnMix") %>% 
  filter(str_detect(mp, "retro")) %>% 
  mutate(#ppnMix = as.numeric(ppnMix),
         om = fct_recode(om, "Reference" = "lowProd",
                         "Mortality Before Fishery" = "fullEnRoute_lowProd"),
         mp = as.factor(case_when(
           str_detect(mp, "0.3fixed") ~ "Low Exp. Rate (0.3)",
           str_detect(mp, "0.7fixed") ~ "High Exp. Rate (0.7)"
         ))) %>% 
  filter(mp %in% c("Low Exp. Rate (0.3)", "High Exp. Rate (0.7)")) %>% 
  mutate(mp = fct_relevel(mp, "High Exp. Rate (0.7)", after = Inf))
  
noEnRoute <- allocationDat %>% 
  filter(om == "Reference") 
png(file = paste(here(),"/figs/april2019Meeting/aggTO_singleStock.png", 
                 sep = ""),
    height = 3, width = 5, units = "in", res = 300)
plotAgTradeoff(noEnRoute, consVar = "medSpawners", 
                    catchVar = "medCatch",
                    facet = "mp", showUncertainty = TRUE,
                    xLab = "Median Catch", yLab = "Median Spawners", 
                    legendLab = "Proportion\nTAC Mix\nFishery",
                    axisSize = 11, dotSize = 4, legendSize = 10,
                    lineSize = 0.5, scaleAxis = "fixed")
dev.off()

lowERDat <- allocationDat %>% 
  filter(mp == "Low Exp. Rate (0.3)")
png(file = paste(here(),"/figs/april2019Meeting/aggTO_singleStockER.png", 
                 sep = ""),
    height = 4, width = 6, units = "in", res = 300)
plotAgTradeoff(lowERDat, consVar = "medSpawners", 
               catchVar = "medCatch",
               facet = "om", shape = NULL, showUncertainty = TRUE,
               xLab = "Median Catch", yLab = "Median Spawners", 
               legendLab = "Proportion\nTAC Mix\nFishery",
               axisSize = 11, dotSize = 4, legendSize = 10,
               lineSize = 0.5, scaleAxis = "fixed")
dev.off()


# CU-specific trade-off plots
summCUs <- cuPar %>% 
  filter(manUnit == "Summ") %>% 
  mutate(abbStkName = abbreviate(stkName, minlength = 4)) %>% 
  select(abbStkName) %>% 
  unlist()

cuVarsToPlot <- c("medSpawners", "ppnYrsLower", "medCatch", 
                  "medTotalER", "ppnYrsSingleOpen", "ppnYrsUpper")
cuDat <- buildDataCU(dirNames = singleDirNames, 
                     cuVars = cuVarsToPlot,
                     keyVarName = "ppnMix", 
                     selectedCUs = summCUs) %>% 
  filter(str_detect(mp, "retro")) %>% 
  mutate(om = fct_recode(om, "Reference" = "lowProd",
                    "Mortality Before Fishery" = "fullEnRoute_lowProd"),
    mp = as.factor(case_when(
      str_detect(mp, "0.3fixed") ~ "Low Exp. Rate (0.3)",
      str_detect(mp, "0.7fixed") ~ "High Exp. Rate (0.7)"
    ))) %>% 
  filter(mp %in% c("Low Exp. Rate (0.3)", "High Exp. Rate (0.7)"), 
         om == "Reference") %>% 
  mutate(mp = fct_relevel(mp, "High Exp. Rate (0.7)", after = Inf)) %>% 
  select(-plotOrder, -muName)

aggCatchList <- lapply(seq_along(summCUs), function (i) {
  noEnRoute %>% 
    filter(var == "medCatch") %>% 
    mutate(cuName = summCUs[i],
           var = fct_recode(var, "medAggCatch" = "medCatch")) %>% 
    select(ppnMix, mp, om, hcr, cuName, var, avg, lowQ, highQ)
}) 
cuDatFull <- do.call(rbind, aggCatchList) %>% 
  rbind(., cuDat) %>% 
  mutate(cuName = fct_recode(cuName, "L. Stuart" = "L.St", "Stellako" = "Stll",
                             "Quesnel" = "Qsnl", "Chilko" = "Chlk", "Harrison" = 
                               "Hrrs"))

cuDLow <- cuDatFull %>% 
  filter(mp == "Low Exp. Rate (0.3)")
cuDHigh <- cuDatFull %>% 
  filter(mp == "High Exp. Rate (0.7)")

png(file = paste(here(),"/figs/april2019Meeting/cuTO_singleStock_highER.png", 
                 sep = ""), height = 4, width = 7, units = "in", res = 300)
plotCUTradeoff(cuDHigh, consVar = "ppnYrsLower", catchVar = "medAggCatch", 
               facet = "cu", panel = "om", showUncertainty = TRUE,
               xLab = "Median Aggregate Catch", 
               # yLab = "Proportion of Years Above Benchmark", 
               legendLab = "Proportion TAC\nin Mixed Stock", main = FALSE, 
               legendSize = 10, lineSize = 0.3, dotSize = 3, axisSize = 9, 
               scaleAxis =  "fixed")
dev.off()

# 
# plotCUTradeoff(cuDHigh, consVar = "medSpawners", catchVar = "medTotalER",
#                facet = "cu", panel = "om", showUncertainty = TRUE,
#                xLab = "Median Catch", main = FALSE, legendSize = 10,
#                lineSize = 0.3, dotSize = 3, axisSize = 11)
# plotCUTradeoff(cuDHigh, consVar = "medSpawners", catchVar = "ppnYrsSingleOpen",
#                facet = "cu", panel = "om", showUncertainty = TRUE,
#                xLab = "Median Catch", main = FALSE, legendSize = 10,
#                lineSize = 0.3, dotSize = 3, axisSize = 11)


#------------------------------------------------------------------------------

### Impacts of outcome uncertainty on Nass chum
nassDirNames <- c("variableOU", "variableOU_lowProd")

agVarsToPlot <-  c("medRecRY", "medSpawners", "ppnCULower", "ppnCUExtant",
                   "medCatch", "ppnYrsEscGoal")
agDat <- buildDataAgg(nassDirNames, agVars =  agVarsToPlot,
                      keyVarName = "ou") %>% 
  mutate(
    om = as.factor(case_when(
      str_detect(om, "medProd_") ~ "Ref. Prod",
      str_detect(om, "lowProd_") ~ "Low Prod")),
    ou = as.factor(as.numeric(as.character(ou)))
  ) %>% 
  mutate_when(
    var %in% c("medRecRY", "medSpawners", "medCatch"),
    list(avg = avg / 1000,
         lowQ = lowQ / 1000,
         highQ = highQ / 1000)
  ) %>% 
  filter(!mp == "0.3")

## Look at tradeoffs between different operating models
png(file = paste(here(),"/figs/april2019Meeting/aggTO_nassOUSpwn_refOMOneMP.png", 
                 sep = ""),
    height = 4.5, width = 6, units = "in", res = 300)
dum <- agDat %>% 
  filter(om == "Ref. Prod",
         mp == "0.1")
plotAgTradeoff(dum, consVar = "medSpawners", catchVar = "medCatch",
               facet = "om", shape = "mp", showUncertainty = FALSE, 
               mainLab = "", legendLab = "Outcome\nUncertainty",
               xLab = "Median Catch (thousands)", 
               yLab = "Median Spawners (thousands)", 
               axisSize = 13, legendSize = 11, lineSize = 0.5,
               scaleAxis = "fixed")
dev.off()

png(file = paste(here(),"/figs/april2019Meeting/aggTO_nassOUSpwn_OneMP.png", 
                 sep = ""),
    height = 4.5, width = 6, units = "in", res = 300)
dum <- agDat %>% 
  filter(mp == "0.1")
plotAgTradeoff(dum, consVar = "medSpawners", catchVar = "medCatch",
               facet = "om", shape = "mp", showUncertainty = FALSE, 
               mainLab = "", legendLab = "Outcome\nUncertainty",
               xLab = "Median Catch (thousands)", 
               yLab = "Median Spawners (thousands)", 
               axisSize = 13, legendSize = 11, lineSize = 0.5,
               scaleAxis = "fixed")
dev.off()

png(file = paste(here(),"/figs/april2019Meeting/aggTO_nassOUSpwn_points.png", 
                 sep = ""),
    height = 4.5, width = 8, units = "in", res = 300)
plotAgTradeoff(agDat, consVar = "medSpawners", catchVar = "medCatch",
               facet = "om", shape = "mp", showUncertainty = FALSE, 
               mainLab = "", legendLab = "Outcome\nUncertainty",
               xLab = "Median Catch (thousands)", 
               yLab = "Median Spawners (thousands)", 
               axisSize = 13, legendSize = 11, lineSize = 0.5,
               scaleAxis = "fixed")
dev.off()

png(file = paste(here(),"/figs/april2019Meeting/aggTO_nassOUSpwn.png", 
                 sep = ""),
    height = 4.5, width = 8, units = "in", res = 300)
plotAgTradeoff(agDat, consVar = "medSpawners", catchVar = "medCatch",
               facet = "om", shape = "mp", showUncertainty = TRUE, 
               mainLab = "", legendLab = "Outcome\nUncertainty",
               xLab = "Median Catch (thousands)", 
               yLab = "Median Spawners (thousands)", 
               axisSize = 13, legendSize = 11, lineSize = 0.5,
               scaleAxis = "fixed")
dev.off()

png(file = paste(here(),"/figs/april2019Meeting/aggTO_nassOUEscGoal.png", 
                 sep = ""),
    height = 4.5, width = 8, units = "in", res = 300)
plotAgTradeoff(agDat, consVar = "ppnYrsEscGoal", catchVar = "medCatch",
               facet = "om", shape = "mp", showUncertainty = TRUE, 
               mainLab = "", legendLab = "Outcome\nUncertainty",
               xLab = "Median Catch", 
               yLab = "Proportion of Years Above\n. Escapement Goal", 
               axisSize = 11, dotSize = 4, legendSize = 10,
               lineSize = 0.5, scaleAxis = "fixed")
dev.off()

