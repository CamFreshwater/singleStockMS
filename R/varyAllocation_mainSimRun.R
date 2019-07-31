#******************************************************************************
# varyAllocation_mainSimRun.R
# Date revised: May 10, 2019
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
# models

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
  filter(nameOM %in% c("ref", "enRoute"),
         benchmark == "stockRecruit")
scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species),
                                                sep = "_"))
# dirNames <- sapply(scenNames, function(x) paste("forePpn", x, 
#                                                 unique(simParTrim$species), 
#                                                 sep = "_"))

# Focal CUs and CU-specific PMs
summCUs <- cuPar %>% 
  filter(manUnit == "Summ") %>% 
  mutate(abbStkName = abbreviate(stkName, minlength = 4)) %>% 
  select(abbStkName) %>% 
  unlist()
  

#------------------------------------------------------------------------------

## Run simulation
# recoverySim(simParTrim[1, ], cuPar, catchDat = catchDat, srDat = srDat,
#             variableCU = FALSE, ricPars, larkPars = larkPars, tamFRP = tamFRP,
#             dirName = "test", nTrials = 5, makeSubDirs = FALSE, random = FALSE)

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
agVarsToPlot <-  c("medRecRY", "medSpawners", "ppnCULower", 
                   "medCatch", "stabilityCatch", "medER")

agDat <- buildDataAgg(dirNames, agVars =  agVarsToPlot, 
                      keyVarName = "ppnMixed") %>%
  mutate(mp = case_when(
    str_detect(mp, "genPA") ~ "Generic",
    str_detect(mp, "TAM") ~ "TAM"
  ), 
  labelledVar = var) %>%
  mutate(mp = as.factor(mp), 
         labelledVar = fct_recode(labelledVar, 
                                  "Aggregate\nReturn Abundance" = "medRecRY",
                                  "Aggregate\nSpawner Abundance" = 
                                    "medSpawners",
                                  "Ppn. CUs Not\nRed (Throughout)" = 
                                    "ppnCULower",
                                  "Aggregate Catch" = "medCatch",
                                  "Catch Stability" = "stabilityCatch",
                                  "Median Exploitation Rate" = "medER")) 

saveRDS(agDat, here::here("outputs", "generatedData", "mainSimAgDat.rds"))

refDat <- agDat %>% 
  filter(om == "ref", 
         mp == "Generic") 

# Fig. 1 - Dot plots
colPal <- viridis::viridis(length(levels(refDat$ppnMixed)), begin = 0, 
                           end = 1)
names(colPal) <- levels(refDat$ppnMixed)

agDotPlotList <- lapply(seq_along(agVarsToPlot), function(x) {
  dum <- refDat %>% 
    filter(var == agVarsToPlot[x])
  p <- ggplot(dum, aes(x = ppnMixed, y = avg, ymin = lowQ, ymax = highQ,
                          fill = ppnMixed)) +
    labs(x = "Proportion Mixed-Stock Fishery", y = unique(dum$labelledVar)) +
    geom_pointrange(shape = 21, fatten = 4, size = 1,
                    position = position_dodge(width = 0.65)) +
    theme_sleekX() +
    scale_fill_manual(guide = FALSE, values = colPal)
  return(p)
})
ggpubr::ggarrange(agDotPlotList[[1]], agDotPlotList[[2]], agDotPlotList[[3]],
                  agDotPlotList[[4]], agDotPlotList[[5]], agDotPlotList[[6]],
                  ncol = 3, nrow = 2, common.legend = TRUE, legend = "right")
#(remove x axis in top row)

# Fig. 2 - Agg trade-off 
consVec <- c("medRecRY", "ppnCULower")
catchVec <- rep(c("medCatch"), times = length(consVec))
conLabels <- c("Median Aggregate Return Abundance", 
               "Proportion of CUs Above Benchmark")


toAgList <- lapply(seq_along(consVec), function(x) {
  p <- plotAgTradeoff(refDat %>% select(-labelledVar), consVar = consVec[x], 
                      catchVar = catchVec[x],
                      facet = NULL, shape = NULL, showUncertainty = TRUE,
                      xLab = "Median Aggregate Catch", yLab = conLabels[x], 
                      legendLab = "Proportion\nTAC Mix\nFishery",
                      axisSize = 11, dotSize = 4, legendSize = 13,
                      lineSize = 0.5, scaleAxis = "fixed")
  return(p)
})
ggpubr::ggarrange(toAgList[[1]], toAgList[[2]],
          ncol = 2, nrow = 1, common.legend = TRUE, legend = "right",
          align = "h", widths = c(1, 1))

# Fig. S2 - Agg trade-off with multiple OMs 
enRouteDat <- agDat %>% 
  filter(mp == "Generic") %>% 
  mutate(om = case_when(
    str_detect(om, "enRoute") ~ "En-Route Mortality\nBefore Fishery",
    str_detect(om, "ref") ~ "En-Route Mortality\nAfter Fishery"
  ))
conLabels2 <- c("Median Aggregate\nReturn Abundance", 
               "Proportion of CUs\nAbove Benchmark")


toAgList <- lapply(seq_along(consVec), function(x) {
  p <- plotAgTradeoff(enRouteDat %>% select(-labelledVar), consVar = consVec[x], 
                      catchVar = catchVec[x],
                      facet = "om", shape = NULL, showUncertainty = TRUE,
                      xLab = "Median Aggregate Catch", yLab = conLabels2[x], 
                      legendLab = "Proportion\nTAC Mix\nFishery",
                      axisSize = 11, dotSize = 4, legendSize = 13,
                      lineSize = 0.5, scaleAxis = "fixed")
  return(p)
})
ggpubr::ggarrange(toAgList[[1]], toAgList[[2]],
                  ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")

# Fig. S3 - Agg trade-off with multiple MPs 
refDat2 <- agDat %>% 
  filter(om == "ref")

toAgList <- lapply(seq_along(consVec), function(x) {
  p <- plotAgTradeoff(refDat2 %>% select(-labelledVar), consVar = consVec[x], 
                      catchVar = catchVec[x],
                      facet = NULL, shape = "mp", showUncertainty = TRUE,
                      xLab = "Median Aggregate Catch", yLab = conLabels[x], 
                      legendLab = "Proportion\nTAC Mix\nFishery",
                      axisSize = 11, dotSize = 4, legendSize = 13,
                      lineSize = 0.5, scaleAxis = "fixed")
  return(p)
})
ggpubr::ggarrange(toAgList[[1]], toAgList[[2]],
                  ncol = 2, nrow = 1, common.legend = TRUE, legend = "right",
                  align = "h", widths = c(1, 1))


## CU-specific plots
cuVarsToPlot <- c("medRecRY", "medSpawners", "medCatch", "ppnYrsLower", 
                  "medTotalER")

cuDat <- buildDataCU(dirNames = dirNames, 
                     cuVars = c(cuVarsToPlot, "medTotalER"),
                     keyVarName = "ppnMixed", 
                     selectedCUs = summCUs) %>% 
  mutate(ppnMixed = as.factor(ppnMixed),
         mp = case_when(
           str_detect(mp, "genPA") ~ "Generic",
           str_detect(mp, "TAM") ~ "TAM"
         ),
         cuName = fct_recode(cuName, "L. Stuart" = "L.St", "Stellako" = "Stll",
                        "Quesnel" = "Qsnl", "Chilko" = "Chlk", "Harrison" = 
                          "Hrrs")) %>% 
  select(-plotOrder, -muName)


#Trade off plots
aggCatchList <- lapply(seq_along(summCUs), function (i) {
  refDat %>% 
    filter(var == "medCatch") %>% 
    mutate(cuName = summCUs[i],
           var = fct_recode(var, "medAggCatch" = "medCatch")) %>% 
    select(ppnMixed, mp, om, hcr, cuName, var, avg, lowQ, highQ)
}) 
cuDatPlot <- do.call(rbind, aggCatchList) %>% 
  rbind(., cuDat) %>% 
  mutate(cuName = fct_recode(cuName, "L. Stuart" = "L.St", "Stellako" = "Stll",
                             "Quesnel" = "Qsnl", "Chilko" = "Chlk", "Harrison" = 
                               "Hrrs")) %>% 
  filter(om == "ref",
         mp == "Generic")

## Fig 3
plotCUTradeoff(cuDatPlot, consVar = "medRecRY", catchVar = "medCatch",
               facet = "cu", panel = NULL, showUncertainty = TRUE,
               xLab = "Median Catch", yLab = "Median Return Abundance",
               legendLab = "Proportion\nMixed\nStock TAC", lineSize = 0.3, 
               main = FALSE)
## Fig S3
plotCUTradeoff(cuDatPlot, consVar = "ppnYrsLower", catchVar = "medCatch",
               facet = "cu", panel = NULL, showUncertainty = TRUE,
               xLab = "Median Catch", 
               yLab = "Mean Proportion Years\nAbove Lower BM",
               legendLab = "Proportion\nMixed\nStock TAC", lineSize = 0.3, 
               main = FALSE)

colPal <- viridis::viridis(length(levels(refDat$ppnMixed)), begin = 0, 
                           end = 1)
names(colPal) <- levels(refDat$ppnMixed)
cuER <- cuDatPlot %>% 
  filter(var == "medTotalER")
ggplot(cuER, 
       aes(x = ppnMixed, y = avg, ymin = lowQ, ymax = highQ, fill = ppnMixed)) +
  labs(x = "Ppn. Mixed-Stock\nFishery", y = "Exp. Rate") +
  geom_pointrange(shape = 21, fatten = 4, size = 1,
                  position = position_dodge(width = 0.65)) +
  scale_fill_manual(guide = FALSE, values = colPal) +
  facet_wrap(~cuName)
#------------------------------------------------------------------------------

## Run simulation with full gradient of productivity and synchrony operating
# models to produce heat maps

simPar <- read.csv(here("data", "manProcScenarios",
                        "fraserMPInputs_varyAllocationProdRange.csv"),
                   stringsAsFactors = F)

## Expand dataset to include full range of productivity regimes 
prodScalars <- seq(0.5, 1.1, by = 0.05)
synchScalars <- seq(0.05, 0.65, by = 0.05)
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

heatMapDat <- buildDataAgg(dirNames, agVars =  agVarsToPlot, 
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

# saveRDS(heatMapDat, here::here("outputs", "generatedData", "heatMapDat.rds"))


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
    labs(x = xLab, y = "Proportion Mixed Stock") +
    theme_sleekX()
}
synchDat <- heatMapDat %>% 
  filter(scenario == "Synch")
prodDat <- heatMapDat %>% 
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

#------------------------------------------------------------------------------

## Misc. supplementary figures

# Generate data
uMsy <- cuPar$uMSY[1]
lowFRP <- cuPar$lowFRP[1]
highFRP <- cuPar$highFRP[1]


ret = seq(0, 2, length.out = 999)
esc = rep(NA, length.out = 999)
hr = rep(NA, length.out = 999)
for (i in 1:length(ret)) {
  if (ret[i] < lowFRP) {
    hr[i] <- 0.1
  }
  if (ret[i] > lowFRP & ret[i] < highFRP) {
    hr[i] <- max(0.1, 
                 uMsy * ((ret[i] - lowFRP) / (highFRP - lowFRP))
                 )
  }
  if (ret[i] > highFRP) {
    hr[i] <- uMsy
  }
  esc[i] <- ret[i] - (ret[i] * hr[i])
}

plot(hr ~ ret, xlim=c(0, 2), ylim=c(0, 1), ylab = "", xlab = "", axes = F, 
     type ="l", lwd = 2)
axis(1, las = 1)
axis(2, las = 1)
mtext("Total Exploitation Rate", 2, line = 2.5)
mtext("Return Abundance (Millions)", 1, line = 2)
abline(v = lowFRP, col = "red", lty = 2, lwd = 2)
abline(v = highFRP, col = "darkgreen", lty = 2, lwd = 2)

