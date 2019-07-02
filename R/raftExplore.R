## Examine why Raft outputs change so dramatically at 100% mixed exploitation
# rate

listOfPackages <- c("here", "parallel", "doParallel", "foreach", "tidyverse",
                    "tictoc", "samSim")
newPackages <- listOfPackages[!(listOfPackages %in% 
                                  installed.packages()[ , "Package"])]
if (length(newPackages)) {
  install.packages(newPackages)
}
lapply(listOfPackages, require, character.only = TRUE)

simPar <- read.csv(here("data", "manProcScenarios",
                        "fraserMPInputs_varyAllocationAdjSRPars.csv"),
                   stringsAsFactors = F)
cuPar <- read.csv(here("data/fraserCUparsNEW.csv"), stringsAsFactors = F)

nTrials <- 150
## Make unique MP vector 
simPar$nameMP <- paste(simPar$propMixHigh, simPar$singleHCR, "_", simPar$canER, 
                       simPar$harvContRule, sep = "")

simParTrim <- simPar %>%
  filter(canER %in% c("0.3", "0.7"))
scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, "stableOut",
                                                # unique(simParTrim$species), 
                                                sep = "_"))

cuVars <- c("medSpawners", "medRecRY", "medCatch", "ppnYrsLower", 
            "medEarlyRecRY", "medEstAlpha", "medTotalER", "ppnYrsSingleOpen")
## Build data
cuData = NULL #construct CU dataframe
for (i in seq_along(dirNames)) {
  #alternatively ID OMs based on prespecified directory
  subDirs <- list.dirs(path = paste(here::here("outputs/simData"),
                                    dirNames[i], sep = "/"),
                       full.names = FALSE, recursive = FALSE)
  
  cuList <- genOutputList(dirNames[i], subDirs, agg = FALSE, 
                          selectedCUs = "Raft")
  nCU <- length(cuList[[1]]$stkName)
  keyVarValue <- sapply(cuList, function(x) unique(x$keyVar))
  cuName <- as.vector(sapply(cuList, function(x) x$stkName))
  singleScenDat = NULL
  for (k in seq_along(cuVars)) {
    keyNames <- sapply(cuList, function(x) x[["keyVar"]])
    dum <- sapply(cuList, function(x) x[[cuVars[k]]])
    # colnames(dum) <- keyNames
    dum2 <- dum %>% 
      reshape2::melt() %>% 
      dplyr::rename("trial" = "Var1", "mpFull" =  "Var2") %>% 
      mutate(var = as.factor(cuVars[k]),
        ppnMix = as.numeric(case_when(
          str_detect(mpFull, "0retro") ~ "0",
          str_detect(mpFull, "0.25retro") ~ "0.25",
          str_detect(mpFull, "0.5retro") ~ "0.5",
          str_detect(mpFull, "0.75retro") ~ "0.75",
          str_detect(mpFull, "1retro") ~ "1"
        )),
        mp = as.numeric(case_when(
          str_detect(mpFull, "0.3fixedER") ~ "0.3",
          str_detect(mpFull, "0.7fixedER") ~ "0.7")),
        om = unique(sapply(cuList, function(x) x$opMod))) %>% 
      select(-mpFull)
    singleScenDat <- rbind(singleScenDat, dum2, row.names = NULL)
  }
  cuData <- rbind(cuData, singleScenDat)
}
raftData <- cuData %>% 
  arrange(ppnMix, mp) %>% 
  mutate(ppnMix = as.factor(ppnMix),
         om = as.factor(om),
         mp = as.factor(mp))

ggplot(raftData %>% filter(mp == "0.3", om == "ref"), 
       aes(x = ppnMix, y = value)) +
  geom_boxplot() +
  facet_wrap(~var, scales = "free_y")

## Only explanation is overcompensation; closures at low abundance ironically
# sufficient harvest to allow stock to rebound