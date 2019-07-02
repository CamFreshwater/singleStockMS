#*******************************************************************************
# frStockRecModels_SummersPooled.r
# Date last revised: 29 April, 2019
# Modified version of frStockRecModels that fits a Ricker model to all the 
# summer data aggregated (trimmed by shortest TS). Goal is to use management 
# reference points estimated from those models to parameterize a generic mixed-
# stock HCR. Data are then pased to cuPars which is trimmed to be used in 
# single stock forward simulations. Finally SR parameters for each CU are
# replace by recursive Bayes estimates for use in forward simulations.
#
# NOTE: Unclear whether fisheries management reference points should also 
# account for non-stationarity. If so, aggregate SR analysis should be updated.
#*******************************************************************************

require(here); require(R2jags); require(ggplot2); require(tidyverse)
require(samSim)

source(here("R/func/stockRecFunctions.R"))


### Read in data 
cuPar <- read.csv(here("data/fraserCUParsNEW.csv"), stringsAsFactors=F)
summCUs <- cuPar %>% 
  filter(manUnit == "Summ") %>% 
  select(stkName)


# Clean data to match format for stockRecFunctions.R
# Use total spawners
summSRDat <- read.csv(here("data/fraserRecDat.csv"), stringsAsFactors=F) %>% 
  select(stk, stkName, yr, totalSpwn, rec) %>% 
  rename(CU = stk, Year = yr, Escape = totalSpwn, Recruit = rec) %>% 
  filter(stkName %in% summCUs$stkName) %>% #collapse to summers
  group_by(Year) %>% 
  summarize(Escape = sum(Escape),
            Recruit = sum(Recruit)) %>% 
  mutate(CU = "xx",
         stkName = "Pooled Summ.", 
         lnRS = log(Recruit / Escape))


ricFName = "fraserRicker_pooledSummers"
# RBayesSingle(dat = summSRDat, Niter=250000, Nthin=50, burnin=50000, 
#              fname = ricFName, model = "ricker", betaPrior = "logN", cap = 1, 
#              tauPrior = 0.001)


modelSummRick <- readRDS(here("outputs/srAnalysis/data", 
                              paste(ricFName, "SimOutput.RDS", 
                                    sep = "")))
ggplot(summSRDat, aes(x = Escape, y = Recruit)) +
  geom_point() +
  stat_function(fun = function(x) x * exp(1.8 - 0.4 * x))


rickerOut <- read.csv(here("outputs/srAnalysis/data", 
                              paste(ricFName, "MCMCPars.csv", sep = ""))) %>% 
  select(alpha, beta0, sigma) %>% 
  gather(key = parameter, value = estimate) %>% 
  group_by(parameter) %>% 
  summarize(mean =  mean(estimate))

summAlpha <- rickerOut[1, 2]
summBeta <- rickerOut[2, 2]

calcUmsy <- function(alpha){ #from Scheuerell 2016
  # as.numeric(alpha * (0.5 - ((0.65 * alpha^1.27) / (8.7 + alpha^1.27))))
  as.numeric(1 - gsl::lambert_W0(exp(1 - alpha)))
}
calcSmsy <- function(alpha, beta) {
  as.numeric((1 - gsl::lambert_W0(exp(1 - alpha))) / beta)
}

## Calculate reference points and add to pooled data
summOnlyPars <- cuPar %>% 
  filter(manUnit == "Summ") %>% 
  mutate(uMSY = calcUmsy(summAlpha),
         lowFRP = 0.4 * calcSmsy(summAlpha, summBeta),
         highFRP = 0.8 * calcSmsy(summAlpha, summBeta))

## Replace original CU pars (estimated with stationary model) with non-
# stationary estimates (note that some values seem off, but summers are ok)
rickFrRecursive <- read.csv(here("data", "old", "recursiveRickerMCMCPars.csv"),
                            stringsAsFactors = F)
larkFrRecursive <- read.csv(here("data", "old", "recursiveLarkinMCMCPars.csv"),
                            stringsAsFactors = F)

#calculate medians
rickFrPar <- rickFrRecursive %>% 
  group_by(stk.name) %>% 
  dplyr::rename(stkName = stk.name) %>% 
  summarize(alpha = median(last.alpha.v), beta0 = median(beta.v), 
            sigma = median(tlsig.v))
larkFrPar <- larkFrRecursive %>% 
  group_by(stk.name) %>% 
  dplyr::rename(stkName = stk.name) %>% 
  summarize(larkAlpha = median(last.alpha.v), larkBeta0 = median(beta0.v),
            larkBeta1 = median(beta1.v), larkBeta2 = median(beta2.v),
            larkBeta3 = median(beta3.v), larkSigma = median(tlsig.v))

#mash together
summOnlyPars <- summOnlyPars %>% 
  select(-alpha, -beta0, -sigma, -larkAlpha, -larkBeta0, -larkBeta1, -larkBeta2,
         -larkBeta3, -larkSigma) %>% 
  left_join(rickFrPar, by = "stkName", suffix = c("", ".1")) %>% 
  left_join(larkFrPar, by = "stkName", suffix = c("", ".1"))
#reorder
summOnlyPars <- summOnlyPars[ , c(1:8, 26:34, 9:25)]

# write.csv(summOnlyPars, here("data", "summOnlyCUPars.csv"), row.names = FALSE)


## Finally trim the MCMC parameters so they're in an appropriate format to pass
# directly to recoverySim
stkId <- cuPar %>%
  select(stkName, stk)
  
ricRecursiveOut <- rickFrRecursive %>% 
  mutate(stkName = stk.name) %>% 
  left_join(stkId, by = "stkName") %>% 
  select(stk, alpha = last.alpha.v, beta0 = beta.v, sigma = tlsig.v) %>% 
  arrange(stk)
write.csv(ricRecursiveOut, here("data", "trimRecursiveRickerMCMCPars.csv"), 
          row.names = FALSE)

larkRecursiveOut <- larkFrRecursive %>% 
  mutate(stkName = stk.name) %>% 
  left_join(stkId, by = "stkName") %>% 
  select(stk, alpha = last.alpha.v, beta0 = beta0.v, 
         beta1 = beta1.v, beta2 = beta2.v, beta3 = beta3.v, 
         sigma = tlsig.v) %>% 
  arrange(stk)
write.csv(larkRecursiveOut, here("data", "trimRecursiveLarkinMCMCPars.csv"),
          row.names = FALSE)

