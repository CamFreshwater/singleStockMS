#*******************************************************************************
# frStockRecModels.r
# Date last revised: 6 Feb, 2019
# Originally this script was generated to fit Ricker and Larkin models to pooled 
# early Shuswap stocks (i.e. Seymour and Scotch). Necessary because we are 
# interested in trends in CUs, but they are treated independently by management 
# and FRSSI. 
# In early Feb AMH noted that her earlier estimates of SR parameters for Cultus 
# were incorrect and should be adjusted (post-1999 data removed). To streamline 
# analysis all CUs were refit using the newly coded JAGS models after confirming
# they produced equivalent parameter estimates to AMH's originals.
#*******************************************************************************

require(here); require(R2jags); require(ggplot2); require(rstan); require(dplyr)

source(here("scripts/func/stockRecFunctions.R"))

# Clean data to match format for stockRecFunctions.R
# Use total
srDat <- read.csv(here("data/fraserDat/camGenDat/fraserRecDat.csv"), 
                  stringsAsFactors=F) %>% 
  select(stk, stkName, yr, totalSpwn, rec) %>% 
  rename(CU = stk, Year = yr, Escape = totalSpwn, Recruit = rec)
rickerPars <- read.csv(here("data/fraserDat/rickerMCMCPars.csv"), 
                       stringsAsFactors=F) 
larkinPars <- read.csv(here("data/fraserDat/larkinMCMCPars.csv"), 
                       stringsAsFactors=F) 


### Compare the JAGS models to FRSSI par estimates to ensure that models
## produce equivalent results
trimDat <- srDat %>% 
  filter(CU %in% c("1", "2"),
         !Year == "2011")  %>% 
  mutate(lnRS = log(Recruit / Escape))
RBayesSingle(dat = trimDat, Niter=250000, Nthin=50, burnin=50000, fname="TEST", 
             model = "larkin")

newLarkin <- read.csv(here("outputs/srAnalysis/data", 
                          paste("TEST","MCMCPars.csv", sep = "")))
oldLarkin <- larkinPars %>% 
  filter(stk %in% c(1,2))

stkId <- unique(newLarkin$CU)
par(mfrow = c(2, 2))
for (i in 1:2) {
  newDat <- newLarkin %>% 
    filter(CU == stkId[i])
  oldDat <- oldLarkin %>% 
    filter(stk == stkId[i])  
  hist(newDat$betas, col = alpha("red", 0.5), main = stkId[i])
  hist(oldDat$beta0, col = alpha("blue", 0.5), add=T)
  hist(newDat$alphas, col = alpha("red", 0.5), main = stkId[i])
  hist(oldDat$alpha, col = alpha("blue", 0.5), add=T)
}
# Larkin fits well

trimRicDat <- srDat %>% 
  filter(CU %in% c("3", "4"),
         !Year == "2011") %>% 
  mutate(lnRS = log(Recruit / Escape))
RBayesSingle(dat = trimRicDat, Niter=250000, Nthin=50, burnin=50000, 
             fname="TESTR", model = "ricker")

newRicker <- read.csv(here("outputs/srAnalysis/data", 
                           paste("TESTR","MCMCPars.csv", sep = "")))
oldRicker <- rickerPars %>% 
  filter(stk %in% c(3, 4)) 

oldRicker %>% 
  group_by(stk) %>% 
  summarize(medB = median(beta0), medA = median(alpha))
newRicker %>% 
  group_by(CU) %>% 
  summarize(medB = median(beta0), medA = median(alpha))

stkId <- unique(newRicker$CU)
par(mfrow = c(2, 2))
for (i in 1:2) {
  newDat <- newRicker %>% 
    filter(CU == stkId[i])
  oldDat <- oldRicker %>% 
    filter(stk == stkId[i])  
  hist(newDat$beta0, col = alpha("red", 0.5), main = stkId[i])
  hist(oldDat$beta0, col = alpha("blue", 0.5), add=T)
  hist(newDat$alpha, col = alpha("red", 0.5), main = stkId[i])
  hist(oldDat$alpha, col = alpha("blue", 0.5), add=T)
}
#Ricker fits well


### Fit all CUs with full dataset (i.e. don't exclude 2010 even though it 
## wasn't included in AMH's analysis); exception is Cultus
srDatTrim <- srDat %>% 
  filter(!(stkName == "Cultus" & Year > 1999)) %>% 
  mutate(lnRS = log(Recruit/Escape))

RBayesSingle(dat = srDatTrim, Niter=250000, Nthin=50, burnin=50000, 
             fname="fraserPooledRicker", model = "ricker")
RBayesSingle(dat = srDatTrim, Niter=250000, Nthin=50, burnin=50000, 
             fname="fraserPooledLarkin", model = "larkin")


modelSummRick <- readRDS(here("outputs/srAnalysis/data", 
                              paste("fraserPooledRicker", "SimOutput.RDS", 
                                    sep = "")))
modelSummLark <- readRDS(here("outputs/srAnalysis/data", 
                              paste("fraserPooledLarkin", "SimOutput.RDS", 
                                    sep = "")))

pooledRicker <- read.csv(here("outputs/srAnalysis/data", 
                              paste("fraserPooledRicker", "MCMCPars.csv", 
                                    sep = ""))) %>% 
  select(CU, alpha, beta0, sigma, deviance) %>% 
  mutate(deviance = -1 * deviance) %>% 
  rename(stk = CU)

pooledLarkin <- read.csv(here("outputs/srAnalysis/data", 
                          paste("fraserPooledLarkin", "MCMCPars.csv", 
                                sep = ""))) %>% 
  select(CU, alpha, beta0, beta1, beta2, beta3, sigma, deviance) %>% 
  mutate(deviance = -1 * deviance) %>% 
  rename(stk = CU)


## Look at medians to check if they seem reasonable
pooledRicker %>% 
  group_by(stk) %>% 
  summarize(medA = median(alpha),
            medB = median(beta0))
pooledLarkin %>% 
  group_by(stk) %>% 
  summarize(medA = median(alpha),
            medB = median(beta0))

write.csv(pooledRicker,
          here("data/fraserDat/camGenDat", "pooledRickerMCMCPars.csv"),
          row.names = FALSE)
write.csv(pooledLarkin,
          here("data/fraserDat/camGenDat", "pooledLarkinMCMCPars.csv"),
          row.names = FALSE)


#### OLD CODE USED TO COMBINE ESHU AND ORIGINAL FILES

### Combine Scotch/Seymour and fit model
# srDatPool <- srDat %>% 
#   filter(Year > 1979, #crop at first year of shorter TS
#          CU == "815") %>% 
#   group_by(Year, stkName, CU) %>% 
#   summarize(Escape = sum(Escape), Recruit = sum(Recruit)) %>% #combine R and S
#   mutate(CU = as.character(CU),
#          lnRS = log(Recruit/Escape)) #recalculate prod
# 
# # Does it look like there is evidence of cycles?
# plot(Recruit ~ Year, type = "l", data = srDatPool)
# 
# RBayesSingle(dat = srDatPool, Niter=250000, Nthin=50, burnin=50000, 
#              fname="earlyShuRicker", model = "ricker")
# RBayesSingle(dat = srDatPool, Niter=250000, Nthin=50, burnin=50000, 
#              fname="earlyShuLarkin", model = "larkin")
# 
# modelSummRick <- readRDS(here("outputs/srAnalysis/data", 
#                               paste("earlyShuRicker","SimOutput.RDS", sep = "")))
# modelSummLark <- readRDS(here("outputs/srAnalysis/data", 
#                               paste("earlyShuLarkin","SimOutput.RDS", 
#                                     sep = "")))
# head(modelSummRick) #DIC = 89.2
# head(modelSummLark) #DIC = 77
# 
# 
# ## Add new MCMC estimates to FRSSI files
# eShuRick <- read.csv(here("outputs/srAnalysis/data", 
#                           paste("earlyShuRicker","MCMCPars.csv", sep = ""))) %>% 
#   select(CU, alpha, beta0, sigma, deviance) %>% 
#   mutate(deviance = -1 * deviance) %>% 
#   rename(stk = CU)
# eShuLark <- read.csv(here("outputs/srAnalysis/data", 
#                           paste("earlyShuLarkin","MCMCPars.csv", sep = ""))) %>% 
#   select(CU, alpha, beta0, beta1, beta2, beta3, sigma, deviance) %>% 
#   mutate(deviance = -1 * deviance) %>% 
#   rename(stk = CU)
# 
# updatedRicker <- rbind(rickerPars, eShuRick)
# updatedLarkin <- rbind(larkinPars, eShuLark)