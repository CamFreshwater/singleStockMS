#*******************************************************************************
# frStockRecModels_Constrained.r
# Date last revised: 2 April, 2019
# Modified version of frStockRecModels.R that runs stock-specific SR models with
# more severe constraints on capacity (top section)
# For beta use a cap on the capacity prior equal to maximum observed spawner
# abundance (assumes that capacity hasn't increased which seems reasonable)
# Based on FRSSI base run uses lognormal prior for capacity and 0.001 for tau 
# (neither parameter value had strong imapct on estimates compared to cap;
# compared in detail at bottom of this script)
#*******************************************************************************

require(here); require(R2jags); require(ggplot2); require(tidyverse)
require(samSim)

source(here("R/func/stockRecFunctions.R"))

# Clean data to match format for stockRecFunctions.R
# Use total spawners
srDatTrim <- read.csv(here("data/fraserRecDat.csv"), 
                  stringsAsFactors=F) %>% 
  select(stk, stkName, yr, totalSpwn, rec) %>% 
  rename(CU = stk, Year = yr, Escape = totalSpwn, Recruit = rec) %>% 
  filter(!(stkName == "Cultus" & Year > 1999)) %>% # trim Cultus
  mutate(lnRS = log(Recruit/Escape))
cuPar <- read.csv(here("data/fraserCUPars.csv"), stringsAsFactors=F)


# Fit models
ricFName = "fraserRicker_pooledConstrained"
RBayesSingle(dat = srDatTrim, Niter=250000, Nthin=50, burnin=50000, 
             fname = ricFName, model = "ricker", betaPrior = "logN", cap = 1, 
             tauPrior = 0.001)
larkFName = "fraserLarkin_pooledConstrained"
RBayesSingle(dat = srDatTrim, Niter=250000, Nthin=50, burnin=50000, 
             fname = larkFName, model = "larkin", betaPrior = "logN", cap = 1, 
             tauPrior = 0.001)


modelSummRick <- readRDS(here("outputs/srAnalysis/data", 
                              paste(ricFName, "SimOutput.RDS", 
                                    sep = "")))
modelSummLark <- readRDS(here("outputs/srAnalysis/data", 
                              paste(larkFName, "SimOutput.RDS", 
                                    sep = "")))

pooledRicker <- read.csv(here("outputs/srAnalysis/data", 
                              paste(ricFName, "MCMCPars.csv", sep = ""))) %>% 
  select(CU, alpha, beta0, sigma, deviance) %>% 
  mutate(deviance = -1 * deviance) %>% 
  rename(stk = CU)

pooledLarkin <- read.csv(here("outputs/srAnalysis/data", 
                              paste(larkFName, "MCMCPars.csv", sep = ""))) %>% 
  select(CU, alpha, beta0, beta1, beta2, beta3, sigma, deviance) %>% 
  mutate(deviance = -1 * deviance) %>% 
  rename(stk = CU)

## Look at medians to check if they seem reasonable
medRick <- pooledRicker %>% 
  group_by(stk) %>% 
  summarize(medA = median(alpha),
            medB = median(beta0),
            medSig = median(sigma))
medLark <- pooledLarkin %>% 
  group_by(stk) %>% 
  summarize(medA = median(alpha),
            medB = median(beta0),
            medB1 = median(beta1),
            medB2 = median(beta2),
            medB3 = median(beta3),
            medSig = median(sigma))

write.csv(pooledRicker, here("data", "constrainedRickerMCMCPars.csv"),
          row.names = FALSE)
write.csv(pooledLarkin, here("data", "constrainedLarkinMCMCPars.csv"),
          row.names = FALSE)

# --------
# Replace median values in FraserCUPars.csv with new MCMC samples and save new 
# copy
newCuPar <- read.csv(here("data/fraserCUpars.csv"), stringsAsFactors = F) %>% 
  mutate(alpha = round(medRick$medA, 2),
         beta0 = round(medRick$medB, 2),
         sigma = round(medRick$medSig, 2),
         larkAlpha = round(medLark$medA, 2),
         larkBeta0 = round(medLark$medB, 2),
         larkBeta1 = round(medLark$medB1, 2),
         larkBeta2 = round(medLark$medB2, 2),
         larkBeta3 = round(medLark$medB3, 2),
         larkSigma = round(medLark$medSig, 2))
write.csv(newCuPar, here("data", "fraserCUParsNEW.csv"),
          row.names = FALSE)

# --------

### Compare different priors with subset of stocks 
ricStks <- cuPar %>% 
  filter(model == "ricker") %>% 
  select(stkName) %>% 
  as.matrix()

trimSRDat <- srDat %>% 
  filter(stkName %in% ricStks) %>%
  mutate(lnRS = log(Recruit / Escape))

# betaPriors <- rep(c("logN", "unif"), each = 3)
caps <- rep(c(1, 1.5, 3), times = 2)
tauPriors <- rep(c(0.01, 0.001), each = 3)
fileNames <- paste("tau", tauPriors, "_cap", caps, sep = "")

# Run model
# for (i in seq_along(betaPriors)) {
#   RBayesSingle(dat = trimSRDat, Niter = 250000, Nthin = 50, burnin = 50000, 
#                fname = fileNames[i], model = "ricker", betaPrior = "logN", 
#                cap = caps[i], tauPrior = tauPriors[i])
# }

parList <- lapply(seq_along(fileNames), function(i) {
  pars <- read.csv(here("outputs/srAnalysis/data", 
                        paste(fileNames[i],"MCMCPars.csv", sep = "")))
  modelSumm <- readRDS(here("outputs/srAnalysis/data", 
                            paste(fileNames[i],"SimOutput.RDS", sep = "")))
  pars <- pars %>% 
    mutate(tauPriors = as.factor(tauPriors[i]),
           cap = as.factor(caps[i]),
           file = as.factor(fileNames[i]),
           stk = as.factor(abbreviate(stk, minlength = 4)))
  tList <- list(pars, modelSumm)
  names(tList) <- c("pars", "modelSummary")
  return(tList)
})

wideParDat <- do.call(rbind, lapply(parList, function(x) x$pars))

parDat <- wideParDat %>% 
  gather(key = parameter, value = value, -CU, -stk, -tauPriors, -file, -cap)
# write.csv(parDat, here("outputs", "generatedData", 
#                        "ricker_betaPriors_mcmcOut.csv"))

betaDat <- parDat %>% 
  filter(parameter == "beta0")

## Compare medians for different treatments
betaDat %>% 
  filter(tauPriors == "0.001",
         !cap == "1.5") %>% 
  group_by(stk, cap) %>% 
  summarize(medBeta = median(value)) %>% 
  group_by(stk) %>% 
  summarize(relDiff =  (max(medBeta) - min(medBeta)) / max(medBeta))
# differences range from 1.5 to 28% decrease with tighter constraint

png(here("outputs/srAnalysis/betaPosterior_capsTauPriors.png"), units = "in", 
    height = 4.5, width = 6.5, res = 250)
ggplot(betaDat, aes(x = tauPriors, y = value, fill = cap)) +
  labs(x = "Prior on Tau", y = "Beta Estimate", 
       fill = "Cap (Scalar\non Escapement)") +
  geom_boxplot(position = position_dodge(width = 0.65)) +
  theme_sleekX() +
  facet_wrap(~stk, scales = "free_y")
dev.off()

alphaDat <- parDat %>% 
  filter(parameter == "alpha")
png(here("outputs/srAnalysis/alphaPosterior_capsTauPriors.png"),  units = "in", 
    height = 4.5, width = 6.5, res = 250)
ggplot(alphaDat, aes(x = tauPriors, y = value, fill = cap)) +
  labs(x = "Prior on Tau", y = "Alpha Estimate", 
       fill = "Cap (Scalar\non Escapement)") +
  geom_boxplot(position = position_dodge(width = 0.65)) +
  theme_sleekX() +
  facet_wrap(~stk, scales = "free_y")
dev.off()



### Explore different options w/ Chilko
chk <- trimSRDat %>% 
  filter(stkName == "Chilko")
chkPars <- parDat %>% 
  filter(stk == "Chlk",
         tauPriors == "0.001") %>% 
  group_by(parameter) %>% 
  summarize(med = median(value))

plot(chk$Recruit ~ chk$Escape, pch = 21)
curve(x * exp(1.85 - 1.29 * x), from = 0, to = max(chk$Escape, na.rm = TRUE),
      add = TRUE, lwd = 2, col = "black")
curve(x * exp(1.85 - (0.75 * 1.29) * x), from = 0, 
      to = max(chk$Escape, na.rm = TRUE), add = TRUE, lwd = 2, col = "red")
curve(x * exp(1.85 - (1.25 * 1.29) * x), from = 0, 
      to = max(chk$Escape, na.rm = TRUE), add = TRUE, lwd = 2, col = "blue")

