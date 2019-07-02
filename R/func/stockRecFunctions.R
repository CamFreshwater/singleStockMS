## Functions to estimate stock recruit pararameters as inputs
## in closed-loop simulation model
## Originally formulated for South Coast chum assessment by Brooke
## Davis starting July 20, 2015; minor edits by Cam Freshwater
## starting Feb 22, 2018, predominantly to file paths within function:
#### Moded: RBayesHier and FitRickerB
### Removed: percentile calcs, BM calcs, and all SR models except hierarchical 
### and basic varying pars
### Note: FitRickerB only partially modified; if used match structure RBayesHier
### Note March 16: hierarchical model makes call to BayesHierOneErr which has obs
### tau set to 0 to be consistent w/ Fraser models


library("dplyr"); library("parallel"); library("doParallel"); library("foreach")



##################
##  RBayesHier  ##
##################

## Mutlistock Bayes with simple hierarchy
## alpha does NOT vary over time
## each CU's alpha value chosen from a global distribution
## NOTE: modified from B. Davis's original script to use JAGS
## model MultiStockHierSimple.txt (i.e. lacks process error variation) and to 
## pass log(R/S), rather than R, as input variable
# dat=trimDat;
# plotDiags=T; writeDat=T; saveMod=T; Niter=125000; Nthin=50; burnin=25000;
# fname="TEST"; Fyear=F; threshSMSY=0.8; SmaxNorm=F; CVcap=0.3; CVmean=F; Rhat=T;
# mult=1.5; savedCores=1; PSF=TRUE

RBayesHier <- function(dat = SRdat, Niter = 75000, Nthin = 50, burnin = 25000, 
                       fname, MuA = 1, TauA = 0.001, mult = 1.5, model = "ar"){
  
  ifelse(!dir.exists(paste(here(),"/outputs/srAnalysis/figs/", sep="")), dir.create(paste(here(),"/outputs/srAnalysis/figs/", sep="")), FALSE)
  ifelse(!dir.exists(paste(here(),"/outputs/srAnalysis/data/", sep="")), dir.create(paste(here(),"/outputs/srAnalysis/data/", sep="")), FALSE)
  
  # prepare data to feed into model
  # number of CUs
  CUs <- unique(dat$CU)
  nCUs <- length(CUs)
  
  inits <- function(){
    list(tau_R=1, MuA_global=1, TauA_global=3, alpha=rep(1,nsites) )
  }
  
  #indices of non NA data points for R/S and for just S
  inds <- which(is.na(dat$lnRS)==F & is.na(dat$Escape)==F)
  inds2 <- which(is.na(dat$Escape)==F)
  # Recruits and spawners data
  dataGood <- dat[inds,] 
  lnRSdat <- dataGood$lnRS
  Sdat <- dataGood$Escape
  escGood <- dat[inds2,] 
  # Years and CUs corresponding to data points, as indices
  # need to change so that missing years are accomodated, differing numbers of years for CUs 
  sampledCUs <- NULL
  for(i in 1:length(inds)){
    sampledCUs[i] <- which(CUs==dataGood$CU[i])
  }
  # need upper bound for Smax prior for each site, number of years each site
  upp <- NULL
  Years <- NULL
  nyears <- NULL
  MuSmax <- NULL
  TauSmax <- NULL
  for(i in 1:length(CUs)){
    #use data frame from ricker fit to use later to initialize model
    cuDat <- dataGood[which(dataGood$CU==CUs[i]),]
    EscDat <- escGood[which(escGood$CU==CUs[i]),]
    #calculate upper bound on Smax prior
    upp[i] <- max(EscDat$Escape)*mult
    nyears[i] <- dim(cuDat)[1]
    Years <- c(Years, 1:nyears[i])
  }    
  
  # convert obs spawners and recruits to matrices
  Smat <- matrix(NA, nrow=nCUs, ncol=max(nyears))
  lnRSmat <- matrix(NA, nrow=nCUs, ncol=max(nyears))
  for(k in 1:length(inds)){
    Smat[sampledCUs[k], Years[k]] <- Sdat[k]
    lnRSmat[sampledCUs[k], Years[k]] <- lnRSdat[k]
  }
  
  
  # create data input for JAGS
  if (model == "singleLarkin"){
    Dlist <- list ("Smat"=Smat, "lnRSmat"=lnRSmat, "nyears"=nyears, "low"=1, 
                   "upp"=upp, "nCUs"=nCUs)
    simoutJAGS <- jags(Dlist, model.file = "C:/github/salmon-sim/scripts/jagsModels/MultiStockHierSimple.txt",
                       parameters.to.save=c("alpha", "beta", "beta1", "beta2", 
                                            "beta3", "tau_R", "sig_R"), 
                       n.chains = 3, n.iter = Niter, n.thin=Nthin, n.burnin=burnin)
  }
  if (model == "standard"){
    Dlist <- list ("Smat"=Smat, "lnRSmat"=lnRSmat, "nyears"=nyears, "MuA"= MuA,
                   "TauA" = TauA, "low"=1, "upp"=upp, "nCUs"=nCUs)
    simoutJAGS <- jags(Dlist, model.file = "C:/github/salmon-sim/scripts/jagsModels/MultiStockHierSimple.txt",
                       parameters.to.save=c("alpha", "beta", "MuA_global", "TauA_global", 
                                            "tau_R"), 
                       n.chains = 3, n.iter = Niter, n.thin=Nthin, n.burnin=burnin)
  }
  if (model == "ar") {
    Dlist <- list("Smat"=Smat, "lnRSmat"=lnRSmat, "nyears"=nyears, "MuA"= MuA,
                  "TauA" = TauA, "low"=1, "upp"=upp, "nCUs"=nCUs)
    simoutJAGS <- jags(Dlist, model.file = "C:/github/salmon-sim/scripts/jagsModels/unifMultiStockHierAR.txt",
                       parameters.to.save=c("alpha", "beta", "rho", "MuA_global", "TauA_global", 
                                            "tau_R"), 
                       n.chains = 3, n.iter = Niter, n.thin=Nthin, n.burnin=burnin)
  }
  if(model == "psf"){
    prCV <- c(10, 10, 10, 1, 10, 1) # vector of priors for PSF model provided by E. Hertz plus priors added for CU 2 (which they dropped)
    prSmax <- c(24524, 25000, 42230, 1118, 20378, 2544)# vector of priors for PSF model provided by E. Hertz
    ### First snippet uses PSF's code which is incompatible with most of post-processing code
    Dlist <- list ("Smat"=Smat, "lnRSmat"=lnRSmat, "nyears"=nyears, "nCUs"=nCUs,
                   "MuA" = MuA, "TauA" = TauA, "prCV"=prCV, "prSmax"=prSmax)
    simoutJAGS <- jags(Dlist, model.file = "C:/github/salmon-sim/scripts/jagsModels/MultiStockHierPSFRep.txt",
                       parameters.to.save=c("alpha", "beta", "MuA_global", "sigA_global",
                                            "TauA_global", "tau_R", "sd"),
                       n.chains = 3, n.iter = Niter, n.thin=Nthin, n.burnin=burnin)
  }
  if(model == "altBeta"){
    prCV <- c(10, 10, 10, 1, 10, 1) # vector of priors for PSF model provided by E. Hertz plus priors added for CU 2 (which they dropped)
    prSmax <- c(24524, 25000, 42230, 1118, 20378, 2544)# vector of priors for PSF model provided by E. Hertz
    ### First snippet uses PSF's code which is incompatible with most of post-processing code
    Dlist <- list ("Smat"=Smat, "lnRSmat"=lnRSmat, "nyears"=nyears, "nCUs"=nCUs,
                   "MuA" = MuA, "TauA" = TauA, "prCV"=prCV, "prSmax"=prSmax)
    simoutJAGS <- jags(Dlist, model.file = "C:/github/salmon-sim/scripts/jagsModels/MultiStockHierSimple_altBeta.txt",
                       parameters.to.save=c("alpha", "beta", "MuA_global", "sigA_global",
                                            "TauA_global", "tau_R", "sd"),
                       n.chains = 3, n.iter = Niter, n.thin=Nthin, n.burnin=burnin)
  } 
  
  # change to bugs output to easily extract values
  simoutBUGS <- simoutJAGS$BUGSoutput
  
  # extract mcmc list to plot usuing coda functions
  ToPlot <- as.mcmc(simoutJAGS)
  # plot traces and posterior densities
  pdf( paste(here(),"/outputs/srAnalysis/figs/",fname,"TraceDens.pdf", sep="") )
  plot(ToPlot)
  dev.off()
  # plot gelman plots
  pdf( paste(here(),"/outputs/srAnalysis/figs/",fname,"Gelman.pdf", sep="" ))
  tryCatch( {
    gelman.plot(ToPlot)
  }, error=function(x=1){plot.new() })
  dev.off()
  
  
  # Re-arrange and store each iteration's parameter estimates
  if(nrow(simoutBUGS$sims.matrix) > 1000){
    simSamples <- sample_n(as.data.frame(simoutBUGS$sims.matrix), size=1000) #subset to 1000 samples equivalent to Fraser outputs generated by FRSSI modelers
  } else{
    simSamples <- simoutBUGS$sims.matrix
  }
  simPars <- NULL
  for(i in seq_along(CUs)){
    CU <- as.character(rep(CUs[i], each=nrow(simSamples)))
    stk <- i
    muAlpha <- simSamples[,"MuA_global"]
    tauAlpha <- simSamples[,"TauA_global"]
    alphas <- simSamples[,paste('alpha', "[", i , "]", sep="")]
    betas <- simSamples[,paste('beta', "[", i , "]", sep="")]
    tauRs <- simSamples[,paste('tau_R', "[", i , "]", sep="")]
    allParsTemp <- cbind(stk, CU, muAlpha, tauAlpha, alphas, betas, tauRs)
    if (model == "ar") {
      rhos <- simSamples[,paste('rho', "[", i , "]", sep="")]
      allParsTemp <- cbind(allParsTemp, rhos)
    }
    simPars <- rbind(simPars, allParsTemp)
  }

  #save output for later
  saveRDS(simoutJAGS, file=paste(here(),"/outputs/srAnalysis/data/",fname,"SimOutput",".rds", sep=""))
  write.csv(simPars, file=paste(here(),"/outputs/srAnalysis/data/",fname,"MCMCPars.csv", sep=""), row.names=FALSE)
  
  # save Gelman - Rubin Statistics, Geweke; need readable data frame with mod, year etc
  RhatVec <- simoutBUGS$summary[,"Rhat"]
  GewekeList <- geweke.diag(as.mcmc(simoutJAGS))
  GewekeTab <- data.frame(Param = names(GewekeList[[1]]$z) ,G1=GewekeList[[1]]$z, G2=GewekeList[[2]]$z, G3=GewekeList[[3]]$z )
  RhatTab <- data.frame(Mod=rep("RBHier", length(RhatVec)), Year = rep(max(dat$Year), length(RhatVec)), Param = names(RhatVec), 
                        Rhat=RhatVec, row.names=NULL)
  Rhats <- inner_join(RhatTab, GewekeTab, by="Param")
  CU <- NULL
  for(i in 1:dim(Rhats)[1]){
    Separated <- unlist(strsplit(as.character(Rhats$Param[i]), split="\\[|\\]"))
    # if there is a second item it will be the number, referring to the CU
    if(length(Separated)==2){
      CU[i] <- as.character(CUs[as.numeric(Separated[2])])
    } else { CU[i] <- unlist(strsplit(gsub("([0-9]+)","~\\1", fname), "~"))[1] }
  }
  Rhats$CU <- as.character(CU)
}
#___________________________________________________________________________________________________________


## Single stock Bayesian Ricker or Larkin models
# Although Fraser doesn't include gappy data left that functionality in just in
# case
# dat=srDatTrim
# Niter=12500; Nthin=50; burnin=2500
# fname="TEST"; savedCores=1; model = "larkin"
# cap = 1; betaPrior = "logN"; tauPrior = 0.001

RBayesSingle <- function(dat = SRdat, Niter = 75000, Nthin = 50, burnin = 25000,
                         fname, model = "ricker", betaPrior = "logN", 
                         cap = 3, tauPrior = 0.001) {
  
  ifelse(!dir.exists(paste(here(),"/outputs/srAnalysis/figs/", sep="")), 
         dir.create(paste(here(),"/outputs/srAnalysis/figs/", sep="")), FALSE)
  ifelse(!dir.exists(paste(here(),"/outputs/srAnalysis/data/", sep="")), 
         dir.create(paste(here(),"/outputs/srAnalysis/data/", sep="")), FALSE)
  
  # prepare data to feed into model
  # number of CUs
  stkNames <- unique(dat$stkName)
  CUs <- unique(dat$CU)
  nCUs <- length(CUs)
  
  inits <- function(){
    list(tau_R = 1, alpha=rep(1, nsites))
  }
  
  #indices of non NA data points for R/S and for just S
  inds <- which(is.na(dat$lnRS)==F & is.na(dat$Escape)==F)
  inds2 <- which(is.na(dat$Escape)==F)
  # Recruits and spawners data
  dataGood <- dat[inds,] 
  lnRSdat <- dataGood$lnRS
  Sdat <- dataGood$Escape
  escGood <- dat[inds2,] 
  # Years and CUs corresponding to data points, as indices
  # need to change so that missing years are accomodated, differing numbers of years for CUs 
  sampledCUs <- NULL
  for(i in 1:length(inds)){
    sampledCUs[i] <- which(CUs==dataGood$CU[i])
  }
  # need upper bound for Smax prior for each site, number of years each site
  upp <- NULL
  Years <- NULL
  nyears <- NULL
  MuSmax <- NULL
  TauSmax <- NULL
  for(i in 1:length(CUs)){
    #use data frame from ricker fit to use later to initialize model
    cuDat <- dataGood[which(dataGood$CU==CUs[i]), ]
    EscDat <- escGood[which(escGood$CU==CUs[i]), ]
    #calculate upper bound on Smax prior
    upp[i] <- max(EscDat$Escape) #
    nyears[i] <- dim(cuDat)[1]
    Years <- c(Years, 1:nyears[i])
  }    
  
  # convert obs spawners and recruits to matrices
  Smat <- matrix(NA, nrow=nCUs, ncol=max(nyears))
  lnRSmat <- matrix(NA, nrow=nCUs, ncol=max(nyears))
  for(k in 1:length(inds)){
    Smat[sampledCUs[k], Years[k]] <- Sdat[k]
    lnRSmat[sampledCUs[k], Years[k]] <- lnRSdat[k]
  }
  
  # create data input for JAGS
  Dlist <- list("Smat" = Smat, "lnRSmat" = lnRSmat, "nyears" = nyears, 
                "upp" = upp, "nCUs" = nCUs, "cap" = cap, "tauPrior" =  tauPrior)
  if (model == "larkin") {
    if (betaPrior != "unif") {
      simoutJAGS <- jags(Dlist,
                         model.file = here("R", "jagsModels", 
                                           "singleStockLarkin_logN.txt"),
                         parameters.to.save=c("alpha", "beta", "beta1", "beta2", 
                                              "beta3", "tau_R", "sig_R"), 
                         n.chains = 3, n.iter = Niter, n.thin=Nthin, 
                         n.burnin=burnin)
    } else {
      simoutJAGS <- jags(Dlist,
                         model.file = here("R", "jagsModels", 
                                           "singleStockLarkin_unif.txt"),
                         parameters.to.save=c("alpha", "beta", "beta1", "beta2", 
                                              "beta3","tau_R", "sig_R"), 
                         n.chains = 3, n.iter = Niter, n.thin=Nthin, 
                         n.burnin=burnin)
    }
  } #end if model = Larkin
  if (model == "ricker") {
    if (betaPrior != "unif") {
      simoutJAGS <- jags(Dlist,
                         model.file = here("R", "jagsModels", 
                                         "singleStockRicker_logN.txt"),
                         parameters.to.save=c("alpha", "beta", "tau_R", "sig_R"), 
                         n.chains = 3, n.iter = Niter, n.thin=Nthin, 
                         n.burnin=burnin)
    } else {
      simoutJAGS <- jags(Dlist,
                         model.file = here("R", "jagsModels", 
                                           "singleStockRicker_unif.txt"),
                         parameters.to.save=c("alpha", "beta", "tau_R", "sig_R"), 
                         n.chains = 3, n.iter = Niter, n.thin=Nthin, 
                         n.burnin=burnin)
    }
  } #end if model = Rciker
  # change to bugs output to easily extract values
  simoutBUGS <- simoutJAGS$BUGSoutput
  
  # extract mcmc list to plot usuing coda functions
  ToPlot <- as.mcmc(simoutJAGS)
  # plot traces and posterior densities
  pdf( paste(here(),"/outputs/srAnalysis/figs/",fname,"TraceDens.pdf", sep="") )
  plot(ToPlot)
  dev.off()
  # plot gelman plots
  pdf( paste(here(),"/outputs/srAnalysis/figs/",fname,"Gelman.pdf", sep="" ))
  tryCatch( {
    gelman.plot(ToPlot)
  }, error=function(x=1){plot.new() })
  dev.off()
  
  # Re-arrange and store each iteration's parameter estimates
  if(nrow(simoutBUGS$sims.matrix) > 1000){
    simSamples <- sample_n(as.data.frame(simoutBUGS$sims.matrix), size=1000) #subset to 1000 samples equivalent to Fraser outputs generated by FRSSI modelers
  } else{
    simSamples <- simoutBUGS$sims.matrix
  }
  #Change header for simSamples since default is different when nCU=1
  if (nCUs == 1) {
    colnames(simSamples) <- paste(colnames(simSamples), "[1]", sep="") 
    simSamples <- simSamples %>% 
      dplyr::rename(deviance = "deviance[1]")
  }
  simPars <- NULL
  for(i in seq_along(CUs)){
    CU <- rep(CUs[i], each=nrow(simSamples))
    alpha <- simSamples[, paste('alpha', "[", i , "]", sep="")]
    beta0 <- simSamples[,paste('beta', "[", i , "]", sep="")]
    sigma <- simSamples[,paste('sig_R', "[", i , "]", sep="")]
    deviance <- simSamples[, "deviance"]
    allParsTemp <- cbind(CU, alpha, beta0, sigma, deviance)
    if (model == "larkin") {
      beta1 <- simSamples[,paste('beta1', "[", i , "]", sep="")]
      beta2 <- simSamples[,paste('beta2', "[", i , "]", sep="")]
      beta3 <- simSamples[,paste('beta3', "[", i , "]", sep="")]
      allParsTemp <- cbind(allParsTemp, beta1, beta2, beta3)
    }
    allParsTemp <- allParsTemp %>% 
      as.data.frame() %>% 
      mutate(stk = stkNames[i])
    simPars <- rbind(simPars, allParsTemp)
  }
  
  #save output for later
  saveRDS(simoutJAGS, file=paste(here(),"/outputs/srAnalysis/data/",
                                 fname,"SimOutput",".rds", sep=""))
  write.csv(simPars, file=paste(here(),"/outputs/srAnalysis/data/",
                                fname,"MCMCPars.csv", sep=""), row.names=FALSE)
}
#___________________________________________________________________________________________________________


#Function to sample vector with replacement then calculate median
genMedian <- function(dat, lengthOut){
  output <- rep(NA, length.out = lengthOut)
  for(i in 1:lengthOut){
    output[i] <- median(sample(dat, length(dat), replace=TRUE))
  }
  return(output)
}
#___________________________________________________________________________________________________________

# Equivalent to RBayesHier, but includes an AR term (necessary given evidence of autocorrelation in north coast chum)
RBayesAR <- function(dat=SRdat, plotDiags=T, writeDat=T, saveMod=T, Niter=75000, Nthin=50, burnin=25000, fname, Fyear=F, threshSMSY=0.8,
                       SmaxNorm=F, CVcap=0.3, CVmean=CVmean, Rhat=F, mult=1.5, savedCores=1,
                     varDist="gamma"){
  
  
  ifelse(!dir.exists(paste(here(),"/outputs/srAnalysis/figs/", sep="")), dir.create(paste(here(),"/outputs/srAnalysis/figs/", sep="")), FALSE)
  ifelse(!dir.exists(paste(here(),"/outputs/srAnalysis/data/", sep="")), dir.create(paste(here(),"/outputs/srAnalysis/data/", sep="")), FALSE)
  
  # prepare data to feed into model
  # number of CUs
  CUs <- unique(dat$CU)
  nCUs <- length(CUs)
  
  # inits <- function(){
  #   list(tau_R=1, MuA_global=1, TauA_global=3, alpha=rep(1,nsites))
  # }
  
  #indices of non NA data points for R/S and for just S
  inds <- which(is.na(dat$lnRS)==F & is.na(dat$Escape)==F)
  inds2 <- which(is.na(dat$Escape)==F)
  # Recruits and spawners data
  dataGood <- dat[inds,] 
  lnRSdat <- dataGood$lnRS
  Sdat <- dataGood$Escape
  escGood <- dat[inds2,] 
  # Years and CUs corresponding to data points, as indices
  # need to change so that missing years are accomodated, differing numbers of years for CUs 
  sampledCUs <- NULL
  for(i in 1:length(inds)){
    sampledCUs[i] <- which(CUs==dataGood$CU[i])
  }
  # need upper bound for Smax prior for each site, number of years each site
  upp<-NULL
  Years <- NULL
  nyears <- NULL
  MuSmax <- NULL
  TauSmax <- NULL
  for(i in 1:length(CUs)){
    #use data frame from ricker fit to use later to initialize model
    cuDat <- dataGood[which(dataGood$CU==CUs[i]),]
    EscDat <- escGood[which(escGood$CU==CUs[i]),]
    #calculate upper bound on Smax prior
    upp[i] <- max(EscDat$Escape)*mult
    nyears[i] <- dim(cuDat)[1]
    Years <- c(Years, 1:nyears[i])
    if( SmaxNorm==T){
      MuSmax[i] <- mean(EscDat$Escape)
      TauSmax[i] <- 1/(log(CVmean^2 + 1))
    }
  }    
  
  # convert obs spawners and recruits to matrices
  Smat <- matrix(NA, nrow=nCUs, ncol=max(nyears))
  lnRSmat <- matrix(NA, nrow=nCUs, ncol=max(nyears))
  for(k in 1:length(inds)){
    Smat[sampledCUs[k], Years[k]] <- Sdat[k]
    lnRSmat[sampledCUs[k], Years[k]] <- lnRSdat[k]
  }
  
  # create data input for JAGS
  Dlist <- list("Smat"=Smat, "lnRSmat"=lnRSmat, "nyears"=nyears, "MuA"=1, "TauA"=0.001, 
                 "low"=1, "upp"=upp, "nCUs"=nCUs)
  if (varDist == "uniform") {
    simoutJAGS <- jags(Dlist, model.file = "C:/github/salmon-sim/scripts/jagsModels/unifMultiStockHierAR.txt",
                      parameters.to.save=c("alpha", "beta", "rho", "MuA_global", "TauA_global", 
                                           "tau_R"), 
                      n.chains = 3, n.iter = Niter, n.thin=Nthin, n.burnin=burnin)
  }                     
  if (varDist == "gamma") {
    simoutJAGS <- jags(Dlist, model.file = "C:/github/salmon-sim/scripts/jagsModels/gammaMultiStockHierAR.txt",
                       parameters.to.save=c("alpha", "beta", "rho", "MuA_global", "TauA_global", 
                                            "tau_R"), 
                       n.chains = 3, n.iter = Niter, n.thin=Nthin, n.burnin=burnin)
  }                     
  
  # change to bugs output to easily extract values
  simoutBUGS <- simoutJAGS$BUGSoutput

  # extract mcmc list to plot usuing coda functions
  if(plotDiags==T){
    ToPlot <- as.mcmc(simoutJAGS)
    # plot traces and posterior densities
    pdf( paste(here(),"/outputs/srAnalysis/figs/",fname,"TraceDens.pdf", sep="") )
    plot(ToPlot)
    dev.off()
    # plot gelman plots
    pdf( paste(here(),"/outputs/srAnalysis/figs/",fname,"Gelman.pdf", sep="" ))
    tryCatch( {
      gelman.plot(ToPlot)
    }, error=function(x=1){plot.new() })
    dev.off()
  }
  
  # Re-arrange and store each iteration's parameter estimates
  if(nrow(simoutBUGS$sims.matrix) > 1000){
    simSamples <- sample_n(as.data.frame(simoutBUGS$sims.matrix), size=1000) #subset to 1000 samples equivalent to Fraser outputs generated by FRSSI modelers
  } else{
    simSamples <- simoutBUGS$sims.matrix
  }
  simPars <- NULL
  for(i in seq_along(CUs)){
    CU <- as.character(rep(CUs[i], each=nrow(simSamples)))
    stk <- i
    muAlpha <- simSamples[,"MuA_global"]
    tauAlpha <- simSamples[,"TauA_global"]
    alphas <- simSamples[,paste('alpha', "[", i , "]", sep="")]
    betas <- simSamples[,paste('beta', "[", i , "]", sep="")]
    rhos <- simSamples[,paste('rho', "[", i , "]", sep="")]
    tauRs <- simSamples[,paste('tau_R', "[", i , "]", sep="")]
    allParsTemp <- cbind(stk, CU, muAlpha, tauAlpha, alphas, betas, rhos, tauRs)
    simPars <- rbind(simPars, allParsTemp)
  }
  
  #save output for later
  if(saveMod==T){
    saveRDS(simoutJAGS, file=paste(here(),"/outputs/srAnalysis/data/",fname,"SimOutput",".rds", sep=""))
    write.csv(simPars, file=paste(here(),"/outputs/srAnalysis/data/",fname,"MCMCPars.csv", sep=""), row.names=FALSE)
  }
  
  # 11/30/15 save Gelman - Rubin Statistics, Geweke; need readable data frame with mod, year etc
  RhatVec <- simoutBUGS$summary[,"Rhat"]
  GewekeList <- geweke.diag(as.mcmc(simoutJAGS))
  GewekeTab <- data.frame(Param = names(GewekeList[[1]]$z) ,G1=GewekeList[[1]]$z, G2=GewekeList[[2]]$z, G3=GewekeList[[3]]$z )
  RhatTab <- data.frame(Mod=rep("RBHier", length(RhatVec)), Year = rep(max(dat$Year), length(RhatVec)), Param = names(RhatVec), 
                        Rhat=RhatVec, row.names=NULL)
  Rhats <- inner_join(RhatTab, GewekeTab, by="Param")
  CU <- NULL
  for(i in 1:dim(Rhats)[1]){
    Separated <- unlist(strsplit(as.character(Rhats$Param[i]), split="\\[|\\]"))
    # if there is a second item it will be the number, referring to the CU
    if(length(Separated)==2){
      CU[i] <- as.character(CUs[as.numeric(Separated[2])])
    } else { CU[i] <- unlist(strsplit(gsub("([0-9]+)","~\\1", fname), "~"))[1] }
  }
  Rhats$CU <- as.character(CU)
}
#___________________________________________________________________________________________________________

# dat=chumDat; plotDiags=T; writeDat=T; saveMod=T; Niter=2000; Nthin=50; burnin=500;
# fname="ncChumStan"; Fyear=F; threshSMSY=0.8; SmaxNorm=F; CVcap=0.3; CVmean=F; Rhat=T;
# mult=1.5; savedCores=1


#Stan version of hierarchical SR model WORK IN PROGRESS
RBayesHierStan <- function(dat=SRdat, plotDiags=T, writeDat=T, saveMod=T, Niter=75000, Nthin=50, burnin=25000, fname, Fyear=F, threshSMSY=0.8,
                       SmaxNorm=F, CVcap=0.3, CVmean=CVmean, Rhat=F, mult=1.5, savedCores=1){
  
  
  ifelse(!dir.exists(paste(here(),"/outputs/srAnalysis/figs/", sep="")), dir.create(paste(here(),"/outputs/srAnalysis/figs/", sep="")), FALSE)
  ifelse(!dir.exists(paste(here(),"/outputs/srAnalysis/data/", sep="")), dir.create(paste(here(),"/outputs/srAnalysis/data/", sep="")), FALSE)
  
  
  
  #how many cores does the computer have?
  # Ncores <- detectCores()
  # cl <- makeCluster(Ncores-savedCores)
  # registerDoParallel(cl)
  
  # prepare data to feed into model
  # number of CUs
  CUs <- unique(dat$CU)
  nCUs <- length(CUs)
  
  inits <- function(){
    list(tau_R=1, MuA_global=1, TauA_global=3, alpha=rep(1,nsites) )
  }
  
  #indices of non NA data points for R/S and for just S
  inds <- which(is.na(dat$lnRS)==F & is.na(dat$Escape)==F)
  inds2 <- which(is.na(dat$Escape)==F)
  # Recruits and spawners data
  dataGood <- dat[inds,] 
  lnRSdat <- dataGood$lnRS
  Sdat <- dataGood$Escape
  escGood <- dat[inds2,] 
  # Years and CUs corresponding to data points, as indices
  # need to change so that missing years are accomodated, differing numbers of years for CUs 
  sampledCUs <- NULL
  for(i in 1:length(inds)){
    sampledCUs[i] <- which(CUs==dataGood$CU[i])
  }
  # need upper bound for Smax prior for each site, number of years each site
  upp<-NULL
  Years <- NULL
  nyears <- NULL
  MuSmax <- NULL
  TauSmax <- NULL
  for(i in 1:length(CUs)){
    #use data frame from ricker fit to use later to initialize model
    cuDat <- dataGood[which(dataGood$CU==CUs[i]),]
    EscDat <- escGood[which(escGood$CU==CUs[i]),]
    #calculate upper bound on Smax prior
    
    
    upp[i] <- max(EscDat$Escape)*mult
    # upp[i] <- mean(EscDat$Escape)
    
    nyears[i] <- dim(cuDat)[1]
    Years <- c(Years, 1:nyears[i])
    if( SmaxNorm==T){
      MuSmax[i] <- mean(EscDat$Escape)
      TauSmax[i] <- 1/(log(CVmean^2 + 1))
    }
  }    
  
  # convert obs spawners and recruits to matrices
  Smat <- matrix(NA, nrow=nCUs, ncol=max(nyears))
  lnRSmat <- matrix(NA, nrow=nCUs, ncol=max(nyears))
  for(k in 1:length(inds)){
    Smat[sampledCUs[k], Years[k]] <- Sdat[k]
    lnRSmat[sampledCUs[k], Years[k]] <- lnRSdat[k]
  }
  
  
  # create data input for JAGS
  Dlist <- list ("Smat"=Smat, "Rmat"=Rmat, "nyears"=nyears, "MuA"=1, "TauA"=0.001, "low"=1,
                 "upp"=upp, "nCUs"=nCUs)
  # Dlist <- list ("Smat"=Smat, "lnRSmat"=lnRSmat, "nyears"=nyears, "upp"=upp, "nCUs"=nCUs)
  simoutJAGS <- jags(Dlist, model.file = "C:/github/salmon-sim/scripts/jagsModels/MultiStockHierSimple.txt",
                     parameters.to.save=c("alpha", "beta", "MuA_global", "TauA_global", "tau_R"), 
                     n.chains = 3, n.iter = Niter, n.thin=Nthin, n.burnin=burnin)
  
  # change to bugs output to easily extract values
  simoutBUGS <- simoutJAGS$BUGSoutput
  
  # extract mcmc list to plot usuing coda functions
  if(plotDiags==T){
    ToPlot <- as.mcmc(simoutJAGS)
    # plot traces and posterior densities
    pdf( paste(here(),"/outputs/srAnalysis/figs/",fname,"TraceDens.pdf", sep="") )
    plot(ToPlot)
    dev.off()
    # plot gelman plots
    pdf( paste(here(),"/outputs/srAnalysis/figs/",fname,"Gelman.pdf", sep="" ))
    tryCatch( {
      gelman.plot(ToPlot)
    }, error=function(x=1){plot.new() })
    dev.off()
  }
  
  # Re-arrange and store each iteration's parameter estimates
  if(Niter>1000){
    simSamples <- sample_n(as.data.frame(simoutBUGS$sims.matrix), size=1000) #subset to 1000 samples equivalent to Fraser outputs generated by FRSSI modelers
  } else{
    simSamples <- simoutBUGS$sims.matrix
  }
  simPars <- NULL
  for(i in seq_along(CUs)){
    CU <- as.character(rep(CUs[i], each=nrow(simSamples)))
    stk <- i
    muAlpha <- simSamples[,"MuA_global"]
    tauAlpha <- simSamples[,"TauA_global"]
    alphas <- simSamples[,paste('alpha', "[", i , "]", sep="")]
    betas <- simSamples[,paste('beta', "[", i , "]", sep="")]
    tauRs <- simSamples[,paste('tau_R', "[", i , "]", sep="")]
    allParsTemp <- cbind(stk, CU, muAlpha, tauAlpha, alphas, betas, tauRs)
    simPars <- rbind(simPars, allParsTemp)
  }
  
  
  #save output for later
  if(saveMod==T){
    saveRDS(simoutJAGS, file=paste(here(),"/outputs/srAnalysis/data/",fname,"SimOutput",".rds", sep=""))
    write.csv(simPars, file=paste(here(),"/outputs/srAnalysis/data/",fname,"MCMCPars.csv", sep=""), row.names=FALSE)
  }
  
  # 11/30/15 save Gelman - Rubin Statistics, Geweke; need readable data frame with mod, year etc
  RhatVec <- simoutBUGS$summary[,"Rhat"]
  GewekeList <- geweke.diag(as.mcmc(simoutJAGS))
  GewekeTab <- data.frame(Param = names(GewekeList[[1]]$z) ,G1=GewekeList[[1]]$z, G2=GewekeList[[2]]$z, G3=GewekeList[[3]]$z )
  RhatTab <- data.frame(Mod=rep("RBHier", length(RhatVec)), Year = rep(max(dat$Year), length(RhatVec)), Param = names(RhatVec), 
                        Rhat=RhatVec, row.names=NULL)
  Rhats <- inner_join(RhatTab, GewekeTab, by="Param")
  CU <- NULL
  for(i in 1:dim(Rhats)[1]){
    Separated <- unlist(strsplit(as.character(Rhats$Param[i]), split="\\[|\\]"))
    # if there is a second item it will be the number, referring to the CU
    if(length(Separated)==2){
      CU[i] <- as.character(CUs[as.numeric(Separated[2])])
    } else { CU[i] <- unlist(strsplit(gsub("([0-9]+)","~\\1", fname), "~"))[1] }
  }
  Rhats$CU <- as.character(CU)
  
  # stopCluster(cl) #end cluster
}
#___________________________________________________________________________________________________________

#Function to sample time series with vector while accounting for autocorrelation
## NOT FUNCTIONAL ##
# bootMed <- function(tsVec, nSamples){
#   # tsb <- as.ts(tsVec)
#   
#   arFit <- auto.arima(tsb, max.p = 25, max.d = 0, max.q = 0, max.P = 0, max.Q = 0, max.D = 0,
#                       ic = "aic", max.order = 25, seasonal = FALSE)
#   modList <- list(order = c(arFit$arma[1], 0, 0), ar = coef(arFit))
#   modResid <- resid(arFit) - mean(resid(arFit))
#   
#   simTS <- function(res, n.sim, ran.args){
#     rg1 <- function(n, res){
#       sample(res, n, replace=TRUE)
#     }
#     ts.orig <- ran.args$ts
#     ts.mod <- ran.args$model
#     mean(ts.orig) + ts(arima.sim(model = ts.mod, n = n.sim, randGen = rg1, res = as.vector(res)))
#   }
#   
#   temp <- tsboot(modResid, median, R = nSamples, sim = "model",
#                  n.sim = length(tsb), orig.t = FALSE, ran.gen = simTS, 
#                  ran.args = list(ts = tsb, model=modList))
# }








### MUST BE EDITED IF USED ###



####################
##  FitRickerB    ##
####################
# fname="RickerB"; plotDiags=T; writeDat=T; saveMod=T; FoldName; Rhat=F;
# burnin = 25000; Niter=75000; Nthin=50; Fyear=F; threshSMSY=0.8; SmaxNorm=F; CVcap = 0.3; CVmean = 0.5;
# mult=1.5; savePosts=F
# Fit basic ricker model in Bayesian context
# 10/8/2015 changed to accomodated different number of years for diff sites
# Fyear=T means will only return "final year" of data; rather than joining with all data
# 10/20/2015 Now also calculates values and posteriors for BMs

FitRickerB<-function(dat, fname="RickerB", plotDiags=T, writeDat=T, saveMod=T, FoldName, Rhat=F,
                     burnin = 25000, Niter=75000, Nthin=50, Fyear=F, threshSMSY=0.8, SmaxNorm=F, CVcap = 0.3, CVmean = 0.5,
                     mult=1.5){
  
  dir.create(paste(here(),"/outputs/srAnalysis/figs/", FoldName, sep=""))
  dir.create(paste(here(),"/outputs/srAnalysis/figs/", FoldName, "/RickerBDiags", sep=""))
  dir.create(paste(here(),"/outputs/srAnalysis/dataOut/", FoldName, "Output", sep=""))
  
  # prepare data to feed into model
  # number of sites
  sites <- as.character(unique(dat$CU))
  nsites <- length(sites)
  
  # Set initial values 
  inits <- function(){
    list(tau_R=rep(1,nsites), tauv=rep(1,nsites), alpha=rep(0.5,nsites) )
  }
  
  #indices of non NA data points
  inds <- which(is.na(dat$Recruit)==F & is.na(dat$Escape)==F) 
  dataGood <- dat[inds,] 
  Rdat <- dataGood$Recruit
  Sdat <- dataGood$Escape
  
  # Years and sites corresponding to data points, as indices
  #Years <- data$Year[inds] - min(data$Year[inds]) +1
  # changed 10/7/2015 so that will work with data gaps
  Sites <- NULL
  for(i in 1:dim(dataGood)[1]){
    Sites[i] <- which(sites==dataGood$CU[i])
    if(is.na(Rdat[i])) {Sdat[i]<-NA}
  }
  # need upper bound for Smax prior for each site, number years for each site, Year indices for obs
  upp<-NULL
  nyears <- NULL
  Years <- NULL
  MuSmax <- NULL
  TauSmax <- NULL
  for(i in 1:nsites){
    #use data frame from ricker fit to use later to initialize model
    SiteDat <- dataGood[which(dataGood$CU==sites[i]),]
    #calculate upper bound on Smax prior *changed multiplier from 1.5 to 2
    upp[i] <- max(SiteDat$Escape, na.rm=TRUE)*mult
    nyears[i] <- dim(SiteDat)[1]
    Years <- c(Years, 1:nyears[i])
    #MuSmax will be capacity, or if NA, mean EScapement
    if( SmaxNorm==T){
      if(is.na(SiteDat$Capacity[1])){
        MuSmax[i] <- mean(SiteDat$Escape)
        TauSmax[i] <- 1/(log(CVmean^2 + 1))
      } else {
        #will all be the same so doesn't matter which use
        MuSmax[i] <- SiteDat$Capacity[1]
        TauSmax[i] <- 1/(log(CVcap^2 + 1))
      }
    }
  }    
  
  
  # create data input for JAGS
  if(SmaxNorm==F){
    Dlist <- list ("Sdat"=Sdat, "Rdat"=Rdat, "nyears"=nyears, 
                   "low"=1, "upp"=upp, "years"=Years, "sites"=Sites, "nsites"=nsites, nobs=length(Rdat))
    simoutJAGS <- jags(Dlist, model.file = paste(here(),"/scripts/jagsModels/SimpleRicker.txt", sep=""),
                       parameters.to.save=c("alpha", "beta", "tauv", "tau_R"), n.chains = 3, 
                       n.iter = Niter, n.thin=Nthin, n.burnin=burnin, inits=inits)
  } else {
    Dlist <- list ("Sdat"=Sdat, "Rdat"=Rdat, "nyears"=nyears, "MuSmax"=MuSmax, "TauSmax"=TauSmax,
                   "years"=Years, "sites"=Sites, "nsites"=nsites, nobs=length(Rdat))
    simoutJAGS <- jags(Dlist, model.file = paste(here(),"/scripts/jagsModels/SimpleRicker.txt", sep=""),
                       parameters.to.save=c("alpha", "beta", "tauv", "tau_R"), n.chains = 3, 
                       n.iter = Niter, n.thin=Nthin, n.burnin=burnin, inits=inits)
  }
  
  
  # change to bugs output to easily extract values
  simoutBUGS <- simoutJAGS$BUGSoutput
  
  # extract mcmc list to plot usuing coda functions
  
  if(plotDiags==T){
    ToPlot <- as.mcmc(simoutJAGS)
    # plot traces and posterior densities
    pdf( paste(here(),"/outputs/srAnalysis/figs/", FoldName, "/RickerBDiags/",fname,"TraceDens.pdf", sep="") )
    plot(ToPlot)
    dev.off()
    # plot gelman plots, Causes errors
    
    pdf( paste(here(),"/outputs/srAnalysis/figs/", FoldName, "/RickerBDiags/" ,fname,"Gelman.pdf", sep="" ))
    tryCatch( {
      gelman.plot(ToPlot)
    }, error=function(x=1){plot.new() })
    dev.off()
    
  }
  
  #extract posterior $ns and concatenate to vectors
  # all vectors will be length 7
  A <- simoutBUGS$median$alpha
  B <- simoutBUGS$median$beta
  tauv <- simoutBUGS$median$tauv
  taur <- simoutBUGS$median$tau_R
  # If only one site the indices on CI's are all messed up, need to add if
  if(nsites > 1){
    alphaCIs <- sapply(paste('alpha', "[", 1:nsites , "]", sep=""), function(x) quantile(simoutBUGS$sims.matrix[,x], c(.025,.975)))
    betaCIs <- sapply(paste('beta', "[", 1:nsites , "]", sep=""), function(x) quantile(simoutBUGS$sims.matrix[,x], c(.025,.975)))
    AUp <- alphaCIs[2,]
    ALow <- alphaCIs[1,]
    BUp <- betaCIs[2,]
    BLow <- betaCIs[1,]
  } else {alphaCIs <- quantile(simoutBUGS$sims.matrix[,'alpha'], c(.025,.975))
  betaCIs <- quantile(simoutBUGS$sims.matrix[,'beta'], c(.025,.975))
  AUp <- alphaCIs[2]
  ALow <- alphaCIs[1]
  BUp <- betaCIs[2]
  BLow <- betaCIs[1]
  }
  # Re-arrange and store each iteration's parameter estimates
  if(Niter>1000){
    simSamples <- sample_n(as.data.frame(simoutBUGS$sims.matrix), size=1000) #subset to 1000 samples equivalent to Fraser
  } else{
    simSamples <- simoutBUGS$sims.matrix
  }
  simPars <- NULL
  for(i in seq_along(sites)){
    stk <- rep(sites[i], each=nrow(simSamples))
    alphas <- simSamples[,paste('alpha', "[", i , "]", sep="")]
    betas <- simSamples[,paste('beta', "[", i , "]", sep="")]
    tauVs <- simSamples[,paste('tauv', "[", i , "]", sep="")]
    tauRs <- simSamples[,paste('tau_R', "[", i , "]", sep="")]
    allParsTemp <- cbind(stk,alphas, betas, tauVs, tauRs)
    simPars <- rbind(simPars, allParsTemp)
  }
  
  #save output for later
  if(saveMod==T){
    saveRDS(simoutJAGS, file=paste(here(),"/outputs/srAnalysis/dataOut/",fname,"simOutput",".rds", sep=""))
    write.csv(simPars, file=paste(here(),"/outputs/srAnalysis/dataOut/",fname,"MCMCPars.csv", sep=""), row.names=FALSE)
  }
  
  # save Gelman - Rubin Statistics, Geweke statistics; need readable data frame with mod, year etc
  RhatVec <- simoutBUGS$summary[,"Rhat"]
  GewekeList <- geweke.diag(as.mcmc(simoutJAGS))
  GewekeTab <- data.frame(Param = names(GewekeList[[1]]$z) ,G1=GewekeList[[1]]$z, G2=GewekeList[[2]]$z, G3=GewekeList[[3]]$z )
  RhatTab <- data.frame(Mod=rep("RickerB", length(RhatVec)), Year = rep(max(dat$Year), length(RhatVec)), Param = names(RhatVec), 
                        Rhat=RhatVec, row.names=NULL)
  Rhats <- inner_join(RhatTab, GewekeTab, by="Param")
  
  # add column with site name
  CU <- NULL
  for(i in 1:dim(Rhats)[1]){
    Separated <- unlist(strsplit(as.character(Rhats$Param[i]), split="\\[|\\]"))
    # if there is a second item it will be the number, referring to the site
    if(length(Separated)==2){
      CU[i] <- as.character(sites[as.numeric(Separated[2])])
    } else { CU[i] <- NA }
  }
  Rhats$CU <- as.character(CU)
  
  #DOn't need all years fed back, if don't want Fyear=T
  if(Fyear==F){
    RBdat <- data.frame("CU"=sites[Sites], "Year"=dataGood$Year, "RickerA"=rep(exp(A),nyears), "RickerB"= rep(B, nyears), 
                        "TauV"= rep(tauv, nyears), "TauR"= rep(taur, nyears), AUp =rep(exp(AUp), nyears), 
                        ALow=rep(exp(ALow), nyears), BUp=rep(BUp,nyears), BLow=rep(BLow, nyears))
    # , SMSYUp=rep(SMSYUp, nyears), SGen=rep(SGen, nyears),
    # LowMin = rep(LowMin, nyears), LowMax = rep(LowMax, nyears), UppMin = rep(UppMin, nyears), UppMax = rep(UppMax, nyears))
    # 
    ## RBdat does not include all years of obs, need to merge will full data set
    RBAllDat <- full_join(dat, RBdat, by=c("CU", "Year"))
    if(writeDat==T){
      write.csv(RBAllDat, paste(here(),"/outputs/srAnalysis/dataOut/",fname,"fitRickerB.csv", sep=""))
    }
    if(Rhat == F){ 
      out <- RBAllDat 
    } else {
      out <- list(RBAllDat, Rhats)
    }
    
  } else {
    #only return final year of data and parameter estimates
    AllDat <- data.frame("CU"=sites, "Year"=rep(max(dataGood$Year), nsites), "RickerA"=exp(A), "RickerB"= B, 
                         "TauV"= tauv, "TauR"= taur, AUp =exp(AUp), ALow=exp(ALow), BUp=BUp, BLow=BLow, Mod=rep("RickerB", nsites))
    # ,
    #                     "SMSYUp"=SMSYUp, "SGen"=SGen, "LowMin"=LowMin, "LowMax"=LowMax, "UppMin"=UppMin, "UppMax"=UppMax, row.names=NULL)
    # #note there might not actually be new data in that year, will have to control for that in data fed in
    if(Rhat == F){ 
      out <- AllDat 
    } else {
      out <- list(AllDat, Rhats) #removed posts
    }
  }
  
  out
  
}
#___________________________________________________________________________________________________________
