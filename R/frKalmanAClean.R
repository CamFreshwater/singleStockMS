#*******************************************************************************
# frKalmanAClean.R
# Date last revised: 29 July, 2019
# To parameterize alterantive static productivity scenarios examine Kalman alpha
# values estimated by C. Holt for subset of Fraser River stocks. For each
# stock calculate difference between max/min alpha and most recent estimate 
# (from frStockRecModels.R).
#*******************************************************************************

library(here); library(tidyverse)

cuPars <- read.csv(here("data", "fraserCUPars.csv"))

# Function from C. Holt to extract first and last with additions for max/min
extractRicaFirstLast <- function(){
  data<-as.tbl(read.csv(here("data", "old",
                             "RickerPars26Oct2017_forBenchmarks.csv"), 
                        col.name=c("Label", "Stock", "StockID", "Model",
                                   "Parameter", "Value",	"SD",	"LCL",
                                   "LQ", "Median",	"UQ",	"UCL"), skip=2))
  
  addChars <- data %>% 
    mutate(chars=as.character(Parameter))
  isTvA <- startsWith(addChars$chars,"alpha[")
  data <- data%>%add_column(isTvA)  # Add variable to distinguish tv alphas from other parameters
  
  x <- data %>% 
    filter(Model=="RickerTVA", isTvA=="TRUE")
  firstA <- x %>% 
    filter(Parameter=="alpha[1]") %>% 
    select(Stock, Median) %>% 
    rename(firstA = Median)
  otherAs <- x %>% 
    group_by(Stock) %>% 
    summarise(lastA = last(Median),
              maxA = max(Median),
              minA = min(Median))
  return(left_join(firstA, otherAs, by="Stock"))
}#End of function

kalAEstimates <- extractRicaFirstLast() %>% 
  mutate(stkName = recode_factor(Stock, "Pitt" = "Upper Pitt")) %>% 
  select(-Stock)

## Merge subset stocks that have kalman a estimates (all Ricker)
alphaEst <- Reduce(function(x, y) merge(x, y, by=c("stkName")), 
                     list(kalAEstimates, cuPars[ , c("stkName", "alpha")])) %>% 
  mutate(maxInc = (maxA - alpha) / alpha,
         maxDec = (alpha - minA) / alpha,
         curChange = (lastA - alpha) / alpha)

hist(alphaEst$maxInc)
hist(alphaEst$maxDec)
mean(alphaEst$maxInc) 
mean(alphaEst$maxDec)
median(alphaEst$maxInc)
median(alphaEst$maxDec)
#use 35% for both to keep balanced
range(alphaEst$curChange)
hist(alphaEst$curChange)
mean(alphaEst$curChange)
