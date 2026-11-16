
## A Solution to Assignment

## Libraries:
library(ggplot2)
library(dplyr)
library(RTMB)

## Load data:
## Path to folder...
sox <- read.csv("../data/sockeye/SRDATA2026 1.csv")
stock <- read.csv("../data/sockeye/FC_StockInfo 2.csv")
covariates <- read.csv("../data/sockeye/FC_Environmental_Data_2026_Std_Offset_Real.csv")
sox <- sox %>% left_join(stock, by = "PopID")
sox <- sox %>% left_join(covariates, by = "yr")
sox <- sox %>% filter(!is.na(spawners), !is.na(recruits))
sox <- sox %>% group_by(PopName) %>% 
  mutate(sd_spawners = sd(spawners), mean_spawners = mean(spawners),
    spawners_scaled = (spawners)/sd_spawners) %>% 
  ungroup()

sox <- sox %>% 
    filter(!is.na(spawners) & !is.na(recruits))

sox <- sox %>% mutate(logRS = log(recruits/spawners))
sox <- sox %>% 
  mutate(stock = as.numeric(factor(PopName)))

nstocks <- max(sox$stock)

X_a <- model.matrix(~ 0 + scale(jlpsst) + scale(pdo) + scale(peak), data = sox)
X_b <- model.matrix(~ 0 + PopName:spawners_scaled, data = sox)

params <- list(alpha0 = numeric(nstocks),
  alpha = numeric(ncol(X_a)), 
  logbeta = numeric(ncol(X_b)), 
  logsigma = numeric(nstocks))

negll <- function(pars){
  getAll(pars, sox, warn = TRUE)
  sigma <- exp(logsigma)
  beta <- exp(logbeta)
  ElogRS <- alpha0[stock] + X_a %*% alpha - X_b %*% beta
  nll <- -sum(dnorm(logRS, ElogRS, sigma[stock], log = TRUE))
  return(nll)
}

obj <- MakeADFun(negll, params, silent = TRUE)
fit <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1000, eval.max = 2000))
sdrep <- sdreport(obj, par.fixed = fit$par)
sumfix <- summary(sdrep, "fixed")

## Check out these alpha effects

