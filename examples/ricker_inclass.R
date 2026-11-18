
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

library(splines2)
K <- 8
knots <- seq(5, 70, length = K-1)

bsMat <- bSpline(sox$yr-min(sox$yr), 
                 knots = knots, 
                 degree = 3,    ## Cubic
                 periodic = FALSE, ## Change later.
                 intercept = FALSE)
class(bsMat) <- "matrix"

bspred <- bSpline(0:74, 
                 knots = knots, 
                 degree = 3,    ## Cubic
                 periodic = FALSE, ## Change later.
                 intercept = FALSE)
class(bspred) <- "matrix"


X_a <- model.matrix(~ 0 + scale(jlpsst) + scale(pdo) + scale(peak), data = sox)
X_b <- model.matrix(~ 0 + PopName:spawners_scaled, data = sox)

params <- list(alpha0 = numeric(nstocks),
  alpha = numeric(ncol(X_a)), 
  b = numeric(ncol(bsMat)),
  logbeta = numeric(ncol(X_b)), 
  logsigma = numeric(nstocks))

negll <- function(pars){
  getAll(pars, sox, warn = TRUE)
  sigma <- exp(logsigma)
  beta <- exp(logbeta)
  ElogRS <- numeric(nrow(sox))
  for( i in seq_along(ElogRS)){
    ElogRS[i] <-  alpha0[stock[i]] + 
      sum(X_a[i,] * alpha) + 
      sum(bsMat[i,] * b) - sum(X_b[i,] * beta)
  }
  pred_wiggle <- bspred %*% b
  REPORT(pred_wiggle)
  REPORT(ElogRS)
  residuals <- (logRS-ElogRS)/sigma[stock]
  REPORT(residuals)
  nll <- -sum(dnorm(logRS, ElogRS, sigma[stock], log = TRUE))
  return(nll)
}

obj <- MakeADFun(negll, params, silent = TRUE, map = list(b = factor(rep(NA, ncol(bsMat)))))
fit <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1000, eval.max = 2000))
sdrep <- sdreport(obj, par.fixed = fit$par)
sumfix <- summary(sdrep, "fixed")

## Check out these alpha effects
est_alpha <- obj$env$parList(fit$par)$alpha
est_alpha_se <- obj$env$parList(sumfix[, "Std. Error"])$alpha
cbind(est_alpha - 1.96*est_alpha_se, est_alpha + 1.96*est_alpha_se)

sox$residuals <- obj$report(fit$par)$residuals
ggplot(data = sox, aes(x = yr, y = residuals)) + 
  geom_point(aes(colour = PopName)) + 
  theme_minimal()

wiggle <- obj$report(fit$par)$pred_wiggle
plot(1:75, wiggle, type ='l')


## Day 2: Temporal Autocorrelation
sox <- sox %>% mutate(time = yr - min(yr) + 1)
nyears <- max(sox$time)
minyear <- min(sox$yr)

negll <- function(pars, data){
  getAll(pars, data, warn = FALSE)
  logRS <- OBS(logRS)
  sigma <- exp(logsigma)
  sigma_re <- exp(logsigma_re)
  phi <- plogis(logitphi)*2-1
  sigma0 <- sigma_re/sqrt(1-phi^2)
  nobs <- nrow(data)
  beta <- exp(logbeta)
  REPORT(alpha_re)
  ElogRS <- alpha_re + alpha0 + beta*spawners_scaled
  REPORT(ElogRS)
  residuals <- (logRS-ElogRS)/sigma
  REPORT(residuals)
  ADREPORT(phi)
  
  nll <- -dnorm(alpha_re[1], 0, sigma0, log = TRUE)
  for(t in 2:nobs){
    nll <- nll - dnorm(alpha_re[t], alpha_re[t-1]*phi, sigma_re, log = TRUE)
  }
  nll <- nll - sum(dnorm(logRS, ElogRS, sigma, log = TRUE))
  return(nll)
}

data <- sox %>% filter(PopName == "Stellako") %>% arrange(time)
params <- list(alpha0 = 0,
  alpha_re = numeric(nrow(data)), 
  logbeta = 0, 
  logsigma = 0,
  logsigma_re = 0,
  logitphi = 0)

comb_fn <- function(f, data){function(pars){f(pars, data)}}
obj <- MakeADFun(comb_fn(negll, data), params, silent = TRUE, random = "alpha_re")
fit <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1000, eval.max = 2000))
plot(alpha_re, type = 'l')
sdrep <- sdreport(obj, par.fixed = fit$par)
sumfix <- summary(sdrep, "fixed")
sumrep <- summary(sdrep, "report")

check <- checkConsistency(obj, par = fit$par, estimate = TRUE)
check


## Multistock
## Path to folder...
library(ggplot2)
library(dplyr)
library(RTMB)

sox <- read.csv("../data/sockeye/SRDATA2026 1.csv")
stock <- read.csv("../data/sockeye/FC_StockInfo 2.csv")
covariates <- read.csv("../data/sockeye/FC_Environmental_Data_2026_Std_Offset_Real.csv")
sox <- sox %>% left_join(stock, by = "PopID")
sox <- sox %>% left_join(covariates, by = "yr")
sox <- sox %>% filter(!is.na(spawners), !is.na(recruits))
sox <- sox %>% group_by(PopName) %>% 
  mutate(sd_spawners = sd(spawners), mean_spawners = mean(spawners),
    spawners_scaled = (spawners)/sd_spawners) %>% 
  ungroup() %>% filter(!PopName %in% c("Cultus", "Nadina", "Harrison"))
sox <- sox %>% mutate(logRS = log(recruits/spawners))
sox <- sox %>% 
  mutate(stock = as.numeric(factor(PopName)))
sox <- sox %>% 
    filter(!is.na(spawners) & !is.na(recruits))

sox <- sox %>% filter(!is.na(jlpsst), !is.na(pdo), !is.na(peak))
sox <- sox %>% mutate(time = yr - min(yr) + 1)
nyears <- max(sox$time)
minyear <- min(sox$yr)
nstocks <- max(sox$stock)

negll_ar <- function(pars){
  
  getAll(pars, sox, warn = FALSE)
  beta <- exp(logbeta)
  
  sigma <- exp(logsigma)
  sigma_re <- exp(logsigma_re)
  sigma_time <- exp(logsigma_time)
  phi <- 2*plogis(logitphi)-1
  ADREPORT(phi)
  lambda <- exp(loglambda)
  
  alpha_stock <- re + alpha0
  ADREPORT(alpha_stock)
  
  mu <- alpha_stock[stock] + X %*% a + re_time[time] - beta[stock]*spawners_scaled
  logRS <- OBS(logRS)
  
  ll <- dnorm(re_time[1], 0, sqrt(sigma_time^2/(1-phi^2)), log = TRUE)
  for( t in 2:nyears ) 
    ll <- ll + dnorm(re_time[t], phi*re_time[t-1], sigma_time, log = TRUE)

  ll <- ll + sum(dnorm(re, 0, sigma_re, log = TRUE))
  ll <- ll + sum(dnorm(logRS, mu, sigma[stock], log = TRUE))
  # ll <- ll + dgamma(sigma_time, 5, 500, log = TRUE)
  ## Ridge Regression
  ll <- ll - lambda*sum(a^2)
  REPORT(mu)
  REPORT(re_time)
  
  return(-ll)
}

X <- model.matrix(~ 0 + scale(maesst) + scale(stn2ja) + scale(peak) + scale(aflow) + scale(pdo), data = sox)
parameters <- list(
  alpha0 = runif(1, 2, 5),
  a = rnorm(ncol(X), 0, 0.01),
  re = rnorm(nstocks),
  logbeta = runif(nstocks, 0, 1),
  logsigma = runif(nstocks, 3, 5),
  logsigma_re = runif(1, 0.5, 2),
  re_time = rnorm(nyears),
  logitphi = 0,
  logsigma_time = 0,
  loglambda = 0
)

obj_ar <- MakeADFun(func = negll_ar, parameters = parameters, 
  random = c("re", "re_time"), silent = TRUE)
fit_ar <- nlminb(obj_ar$par, obj_ar$fn, obj_ar$gr, control = list(eval.max = 1000, iter.max = 1000))
fit_ar

par_fit <- obj_ar$env$parList(fit_ar$par)
# which.min(par_fit$logbeta)

sdr <- sdreport(obj_ar, par = fit_ar$par)
sum_ad <- summary(sdr, "report")
sum_fe <- summary(sdr, "fixed")

est_fit <- obj_ar$env$parList(sum_fe[,1])
sd_fit <- obj_ar$env$parList(sum_fe[,2])

ci <- est_fit$logitphi + c(-1, 1) * sd_fit$logitphi*1.96
ci_phi <- plogis(ci)*2-1
ci_phi2 <- sum_ad[rownames(sum_ad) == "phi","Estimate"] + c(-1.96, 1.96) * sum_ad[rownames(sum_ad) == "phi","Std. Error"]
ci_phi
ci_phi2

prof <- TMB::tmbprofile(obj_ar, "logitphi", trace = FALSE)
plot(prof)
plogis(confint(prof))*2-1

# check <- checkConsistency(obj_ar, n = 100)

re_time <- obj_ar$report()$re_time
ggplot() + geom_line(aes(x = minyear+0:(nyears-1), y = re_time)) + 
  theme_minimal() + geom_hline(yintercept = 0, col = 'red', linetype = 'dashed') + 
  xlab("Year") + ylab("Productivity Time Effect")
