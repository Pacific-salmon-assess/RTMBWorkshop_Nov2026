
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
