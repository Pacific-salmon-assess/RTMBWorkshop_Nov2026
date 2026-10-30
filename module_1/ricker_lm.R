
## Solution to Assignment

## Libraries:
library(ggplot2)
library(dplyr)

## Load data:
## Path to folder...
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
sox <- read.csv("../data/sockeye/SRDATA2026 1.csv")
stock <- read.csv("../data/sockeye/FC_StockInfo 2.csv")
sox <- sox %>% left_join(stock, by = "PopID")

## Option 1:
## Nested function for local environment
my_lm <- function(formula, data){
  ## Design Matrix:
  X <- model.matrix(formula, data)
  ## Response Variable:
  fr <- model.frame(formula, data = data)
  y <- model.response(fr)

  npars <- ncol(X)
  pars_init <- numeric(npars+1)
  pars_init[npars+1] <- log(sd(y))
  n <- nrow(X)
  
  ## Nested negative log-likelihood
  ## Parent environment accessible to function this way.
  nll <- function(pars){
    beta <- pars[1:npars]
    sigma <- exp(pars[npars+1])
    mu <- X %*% beta
    -sum(dnorm(y, mu, sigma, log = TRUE))
  }

  fit <- nlminb(pars_init, nll, control = list(rel.tol = 1e-12, abs.tol = 1e-12))
  par_values <- fit$par
  beta_fit <- par_values[1:npars]
  if(fit$convergence != 0) cat(" [Warning]", fit$message, ".\n")
  names(beta_fit) <- colnames(X)
    return(list(
      beta = beta_fit, 
      sigma = sqrt(n/(n-npars))*exp(par_values[npars+1]), 
      ll = -fit$obj))
}

fit1 <- my_lm(formula = log(recruits/spawners) ~ 0 + PopName + PopName:spawners, data = sox)
fit <- lm(log(recruits/spawners) ~ 0 + PopName + PopName:spawners, data = sox)
range(coef(fit) - fit1$beta) ## Not ideal...
sigma(fit) - fit1$sigma  ## A little bit different...
logLik(fit) - fit1$ll  ## Pretty much same log-likelihood, must be a very flat bumpy surface...

## Let's scale things and see if it fits easier.
sox <- sox %>% 
  mutate(logRS = log(recruits/spawners)) %>% 
  group_by(PopName) %>% 
  mutate(spawners_scaled = scale(spawners))
fit1 <- my_lm(formula = logRS ~ 0 + PopName + PopName:spawners_scaled, data = sox)
fit <- lm(logRS ~ 0 + PopName + PopName:spawners_scaled, data = sox)
range(coef(fit) - fit1$beta) 
## Much better! Numerical issues... Some beta values can be very small. Better to scale them


## Non-nested functions and efficiently passing data
## without creating copies.
nll <- function(pars){
  beta <- pars[1:npars]
  sigma <- exp(pars[npars+1])
  mu <- X %*% beta
  -sum(dnorm(y, mu, sigma, log = TRUE))
}

my_lm2 <- function(formula, data){
  ## Design Matrix:
  X <- model.matrix(formula, data)
  ## Response Variable:
  fr <- model.frame(formula, data = data)
  y <- model.response(fr)

  npars <- ncol(X)
  pars_init <- numeric(npars+1)
  
  ## *** This function will now inherit access to everything
  ## contained in the my_lm function environment.
  environment(nll) <- environment()   

  fit <- nlminb(pars_init, nll)
  par_values <- fit$par
  beta_fit <- par_values[1:npars]
  names(beta_fit) <- colnames(X)
    return(list(beta = beta_fit, sigma = exp(par_values[npars+1])))
}

fit2 <- my_lm2(formula = log(recruits/spawners) ~ 0 + PopName + PopName:spawners, data = sox)