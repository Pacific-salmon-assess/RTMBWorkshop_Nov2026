## Multistock and age Ricker Model
## Based on RTMB workshop

library(ggplot2)
library(dplyr)
library(RTMB)
library(tidyr)

## Read in the data:
sox <- read.csv("../data/sockeye/SRDATA2026 1.csv")
stock <- read.csv("../data/sockeye/FC_StockInfo 2.csv")
covariates <- read.csv("../data/sockeye/FC_Environmental_Data_2026_Std_Offset_Real.csv")
sox <- sox %>% left_join(stock, by = "PopID")
sox <- sox %>% left_join(covariates, by = "yr")
sox <- sox %>% filter(!is.na(spawners), !is.na(recruits))

## Remove some problematic stocks that we identified as having different, or poorly estimated smax.
sox <- sox %>%  filter(!PopName %in% c("Cultus", "Nadina", "Harrison"))

## Scale the spawners, but didn't centre them. May not be necessary anymore now that we removed the more extreme cases.
sox <- sox %>% group_by(PopName) %>% 
  mutate(sd_spawners = sd(spawners), mean_spawners = mean(spawners),
    spawners_scaled = (spawners)/sd_spawners) %>% 
  ungroup()

## Log Recruits, including ages
sox <- sox %>% mutate(logRS = log(recruits/spawners)) %>%
  mutate(logRS4 = log((age4 + 1e-6)/spawners), logRS5 = log((age5 + 1e-6)/spawners))
## Make a stock index and remove NA spawners ( we will deal with NA recruits in model )
sox <- sox %>% 
  mutate(stock = as.numeric(factor(PopName))) %>%
  filter(!is.na(spawners))

# sox <- sox %>% filter(!is.na(jlpsst), !is.na(pdo), !is.na(peak))
sox <- sox %>% mutate(time = yr - min(yr) + 1)
nyears <- max(sox$time)
minyear <- min(sox$yr)
nstocks <- max(sox$stock)

## Correlation within year between ages, varying by stock, with stock random effect. (Full model)
negll <- function(pars, data){
  
  getAll(pars, data, warn = FALSE)
  beta <- exp(logbeta)
  
  sigma_4 <- exp(logsigma_4)
  sigma_5 <- exp(logsigma_5)
  sigma_re <- exp(logsigma_re)
  sigma_time <- exp(logsigma_time)
  
  phi_time <- 2*plogis(logitphi_time)-1
  phi_age <- 2*plogis(logitphi_age)-1

  REPORT(phi_time)
  REPORT(phi_age)
  
  ## Age 4 is base productivity
  mu_age4 <- alpha0 + re_4[stock] + X %*% a + re_time[time] - beta[stock]*spawners
  ## Age 5 is shared but with a shift per stock
  mu_age5 <- mu_age4 + alpha5 + re_5[stock]
  mu_age5_cond <- mu_age5 + phi_age[stock]*sigma_5[stock]/sigma_4[stock]*(logRS4 - mu_age4)

  logRS4 <- OBS(logRS4)
  logRS5 <- OBS(logRS5)

  nyears <- max(time)
  ll <- dnorm(re_time[1], 0, sqrt(sigma_time^2/(1-phi_time^2)), log = TRUE)
  for( t in 2:nyears ) 
    ll <- ll + dnorm(re_time[t], phi_time*re_time[t-1], sigma_time, log = TRUE)

  ll <- ll + sum(dnorm(re_4, 0, sigma_re[1], log = TRUE))
  ll <- ll + sum(dnorm(re_5, 0, sigma_re[2], log = TRUE))

  mis4 <- is.na(logRS4)
  mis5 <- is.na(logRS5) | is.na(logRS4)
  mis5_only <- is.na(logRS5) & !is.na(logRS4)

  ll <- ll + sum(dnorm(logRS4[!mis4], mu_age4[!mis4], sigma_4[stock[!mis4]], log = TRUE))
  ll <- ll + sum(dnorm(logRS5[!mis5], mu_age5_cond[!mis5], sigma_5[stock[!mis5]], log = TRUE))
  
  ## Let's predict the ages and report them:
  ## Use the ad_context trick to do some data processing etc within the object.
  ## Won't help me to do prediction intervals. Will have to add to it for that
  if (!RTMB:::ad_context()) {
    predict <- NULL
    if(any(mis4)){
      predict <- rbind(predict, data.frame(year = yr[mis4], stock = stock[mis4], PopName = PopName[mis4], spawners = spawners[mis4],
        Expected = exp(mu_age4[mis4] + log(spawners[mis4])), age = 4, sigma = sigma_4[stock[mis4]]))
    }
    if(any(mis5_only)){
      predict <- rbind(predict, data.frame(year = yr[mis5_only], stock = stock[mis5_only], PopName = PopName[mis5_only], spawners = spawners[mis5_only],
        Expected = exp(mu_age5_cond[mis5_only] + log(spawners[mis5_only])), age = 5, sigma = sigma_5[stock[mis5_only]]))
    }
    REPORT(predict)
  }
  pred_mu4 <- mu_age4[mis4]
  pred_mu5 <- mu_age5_cond[mis5_only]
  ADREPORT(pred_mu4)
  ADREPORT(pred_mu5)

  return(-ll)
}

X <- model.matrix(~ 0 + scale(jnpsst) + scale(maesst) + scale(stn2ja) + scale(stn2js) + scale(peak) + scale(pdo) + scale(aflow) + scale(mflow) + scale(jflow), data = sox)
parameters <- list(
  alpha0 = runif(1, 2, 5),
  alpha5 = -1,
  a = rnorm(ncol(X), 0, 0.01),
  re_4 = rnorm(nstocks),
  re_5 = rnorm(nstocks),  
  logbeta = runif(nstocks, 0, 1),
  logsigma_4 = runif(nstocks, 3, 5),
  logsigma_5 = runif(nstocks, 3, 5),  
  logsigma_re = runif(2, 0.5, 2),
  re_time = rnorm(max(sox$time)),
  logitphi_time = 0,
  logitphi_age = numeric(nstocks),  
  logsigma_time = 0
)
join_fun <- function(f, data) function(pars) f(pars, data)
obj <- MakeADFun(func = join_fun(negll, sox), parameters = parameters, random = c("re_4", "re_5", "re_time"), silent = TRUE)
fit <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 2000, iter.max = 1000))
fit$convergence

## Some stocks have strong age relationships and some do not.
obj$report()$phi_age

par_fit <- obj$env$parList(fit$par)

report_values <- obj$report()
predicted <- report_values$predict

## Function to do prediction intervals: Basic Idea: SE of estimate + observation level noise.
## Will need to ADREPORT the means.
forecast <- function(obj){
  repvals <- obj$report()
  sdr <- sdreport(obj)
  sdr_rep <- summary(sdr, "report")
  pred <- repvals$predict
  mu4 <- as.data.frame(sdr_rep[rownames(sdr_rep) == "pred_mu4",])
  mu5 <- as.data.frame(sdr_rep[rownames(sdr_rep) == "pred_mu5",])
  pred <- cbind(pred, rbind(mu4, mu5)) 
  pred$se_pred <- sqrt(pred$`Std. Error`^2 + pred$sigma^2)  ## This bit makes it painfully wide prediction interval.
  
  ## Will use simulation to deal with transformation bias and also for se to combine levels while 
  ## properly accounting for confidence levels
  simmat <- replicate(10000, exp(rnorm(nrow(pred), mean = pred$Estimate, sd = pred$se_pred))*pred$spawners)
  out <- NULL
  for( i in 1:nrow(pred) ){
    simest <- simmat[i,]
    out <- rbind(out, data.frame(lowpred = quantile(simest, 0.025), uppred = quantile(simest, 0.975),
                        meanpred = mean(simest), medianpred = median(simest), sepred = sd(simest)))
  }
  pred <- cbind(pred, out)
  
  stocks <- unique(pred$stock)
  pred_recruits <- NULL
  for(i in seq_along(stocks)){
    idx <- which(pred$stock == stocks[i])
    simest <- colSums(simmat[idx, ])
    pred_recruits <- rbind(pred_recruits, data.frame(stock = stocks[i], lowpred = quantile(simest, 0.025), uppred = quantile(simest, 0.975),
                        meanpred = mean(simest), medianpred = median(simest), sepred = sd(simest)))
  }
  return(list(pred = pred, pred_recruits = pred_recruits))
}

## Will try and loop through: (Brood years)
years <- 2000:2021
pred <- NULL
pred_recruits <- NULL
for( i in seq_along(years) ){
  dat <- sox %>% filter(yr <= years[i])
  dat$logRS4[dat$yr == years[i]] <- NA
  dat$logRS5[dat$yr == years[i]] <- NA
  dat$logRS5[dat$yr == years[i]-1] <- NA
  ## Note that the final year is missing, but we will infer the last AR1 term to be used in prediction.

  ## Subset some inputs for the 
  nyears <- max(dat$time)
  parameters$re_time <- rnorm(nyears)
  
  ## This is still passed via environment
  X <- model.matrix(~ 0 + scale(jnpsst) + scale(maesst) + scale(stn2ja) + scale(stn2js) + scale(peak) + scale(pdo) + scale(aflow) + scale(mflow) + scale(jflow), data = dat)
  parameters$a <- numeric(ncol(X))  ## Just in case I change X here.

  obj <- MakeADFun(join_fun(negll, dat), parameters, random = c("re_4", "re_5", "re_time"), silent = TRUE)
  fit <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 2000, iter.max = 1000))
  predi <- forecast(obj)
  predi$pred$predict_year <- years[i] + 4
  predi$pred_recruits$predict_year <- years[i] + 4
  pred <- rbind(pred, predi$pred)
  pred_recruits <- rbind(pred_recruits, predi$pred_recruits)
}
pred_recruits <- pred_recruits %>% left_join(sox %>% distinct(stock, PopName))

recruits <- sox %>% 
  pivot_longer(c(age4, age5)) %>% 
  mutate(year = ifelse(name == "age4", yr + 4, yr + 5)) %>%
  group_by(year, PopName) %>% summarize(recruits = sum(value), .groups = "drop")

# pred_r <- pred %>%  
  # group_by(predict_year, PopName) %>% summarize(recruits = sum(Expected), .groups = "drop")

ggplot(recruits, aes(x = year, y = recruits)) + 
  geom_line() + theme_bw() + xlab("Return Year") + ylab("Return (x millions)") + 
  geom_line(data = pred_recruits, aes(x = predict_year, y = medianpred), colour = "red", linetype = 1) +
  # geom_line(data = pred_recruits, aes(x = predict_year, y = meanpred), colour = "blue", linetype = 2) +  
  facet_wrap(~ PopName, scale = "free_y")

## Plot with prediction uncertainty... Yikes! Could also plot with just expected uncertainty, without additional sigma...
ggplot(recruits, aes(x = year, y = recruits)) + 
  geom_line() + theme_bw() + xlab("Return Year") + ylab("Return (x millions)") + 
  geom_line(data = pred_recruits, aes(x = predict_year, y = medianpred), colour = "red") +
  geom_ribbon(data = pred_recruits, aes(x = predict_year, y = medianpred, ymin = lowpred, ymax = uppred), fill = 'red', alpha = 0.1) +
  facet_wrap(~ PopName, scale = "free_y")