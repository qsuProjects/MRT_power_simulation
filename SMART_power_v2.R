

library(tidyverse)
library(lmerTest)
library(gtsummary)


# Setup for patient i from CDCES j
set.seed(13507)

data.frame(
  # One your of follow up
  study_week = 12:52) %>%
  mutate(
    # tide week of the month
    tide_week = (study_week - 1) %% 4 + 1,
    # default CPM week of the month (equally distributed)
    default_rmp_week = sample(1:4, 1),
    # RMP happening that week (Y/N), by default?
    default_rmp = as.numeric(tide_week == default_rmp_week),
    # eligible for randomization if:
    #   week following a default RPM week,
    #     ie, tide_week == default_rmp_week + 1 (if 5 -> 1)
    #.    ie, tide_week == default_rmp_week %% 4 + 1
    elig_rand_1 = as.numeric(tide_week == default_rmp_week %% 4 + 1)
  ) %>% as_tibble()




# TIR distribution for weekly reviewed pts: Beta(10,5)
a = 10
b_monthly = 3.8
a / (a + b_monthly)
pbeta(.65, a, b_monthly) # prob of TIR < 65%

# TIR distribution for monthly reviewed pts: Beta(10,2.5)
b_weekly = 2.5
a / (a + b_weekly)
pbeta(.65, a, b_weekly) # prob of TIR < 65%

p = seq(0,1, length=100)
plot(p, dbeta(p, 10, 2.5), ylab='density', type ='l', col='red')
lines(p, dbeta(p, 10, 5), col='blue')
legend('topleft', c('Beta(10, 2.5)','Beta(10, 5)'),
       lty=c(1,1),col=c('red', 'blue'))


# Generate pt's trajectory
generate_4t_pt <- function() {
  trt <- 'Monthly' # everyone starts month 3 on monthly cadence
  week <- 12 # microrandomization starts at 12 weeks (6 months)
  weekk_tir <- rbeta(1, a, b_monthly)
  rand_eligible <- 1*(last_wk_tir < .65) # 1st eligible if not doing well
  rand_week <- 0 # randomization performed at every 4-weeks
  
  dt <- data.frame(trt, week, last_wk_tir, rand_week, rand_eligible)
  
  for (i in 1:25) {
    low_tir <- (last_wk_tir < .65)
    high_tir <- (last_wk_tir > .75)
    
    week <- week + 1
    rand_week <- (week %% 4) == 0 # randomization happen if week mod 4 = 0
    
    if (trt == 'Monthly') {
      rand_eligible <-  rand_week * low_tir
      # last_wk_tir <- rbeta(1, a, b_monthly)
    } else {
      rand_eligible <-  rand_week * high_tir
      # last_wk_tir <- rbeta(1, a, b_weekly)
    }
    
    if (rand_eligible == 1) {
      # randomize
      trt <- sample(c('Monthly', 'Weekly'), 1)
      # generate next week's TIR
      last_wk_tir <- rbeta(1, a, ifelse(trt == 'Monthly', b_monthly, b_weekly))
    } else {
      # generate next week's TIR
      last_wk_tir <- rbeta(1, a, ifelse(trt == 'Monthly', b_monthly, b_weekly))
    }
    
    # update
    dt <- dt %>% bind_rows(data.frame(trt, week, last_wk_tir, rand_week, rand_eligible))
  }
  
  return(dt)
}


# #   AND
# #   at risk in the previous week. ie, lag(any_risk_flag) == 1
# elig_rand_2 = (lag(any_risk_flag, default = FALSE) == 1),
# elig_rand = elig_rand_1 * elig_rand_2,
# tir = rbeta(nrow(.), 10, 5)







# Setup for patient i from CDCES j
set.seed(13507)

data.frame(
  # One your of follo up
  study_week = 12:52) %>%
  mutate(
    # tide week of the month
    tide_week = (study_week - 1) %% 4 + 1,
    # default CPM week of the month (equally distributed)
    default_rmp_week = sample(1:4, 1),
    # RMP happening that week (Y/N), by default?
    default_rmp = as.numeric(tide_week == default_rmp_week),
    # risk flag (Y/N), i.e., if any risk flag = 1, otherwise = 0
    any_risk_flag = rbinom(nrow(.), 1, prob = .5), # can change p later
    # eligible for randomization if:
    #   week following a default RPM week,
    #     ie, tide_week == default_rmp_week + 1 (if 5 -> 1)
    #.    ie, tide_week == default_rmp_week %% 4 + 1
    elig_rand_1 = (tide_week == default_rmp_week %% 4 + 1),
    #   AND
    #   at risk in the previous week. ie, lag(any_risk_flag) == 1
    elig_rand_2 = (lag(any_risk_flag, default = FALSE) == 1),
    elig_rand = elig_rand_1 * elig_rand_2,
    tir = rbeta(nrow(.), 10, 5)
) %>% as_tibble()
