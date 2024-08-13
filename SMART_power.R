
library(tidyverse)
library(lmerTest)
library(gtsummary)


# Poulation parametres

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

# set.seed(4560798)
# generate_4t_pt()


set.seed(4560798)
dt <- lapply(1:100, \(x) {
  generate_4t_pt()
}) %>% bind_rows(.id = 'pid')

# compare weekly vs monthly TIRs on a LMM
lmer(last_wk_tir ~ trt + (1 | pid), data = dt %>% filter(rand_eligible == 1)) %>%
  tbl_regression()




expand_grid(id = 1, week = 1:48) %>% 
  mutate(elig_ctrl = case_when(
    week <= 12 ~ 1,
    week %% 4 == 0 ~ 1,
    .default = 0),
    
    at_risk = rbinom(nrow(.), 1, p=.20), # 20% of pt with bad glycemic control
    
    elig_rnd_1 = 0 + (week > 12),
    elig_rnd_2 = 0 + (week %% 4 == 1),
    elig_rnd_3 = at_risk,
    
    elig_rnd = elig_rnd_1 * elig_rnd_2 * elig_rnd_3) %>% 
  print(n=50)
  

















