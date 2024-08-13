

library(tidyverse)
library(lmerTest)
library(gtsummary)
library(gtExtras)
library(gt)

skip.seed.streams <- function(n) {
  x <- .Random.seed
  for (i in seq_len(n))
    x <- parallel::nextRNGStream(x)
  assign('.Random.seed', x, pos=.GlobalEnv)
}

# start block
RNGkind("L'Ecuyer-CMRG")

sim_4t_microrand <- function(N, capacity = Inf, return_dt_sim = FALSE) {
  # N: scalar or vector of sample sizes to be simulated
  
  dt_aux0 <- lapply(N, \(x){
    data.frame(sample_size = x, pid = 1:x)
  }) %>% bind_rows()
  
  # Generate baseline data
  dt_aux <- dt_aux0 %>% 
    group_split(sample_size, pid) %>% 
    lapply(\(x) {
      dt <- data.frame(
        sample_size = x$sample_size[1],
        pid = x$pid[1],
        study_week = 13:52)
      
      bind_rows(
        dt %>% mutate(set = "A1", sd_g = .20, sd_e = .20/2, b3 = .005),
        dt %>% mutate(set = "B1", sd_g = .16, sd_e = .16/2, b3 = .005),
        dt %>% mutate(set = "A2", sd_g = .20, sd_e = .20/2, b3 = .0025),
        dt %>% mutate(set = "B2", sd_g = .16, sd_e = .16/2, b3 = .0025),
        dt %>% mutate(set = "A3", sd_g = .20, sd_e = .20/2, b3 = .001),
        dt %>% mutate(set = "B3", sd_g = .16, sd_e = .16/2, b3 = .001),
        dt %>% mutate(set = "A4", sd_g = .20, sd_e = .20/2, b3 = .0005),
        dt %>% mutate(set = "B4", sd_g = .16, sd_e = .16/2, b3 = .0005),
        dt %>% mutate(set = "A5", sd_g = .20, sd_e = .20/2, b3 = .0001),
        dt %>% mutate(set = "B5", sd_g = .16, sd_e = .16/2, b3 = .0001),
        dt %>% mutate(set = "C1", sd_g = .20, sd_e = .20/2, b3 = 0),
        dt %>% mutate(set = "C2", sd_g = .16, sd_e = .16/2, b3 = 0)) %>%
        group_by(set) %>% 
        mutate(
          time = study_week - 13,
          tide_week = (study_week - 12 + sample(1:4, 1) - 1) %% 4 + 1,
          
          gamma = rnorm(1, 0, sd_g^2),
          error = rnorm(length(study_week), 0, sd_e^2),
          
          treatment = 'Default',
          TIR = .75 + gamma - .005*time + error,
          iter = cumsum(as.numeric(tide_week == 1)),
          rand = 0) %>% 
        ungroup()
    }) %>%
    bind_rows() %>%
    filter(tide_week <= 2)
  
  # Iterate study week
  dt_sim <- dt_aux %>% 
    group_split(sample_size, set) %>% 
    lapply(\(x) {
      aux_risk <- x %>% 
        filter(tide_week == 1, TIR < .65) 
      
      if (nrow(aux_risk) == 0) return()
      
      min_k <- min(aux_risk$iter)
      dt_iter <- x
      for (k in min_k:10) {
        dt_iter <- dt_iter %>% 
          mutate(
            elig = ifelse(tide_week == 1, as.numeric(TIR < .65), NA),
            rand_trt = ifelse((iter == k)*elig == 1, sample(c(0, 1), nrow(.), replace = T), NA)) %>% 
          group_by(pid) %>% 
          mutate(
            rand_trt = ifelse(is.na(rand_trt), '', rand_trt),
            rand_trt = paste0(lag(rand_trt, n = 1, default = ''),
                              lag(rand_trt, n = 2, default = '')),
            rand_trt = ifelse(rand_trt == '', NA, as.numeric(rand_trt))) %>%
          ungroup() %>%
          mutate(
            rand = rand + as.numeric(!is.na(rand_trt)),
            treatment = ifelse(!is.na(rand_trt) & rand_trt == 1, 'Addon', treatment),
            rand_TIR = .75 + gamma - .005*time + b3*time*rand_trt + error,
            TIR = ifelse(!is.na(rand_trt) & rand_trt == 1, rand_TIR, TIR))
        
        # print(dt_iter, n=100)
      }
      
      dt_iter %>%
        mutate(
          sample_size = x$sample_size[1],
          set = x$set[1],
          rand_trt = NULL,
          rand_TIR = NULL) %>% 
        filter(iter > 1, tide_week == 1, rand == 1)
    }) %>%
    bind_rows()
  
  # Impose capacity
  if (is.finite(capacity)) {
    dt_sim <- dt_sim %>% 
      # filter(sample_size==200, set=='B', time == 38) %>% 
      arrange(sample_size, set, time, treatment) %>% 
      group_by(sample_size, set, time) %>% 
      mutate(
        addon_cum = cumsum(as.numeric(treatment == 'Addon')),
        addon_cum = ifelse(treatment == 'Addon', addon_cum, NA)) %>% 
      ungroup() %>% 
      mutate(
        exceed_cap = ifelse(addon_cum > capacity, 1, 0),
        treatment = ifelse(!is.na(exceed_cap) & exceed_cap == 1, 'Default', treatment),
        TIR = ifelse(!is.na(exceed_cap) & exceed_cap == 1, .75 + gamma - .005*time + error, TIR)) %>% 
      arrange(sample_size, set, pid, time)
  }
  
  # Summary
  dt_summ <- dt_sim %>% 
    count(sample_size, set, time, treatment) %>% 
    pivot_wider(names_from = treatment, values_from = n) %>% 
    mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .)),
           atRisk = Addon + Default)

  if(return_dt_sim) return(dt_sim)
  
  # Fit LMM
  aux_time <- seq(8, 40, 8)
  safe_lmer <- purrr::possibly(purrr::quietly(lmer), otherwise = NULL)
  safe_broom_tidy <- purrr::possibly(\(m) {
    broom.mixed::tidy(m) %>%
      filter(term == 'time:trt') %>% 
      select(-group)},
    otherwise = tibble())
  
  dt_reg <- lapply(aux_time, \(t) {
    dt_mod = dt_sim %>% 
      filter(time <= t) %>% 
      mutate(trt = as.numeric(treatment == 'Addon')) %>% 
      nest_by(sample_size, set, b3) %>% 
      mutate(
        up_to_time = t,
        n_clust = length(unique(data$pid)),
        n_obs = nrow(data),
        mod = list(safe_lmer(TIR ~ time + time*trt + (1|pid), data = data)$result)) %>% 
      reframe(up_to_time, n_clust, n_obs, safe_broom_tidy(mod))
  }) %>% bind_rows() %>% filter(!is.na(term))
  
  return(list(lmer_reg = dt_reg, summary = dt_summ))
}


# sim_4t_microrand(N = 50)
# sim_4t_microrand(N = c(50, 100, 200))

# dt_sim = sim_4t_microrand(N = c(50, 100, 200), return_dt_sim = T)
# 
# dt_summ <- dt_sim %>% 
#   count(sample_size, set, time, treatment) %>% 
#   pivot_wider(names_from = treatment, values_from = n) %>% 
#   mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .)),
#          atRisk = Addon + Default)
# 
# dt_sim %>% 
#   filter(sample_size==200, set=='B', time == 38) %>% 
#   arrange(sample_size, set, time, treatment) %>% 
#   group_by(sample_size, set, time) %>% 
#   mutate(
#     addon_cum = cumsum(as.numeric(treatment == 'Addon')),
#     addon_cum = ifelse(treatment == 'Addon', addon_cum, NA)) %>% 
#   ungroup() %>% 
#   mutate(
#     exceed_cap = ifelse(addon_cum > capacity, 1, 0),
#     treatment = ifelse(!is.na(exceed_cap) & exceed_cap == 1, 'Default', treatment),
#     TIR = ifelse(!is.na(exceed_cap) & exceed_cap == 1, .75 + gamma - .005*time + error, TIR)) %>% 
#   arrange(sample_size, set, pid, time)
# 
# 
# dt_summ %>% 
#   filter(sample_size==200, set=='B', time == 38)





set.seed(1234556)
set.seed(365342)

N = c(25, 50, 100, 200)
B = 2000
# N = c(25, 50)
# B = 10

skip.seed.streams(10)
list_reg_all <- seq(1, B) %>% 
  parallel::mclapply(\(x) {
    sim_list <- sim_4t_microrand(N = N, capacity = 10)
    
    list(sim_list$lmer_reg %>% mutate(r = x),
         sim_list$summary %>% mutate(r = x))
  }, mc.cores = 10, mc.set.seed = T)

# saveRDS(list_reg_all, 'list_reg_all_B2000_cap10_v3.rds')


list_reg_all <- readRDS("list_reg_all_B2000_cap10_v3.rds")

dt_reg_all <- lapply(list_reg_all, '[[', 1) %>% bind_rows() %>%
  filter(complete.cases(.))

dt_risk_all <- lapply(list_reg_all, '[[', 2) %>% bind_rows()

dt_reg_all %>%
  group_by(sample_size, set, up_to_time) %>%
  count() %>%
  arrange(up_to_time) %>% 
  print(n=500)

# compute metrics
dt_res <- dt_reg_all %>% 
  filter(term == 'time:trt') %>% 
  group_by(sample_size, set, up_to_time) %>% 
  add_count(name = 'n_boot') %>% 
  mutate(bias = estimate - b3,
         mse = bias^2,
         empirical_stderr = estimate - mean(estimate, na.rm = T),
         prop_rej = as.numeric(p.value < 0.05),
         ll = estimate - 1.96*std.error,#/sqrt(n),
         ul = estimate + 1.96*std.error,#/sqrt(n),
         cover = as.numeric((ll <= b3)*(b3 <= ul))) %>% 
  summarise(across(c(n_boot, bias, empirical_stderr, mse, cover, prop_rej),
                   ~ mean(., na.rm = T)), .groups = 'drop')

# dt_res %>% 
#   arrange(set, sample_size) %>% 
#   print(n = 100)

dt_risk_all %>% 
  group_by(sample_size, set, time) %>% 
  summarise(across(c(Default, Addon, atRisk),
                   list(m = mean, sd = sd)), .groups = 'drop') %>% 
  select(sample_size, set, time,
         `Avg. # addon RPM` = Addon_m,
         `Avg. # at risk` = atRisk_m) %>% 
  pivot_longer(-c(sample_size, set, time)) %>% 
  ggplot(aes(x = time, y = value, color = set, linetype = name)) +
  geom_hline(yintercept = 10, linetype = 'dashed') +
  annotate('text', x = 4, y = 12, label = 'capacity') +
  geom_line() +
  scale_x_continuous(limits = c(1,39)) +
  labs(x = 'Time (weeks) after randomization start',
       y = 'Average number of patients at risk',
       color = 'Scenario',
       linetype = '') +
  facet_wrap(~ sample_size) +
  theme_bw() +
  theme(legend.position = 'top')


dt_res %>% 
  filter(n_boot >= 50) %>%
  # filter(str_detect(set, 'A')) %>% 
  mutate(set_group = str_sub(set, 1, 1),
         set_nbr = str_sub(set, 2, 2)) %>% 
  ggplot(aes(x = up_to_time, y = prop_rej, color = set_nbr)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = .8, linewidth = .5, linetype = 'dashed') +
  facet_grid(set_group ~ sample_size) +
  scale_color_brewer(palette = 'Dark2') +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = seq(8, 40, 8)) +
  labs(x = 'Study week', y = 'Power') +
  theme_bw()
  
  
  
  
  
  
  
  
  
  








