#----------------------------

# Prospective Power Analysis

#----------------------------

# this script shows you how to implement prospective power calculations
# using fake-data simulations

# load required packages
list.of.packages <- c("tidyverse", "here", "AER", "broom")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(here) # for file paths management
library(tidyverse) # for data manipulation and visualization
library(AER) # for running iv model
library(broom) # for cleaning statistical model output

##-------------------

## Standard Approach

##-------------------

# power formula for complete experiment
power_calculator <- function(mu_t, mu_c, sigma, alpha = 0.05, N) {
  lowertail <- (abs(mu_t - mu_c) * sqrt(N)) / (2 * sigma)
  uppertail <- -1 * lowertail
  beta <- pnorm(lowertail - qnorm(1 - alpha / 2), lower.tail = TRUE) + 1 - pnorm(uppertail - qnorm(1 - alpha / 2), lower.tail = FALSE)
  return(beta)
} 

# computing power of our agricultural experiment
power_calculator(mu_t = 15, mu_c = 14, sigma = 5, alpha = 0.05, N = 100)


##-----------------------------------------------------------

## Simulation-Based Approach: Health Effects of Forest Fires

##-----------------------------------------------------------

# generate IV data
generate_data_IV <- function(N,
                             param_fire,
                             sigma_u,
                             sigma_ee,
                             sigma_ep,
                             alpha,
                             gamma,
                             treatment_effect,
                             iv_strength,
                             ovb_intensity) {
  # first create day index
  data <- tibble(id = 1:N) %>%
    mutate(
      # create fire radiative power from a gamma distribution
      fire = rgamma(N, shape = param_fire[1], scale = param_fire[2])*50000,
      # create the unobserved variable
      u = rnorm(nrow(.), 0, sigma_u),
      # random noise for the emergency model
      e_e = rnorm(nrow(.), 0, sigma_ee),
      # random noise for the air pollution model
      e_p = rnorm(nrow(.), 0, sigma_ep),
      # create air pollution variable
      pollution = gamma + iv_strength * fire + ovb_intensity * u + e_p,
      # create emergency admission variable
      emergency = alpha + treatment_effect * pollution + ovb_intensity * u + e_e
    )
  
  return(data)
}

# set baseline parameters
baseline_param_IV <- tibble(
  N = 3650,
  param_fire = list(c(1, 2)),
  sigma_u = 1,
  sigma_ee = 100,
  sigma_ep = 5,
  alpha = 400,
  gamma = 30,
  treatment_effect = 1,
  iv_strength = 10/100000, 
  ovb_intensity = 0.4
)

# set seed
set.seed(42)

# simulate one dataset
data_fire <- baseline_param_IV %>%
  pmap_dfr(generate_data_IV)

# plot the density distribution of the main variables
data_fire %>%
  pivot_longer(cols = c(fire, pollution, emergency),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = case_when(metric == "fire" ~ "FRP (MV)",
                            metric == "pollution" ~ "PM10 (µg/m³)",
                            metric == "emergency" ~ "Emergency Admission (N)")) %>%
  ggplot(.,
         aes(x = value)) +
  geom_density() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), labels = scales::scientific) +
  facet_wrap( ~ metric, nrow = 1, scales = "free") +
  xlab("Value of the Variable") + ylab("") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# check effect of forest fires on air pollution
data_fire %>%
  mutate_at(vars(fire, pollution), ~ scale(.)) %>%
  lm(pollution ~ fire + u, data = .) %>%
  broom::tidy(., conf.int = TRUE) %>%
  filter(term == "fire")

# check effect of instrumented air pollution on emergency admission
data_fire %>%
  mutate_at(vars(fire, pollution), ~ scale(.)) %>%
  AER::ivreg(data = .,
             formula = log(emergency) ~ pollution | fire) %>%
  broom::tidy(., conf.int = TRUE) %>%
  filter(term == "pollution")

# create function to estimate IV model
estimate_IV <- function(data) {
  reg_IV <- AER::ivreg(data = data,
                       formula = emergency ~ pollution | fire)
  
  fstat_IV <- summary(reg_IV,
                      diagnostics = TRUE)$diagnostics["Weak instruments", "statistic"]
  
  reg_IV <- reg_IV %>%
    broom::tidy(., conf.int = TRUE) %>%
    mutate(fstat = fstat_IV) %>%
    filter(term == "pollution") %>%
    rename(p_value = p.value, se = std.error, conf_low = conf.low, conf_high = conf.high) %>%
    select(estimate, p_value, se, conf_low, conf_high, fstat) %>%
    
    return(reg)
}

# define compute_sim_IV
compute_sim_IV <- function(N,
                           param_fire,
                           sigma_u,
                           sigma_ee,
                           sigma_ep,
                           alpha,
                           gamma,
                           treatment_effect,
                           iv_strength,
                           ovb_intensity) {
  generate_data_IV(
    N = N,
    param_fire = param_fire,
    sigma_u = sigma_u,
    sigma_ee = sigma_ee,
    sigma_ep = sigma_ep,
    alpha = alpha,
    gamma = gamma,
    treatment_effect = treatment_effect,
    iv_strength = iv_strength,
    ovb_intensity = ovb_intensity
  ) %>%
    estimate_IV() %>%
    mutate(
      iv_strength = iv_strength,
      ovb_intensity = ovb_intensity,
      true_effect = treatment_effect
    )
} 

# set vector of iv strengths
vect_iv_strength <- c(1, 2, 3, 4, 5, 10) / 100000

# number of iterations
n_iter <- 500

# create parameters of the simulations
param_IV <- baseline_param_IV %>%
  select(-iv_strength) %>%
  crossing(vect_iv_strength) %>%
  rename(iv_strength = vect_iv_strength) %>%
  crossing(rep_id = 1:n_iter) %>%
  select(-rep_id)

# load simulations data
sim_IV <- readRDS(here::here("data", "data_sim_power_fire.rds"))

# function to summarise relevant metrics
summarise_simulations <- function(data) {
  data %>%
    mutate(significant = (p_value <= 0.05)) %>% 
    group_by(ovb_intensity, iv_strength) %>%
    summarise(
      power = mean(significant, na.rm = TRUE)*100, 
      type_m = mean(ifelse(significant, abs(estimate/true_effect), NA), na.rm = TRUE),
      mean_fstat = mean(fstat, na.rm = TRUE),
      mean_width_ci = mean(abs(conf_high - conf_low)),
      .groups	= "drop"
    ) %>% 
    ungroup()
} 

# run summarise_simulations
summary_sim_IV <- summarise_simulations(sim_IV)

# display simulation results
summary_sim_IV %>%
  select(-ovb_intensity) %>%
  mutate_at(vars(power, type_m, mean_width_ci), ~ round(., 1)) %>%
  mutate(mean_fstat = round(mean_fstat, 0)) %>%
  rename(
    "Strength of the IV" = iv_strength,
    "Power (%)" = power,
    "Average Type M Errror" = type_m,
    "Mean F-Statistic" = mean_fstat,
    "Mean Width of 95% CI" = mean_width_ci
  )