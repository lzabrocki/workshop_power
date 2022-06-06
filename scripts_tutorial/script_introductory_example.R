#--------------------------------------------------

# Why P-Values Are Misleading in Low Power Studies

#--------------------------------------------------

# this script shows you why statistically significant estimates
# are inflated in low power studies

# load required packages
list.of.packages <- c("tidyverse", "here")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(here) # for file paths management
library(tidyverse) # for data manipulation and visualization

##----------------------

## Simulating the Data

##----------------------

# function to create the science table
function_science_table <- function(sample_size) {
  # define field index and y(0)
  tibble(field = c(1:sample_size),
         y_0 = rnorm(sample_size, 14, 5)) %>%
    # define y(1) by adding treatment effect size
    mutate(y_1 = y_0 + rnorm(sample_size, 1, 3)) %>%
    # round values of potential outcomes
    mutate_at(vars(y_0, y_1), ~ round(., 1))
}

# set seed
set.seed(42)

# create a science table for 100 fields
data_science <- function_science_table(50)

# plot the distributions of potential outcomes
data_science %>%
  rename("Y(0)" = y_0, "Y(1)" = y_1) %>%
  pivot_longer(cols = -c(field),
               names_to = "Potential Outcomes",
               values_to = "values") %>%
  group_by(`Potential Outcomes`) %>%
  mutate(mean_outcomes = mean(values)) %>%
  ggplot(. , aes(x = `Potential Outcomes`, y = values, colour = `Potential Outcomes`)) +
  geom_boxplot()+
  geom_jitter(alpha = 0.8) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  xlab("") + ylab(expression(paste("Yields of Watermelons (", "kg/", m^2, ")"))) +
  theme_minimal() +
  labs(colour = "Potential Outcomes:")

# compute the true value of the ate
ate <- data_science %>%
  mutate(tau = y_1 - y_0) %>%
  summarise(mean(tau) %>% round(., 1)) %>%
  pull()

# display ate value
ate

##-----------------------------------------

## Running One Iteration of the Experiment

##-----------------------------------------

# add the treatment vector to the science table
data_science <- data_science %>%
  mutate(w = c(rep(0, n() / 2), rep(1, n() / 2))) 

# load one iteration of the experiment
data_science_one_iteration <- readRDS(here::here("data", "data_science_one_iteration.rds"))

# express observed yields
data_science_one_iteration <- data_science_one_iteration %>%
  mutate(y_obs = w * y_1 + (1 - w) * y_0)

# compute estimate for ate
estimate_ate <- data_science_one_iteration %>%
  group_by(w) %>%
  summarise(group_mean = mean(y_obs)) %>%
  summarise(group_mean[2] - group_mean[1]) %>%
  pull()

# compute associated standard error
standard_error <- data_science_one_iteration %>%
  group_by(w) %>%
  mutate(
    n_obs = n(),
    average_y_obs = mean(y_obs),
    squared_difference = (y_obs - average_y_obs) ^
      2,
    variance_group = sum(squared_difference) / (n_obs - 1)
  ) %>%
  summarise(variance_component = mean(variance_group / n_obs)) %>%
  summarise(variance = sum(variance_component),
            standard_error = sqrt(variance)) %>%
  pull(standard_error)

# return ate estimate and 95% confidence interval
results_one_iteration <- tibble(
  estimate_ate = estimate_ate,
  lower_95_ci = estimate_ate - 1.96 * standard_error,
  upper_95_ci = estimate_ate + 1.96 * standard_error
) %>%
  mutate_all( ~ round(., 1))

# display results
results_one_iteration

##---------------------------------------

## Replicating Many Times the Experiment

##---------------------------------------

# load 1000 iterations of the experiment
data_results_1000_experiments <- readRDS(here::here("data", "data_results_1000_experiments.rds"))

# add iteration index
data_results_1000_experiments <- data_results_1000_experiments %>%
  mutate(iteration = 1:1000) %>%
  mutate(significant = ifelse(significant == 1, "True", "False"))

# plot estimates
ggplot(
    data_results_1000_experiments,
    aes(x = iteration, y = estimate_ate, color = significant)
  ) +
  geom_hline(yintercept = mean(data_results_1000_experiments$estimate_ate),
             size = 0.25,
             colour = "black") +
  geom_point(shape = 16,
             size = 1,
             alpha = 0.8) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 14)) +
  xlab("Iteration") + ylab("Estimate") +
  theme_minimal() +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3)))


# alternative graph with confidence intervals
data_results_1000_experiments %>%
  slice(1:100) %>%
  ggplot(
    .,
    aes(
      x = iteration,
      colour = significant,
      y = estimate_ate,
      ymin = lower_95_ci,
      ymax = upper_95_ci
    )
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_hline(yintercept = ate, colour = "black", linetype = "dashed") +
  geom_pointrange(lwd = 0.3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 15)) +
  coord_flip() +
  xlab("") + ylab("Estimate") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


##---------------------------------------------------------------

## Defining and Computing Statistical Power, Type M and S Errors

##---------------------------------------------------------------

# retrieve standard error of the first experiment iteration
s <- (results_one_iteration$upper_95_ci - results_one_iteration$estimate_ate) / 1.96

# compute power using its definition
(pnorm(-1.96 - 1.3 / s) + 1 - pnorm(1.96 - 1.3 / s)) * 100

# compute power using simulation results
data_results_1000_experiments %>%
  summarise(statistical_power = round(sum(significant == "True") / n() * 100, 0))

# compute average type m error
data_results_1000_experiments %>%
  filter(significant == "True") %>%
  summarise(exageration_factor = mean(abs(estimate_ate) / ate) %>% round(., 1))

# compute probability to make a type s error
data_results_1000_experiments %>%
  filter(significant == "True") %>%
  summarise(probability_type_s_error = round(sum(estimate_ate < 0) / n() * 100, 1))

##-------------------------------------------

## How Type M and S Errors Evolve with Power

##-------------------------------------------

# load simulations of the experiment
data_power_sample_size_sim <- readRDS(here::here("data", "data_power_sample_size_sim.rds"))

# compute power, type m and s errors
summary_power_sample_size_sim <- data_power_sample_size_sim %>%
  group_by(sample_size) %>%
  summarise(
    power = mean(significant == 1, na.rm = TRUE) * 100,
    type_m = mean(ifelse(
      significant == 1, abs(estimate_ate / ate), NA
    ), na.rm = TRUE),
    type_s = sum(ifelse(
      significant == 1, sign(estimate_ate) != sign(ate), NA
    ), na.rm = TRUE) / n() * 100
  )

# make the graph
summary_power_sample_size_sim %>%
  pivot_longer(cols = -c(sample_size, type_s),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = ifelse(metric == "power", "Power (%)", "Type M Error")) %>%
  ggplot(.,
         aes(x = sample_size, y = value)) +
  geom_point() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 14)) +
  facet_wrap( ~ metric, nrow = 1, scales = "free_y") +
  xlab("Sample Size") + ylab("") +
  theme_minimal()
