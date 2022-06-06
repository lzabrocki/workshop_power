#------------------------------

# Retrospective Power Analysis

#------------------------------

# this script shows you how to implement prospective power calculations
# using fake-data simulations

# load required packages
list.of.packages <- c("tidyverse", "here", "retrodesign")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(here) # for file paths management
library(tidyverse) # for data manipulation and visualization
library(retrodesign) # for retrospective power analysis

##-------------------

## Heat and Learning

##-------------------

# compute power for a reduction of 25% of the point estimate
retro_design(A = -0.0023*0.75, s = 0.00042, alpha = 0.05)

# compute power for lower bound of 95% ci
retro_design(A = -0.0015,
             s = 0.00042,
             alpha = 0.05) %>% as_tibble() %>% mutate_all( ~ round(., 2)) %>% kable(align = c(rep("c", 4)))

# compute power for an effect size of -0.1%
retro_design(A = -0.0003,
             s = 0.00042,
             alpha = 0.05)

# compute the power, type m and s errors for a range of effect sizes
data_retro <-
  retro_design(as.list(seq(-0.0023, -0.0003, by = 0.0001)), 0.00042) %>%
  unnest(cols = c(effect_size, power, type_s, type_m)) %>%
  mutate(power = power * 100,
         type_s = type_s * 100) %>%
  rename(
    "Statistical Power (%)" = power,
    "Type-S Error (%)" = type_s,
    "Type-M Error (Exaggeration Ratio)" = type_m
  ) %>%
  pivot_longer(
    cols = -c(effect_size),
    names_to = "statistic",
    values_to = "value"
  )

# make graph
ggplot(data_retro, aes(x = effect_size, y  = value)) +
  geom_line() +
  geom_vline(xintercept = -0.0023, linetype = "dashed") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  facet_wrap( ~ statistic, scales = "free") +
  xlab("Hypothetical True Effect Size") + ylab("") +
  theme_minimal()

##-------------------------------------

## Price Elasticity of Gasoline Demand

##-------------------------------------

# compute power for an effect size of -0.31
retro_design(A = -0.31,
             s = 0.39,
             alpha = 0.05) %>% as_tibble() 


