# load function to:
# (i) run a complete randomized experiment,
# (ii) estimate the average treatment effect
# (iii) compute neyman's 95% confidence interval
# see chapter 6 of imbens and rubin (2015) for further details
# the function's argument takes a science table

function_complete_experiment <- function(data) {
  # allocate treatment and express observed yields
  data <-  data %>%
    mutate(w = sample(w, replace = FALSE),
           y_obs = w * y_1 + (1 - w) * y_0)
  
  # compute estimate for ate
  estimate_ate <- data %>%
    group_by(w) %>%
    summarise(group_mean = mean(y_obs)) %>%
    summarise(group_mean[2] - group_mean[1]) %>%
    pull()
  
  # compute associated standard error
  standard_error <- data %>%
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
  tibble(
    estimate_ate = estimate_ate,
    lower_95_ci = estimate_ate - 1.96 * standard_error,
    upper_95_ci = estimate_ate + 1.96 * standard_error
  )
}