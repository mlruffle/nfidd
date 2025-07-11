####Getting Set Up

options(repos = c(
  "CRAN" = "https://cloud.r-project.org",
  "stan-dev" = "https://stan-dev.r-universe.dev",
  "epiforecasts" = "https://epiforecasts.r-universe.dev"
))
install.packages("nfidd", dependencies = TRUE)

library(nfidd)

cmdstanr::install_cmdstan()

###Going Further Challenge 1 - using log normal rather than gamma
### from Probability distributions and parameter estimation

## load lognormal model from the nfidd package
mod <- nfidd_cmdstan_model("lognormal")
mod

lognormals <- rlnorm(100, meanlog = 1.57, sdlog = 0.28)

stan_data <- list(
  n = length(lognormals),
  y = lognormals
)
logn_fit <- nfidd_sample(mod, data = stan_data)
logn_fit

## Extract posterior draws
logn_posterior <- as_draws_df(logn_fit$draws())
head(logn_posterior)

## Generate posterior predictive samples
logn_ppc <- sapply(seq_along(lognormals), function(i) {
  rlnorm(n = length(lognormals),
         meanlog = logn_posterior$meanlog[i],
         sdlog = logn_posterior$sdlog[i])
})

## Plot posterior predictive check
ppc_dens_overlay(y = lognormals, yrep = logn_ppc)



###Going Further Challenge 2 - using gamma rather than log-normal
### from delay-distributions

mod <- nfidd_cmdstan_model("gamma")
mod

## Specify the time from onset to hospitalisation
df_onset_to_hosp <- df |>
  mutate(onset_to_hosp = hosp_time - onset_time) |>
  # exclude infections that didn't result in hospitalisation
  drop_na(onset_to_hosp)
## Use the data to sample from the model posterior
res <- nfidd_sample(
  mod,
  data = list(
    N = nrow(df_onset_to_hosp),
    y = df_onset_to_hosp$onset_to_hosp
  )
)

res$summary() #lp_ values are negative
res #ess not great



#### Simulating Censoring with Floor Division

### Formula rounds values down to the nearest interval boundary:
## Example: Daily censoring (interval = 1)
times <- c(0.2, 0.8, 1.3, 1.9, 2.1, 2.7)
censored <- floor(times / 1) * 1
# Result: [0, 0, 1, 1, 2, 2]
# demonstrates why the rounding works because our real data is 'times' but we observed 'results'

## Weekly censoring (interval = 7)
days <- c(2, 8, 12, 15, 20)
weekly_censored <- floor(days / 7) * 7
# Result: [0, 7, 7, 14, 14]

### For censoring reporting delays:
## Simulate continuous delays
true_delays <- rgamma(100, shape = 2, rate = 0.5)
## Apply daily censoring
daily_censored <- floor(true_delays / 1) * 1
## Apply weekly censoring
weekly_censored <- floor(true_delays / 7) * 7

## Generate example data
set.seed(123)
continuous_times <- runif(50, 0, 10)
daily_censored <- floor(continuous_times / 1) * 1
weekly_censored <- floor(continuous_times / 7) * 7

## Plot
data.frame(
  continuous = continuous_times,
  daily = daily_censored,
  weekly = weekly_censored
) |>
  tidyr::pivot_longer(everything(), names_to = "type", values_to = "time") |>
  ggplot(aes(x = time, fill = type)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 20) +
  facet_wrap(~type, ncol = 1) +
  labs(title = "Effect of Different Censoring Intervals")

###Example of weekly censored data
# Use the floor() function to round down to integers
df_dates_week <- df |>
  mutate(
    infection_time = floor(infection_time/7)*7,
    onset_time = floor(onset_time/7)*7,
    hosp_time = floor(hosp_time/7)*7
  )
head(df_dates_week)

df_dates_week <- df_dates_week |>
  mutate(
    incubation_period = onset_time - infection_time,
    onset_to_hosp = hosp_time - onset_time
  )
summary(df$onset_time - df$infection_time)
summary(df_dates_week$incubation_period)
