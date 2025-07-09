library(cmdstanr)
library(posterior)
library(ggplot2)
library(tidybayes)
library(dplyr)

# ---- 1. Simulate or load your case data ----
# Example data
# Replace with your actual data
set.seed(123)
days <- 0:141
infections <- rpois(length(days), lambda = 5)
inf_ts <- data.frame(day = days, infections = infections)

# ---- 2. Discretize the generation time ----
censored_delay_pmf <- function(rdist, max, ...) {
  breaks <- 0:max
  probs <- diff(sapply(breaks, function(b) mean(rdist(100000, ...)<b)))
  return(probs)
}

gen_time_pmf <- censored_delay_pmf(rgamma, max = 14, shape = 4, rate = 1)
gen_time_pmf <- gen_time_pmf[-1]
gen_time_pmf <- gen_time_pmf / sum(gen_time_pmf)

# ---- 3. Write the Stan model to a temporary file ----
stan_code <- '
functions {
  array[] real renewal(int I0, array[] real R, array[] real gen_time_pmf) {
    int max_gen_time = size(gen_time_pmf);
    int times = size(R);
    array[times] real I;
    array[times + 1] real full_I;

    full_I[1] = I0;
    for (i in 2:(times + 1)) {
      full_I[i] = 0;
    }

    for (t in 1:times) {
      real sum_infectious = 0;
      int first_index = (t - max_gen_time + 1 > 1) ? t - max_gen_time + 1 : 1;
      for (s in first_index:t) {
        int gen_index = t - s + 1;
        sum_infectious += full_I[s] * gen_time_pmf[gen_index];
      }
      full_I[t + 1] = sum_infectious * R[t];
      I[t] = full_I[t + 1];
    }
    return I;
  }
}

data {
  int n;                         // number of days
  int I0;                        // initial infections
  array[n] int obs;              // observed infections
  int gen_time_max;              // maximum generation time
  array[gen_time_max] real gen_time_pmf;  // generation time PMF
}

parameters {
  array[n] real<lower=0> R;      // reproduction number
}

transformed parameters {
  array[n] real infections = renewal(I0, R, gen_time_pmf);
}

model {
  R ~ normal(1, 1);              // prior
  obs ~ poisson(infections);    // likelihood
}
'

stan_file <- write_stan_file(stan_code)

# ---- 4. Compile the model ----
model <- cmdstan_model(stan_file)

# ---- 5. Prepare data for Stan ----
data_list <- list(
  n = nrow(inf_ts) - 1,
  obs = inf_ts$infections[-1],
  I0 = inf_ts$infections[1],
  gen_time_max = length(gen_time_pmf),
  gen_time_pmf = gen_time_pmf
)

# ---- 6. Fit the model ----
fit <- model$sample(
  data = data_list,
  parallel_chains = 4,
  iter_sampling = 500,
  iter_warmup = 500,
  refresh = 100
)
#fit <- nfidd_sample(model, data = data)

# ---- 7. Posterior draws and tidy plot (sampled Rt lines) ----

# Extract 100 sample trajectories of Rt
r_posterior <- fit |>
  gather_draws(R[infection_day]) |>
  ungroup() |>
  mutate(infection_day = infection_day - 1) |>  # shift so day 0 aligns with your data
  filter(.draw %in% sample(unique(.draw), 100))  # sample 100 trajectories

# Plot
ggplot(
  data = r_posterior,
  aes(x = infection_day, y = .value, group = .draw)
) +
  geom_line(alpha = 0.1) +
  labs(
    title = "Estimated Rt",
    subtitle = "Model: renewal equation from infections",
    x = "Day", y = "R(t)"
  ) +
  theme_minimal()
