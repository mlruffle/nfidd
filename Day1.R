#Getting Set Up

options(repos = c(
  "CRAN" = "https://cloud.r-project.org",
  "stan-dev" = "https://stan-dev.r-universe.dev",
  "epiforecasts" = "https://epiforecasts.r-universe.dev"
))
install.packages("nfidd", dependencies = TRUE)

library(nfidd)

cmdstanr::install_cmdstan()
