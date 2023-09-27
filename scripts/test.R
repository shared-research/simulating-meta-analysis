# Test Simulations

library(metafor)
devtools::load_all()
seed <- 2023

# Equal-Effects Model -----------------------------------------------------

set.seed(seed)
k <- 1e3 # number of studies
n <- 1e3 # number of participants per group, per study
theta <- 0.3 # real effect size

sim <- make_data(k = k, nt = n, nc = n, es = theta)
sim <- sim_studies(sim$es, sim$nc, sim$nt, data = sim) # simulated data
res <- rma(yi = yi, vi = vi, method = "EE", data = sim)

check_sim(fixed = list(theta = theta), 
          results = list(theta = res$b[[1]]),
          name = "Equal-Effects Model")

# Random-effects Model ----------------------------------------------------

set.seed(seed)
k <- 1e3 # number of studies
n <- 1e3 # number of participants per group, per study
mu <- 0.3 # real effect size
tau2 <- 0.2 # the effect size heterogeneity

sim <- make_data(k = k, nc = n, nt = n, mu = mu)

# simulate the random-effects adjustment
sim$deltai <- rnorm(k, 0, sqrt(tau2))

# adding the by-study adjustment to the average effect
sim$es <- sim$mu + sim$deltai 

# now we are using mu_deltai and no longer theta
sim <- sim_studies(sim$es, sim$nc, sim$nt, data = sim)
res <- rma(yi = yi, vi = vi, method = "REML", data = sim)

check_sim(fixed = list(theta = theta,
                       tau2 = tau2), 
          results = list(theta = res$b[[1]],
                         tau2 = res$tau2),
          name = "Random-effects Model")

# Random-effects model fixing I2 ------------------------------------------

set.seed(seed)
theta <- 0.3
k <- 1e3 # number of studies
mu <- 0.3 # real effect size
n <- 1e3 # sample size per group, per study
I2 <- 0.6 # desired I2 value

v <- 1/n + 1/n # typical within study variance
tau2 <- -((I2*v)/(I2 - 1))

sim <- make_data(k = k, nc = n, nt = n, mu = mu)

sim$deltai <- rnorm(k, 0, sqrt(tau2))
sim$es <- sim$mu + sim$deltai 
sim <- sim_studies(sim$es, sim$nc, sim$nt, data = sim)
res <- rma(yi, vi, method = "REML", data = sim)

check_sim(fixed = list(theta = theta,
                       tau2 = tau2,
                       I2 = I2*100), 
          results = list(theta = res$b[[1]],
                         tau2 = res$tau2,
                         I2 = res$I2),
          name = "Random-effects Model fixing I2")

# Meta-regression with categorical moderator ------------------------------

set.seed(seed)
k <- 1e3 # the total number of studies
b0 <- 0.1 # intercept, the effect size of the lab-based studies
b1 <- 0.2 # the difference between the two levels of the moderator
tau2r <- 0.1 # the residual heterogeneity
n <- 1e3 # the sample size per group, per study

sim <- make_data(k = k, nc = n, nt = n, exp = rep(c("lab", "online"), each = k/2))

sim$deltai <- rnorm(k, 0, sqrt(tau2r)) # the by-study residual adjustment
sim$es <- b0 + sim$deltai + b1*ifelse(sim$exp == "lab", 0, 1)

sim <- sim_studies(sim$es, sim$nc, sim$nt, data = sim)
res <- rma(yi, vi, mods = ~exp, method = "REML", data = sim)

check_sim(fixed = list(b0 = b0,
                       b1 = b1,
                       tau2r = tau2r), 
          results = list(b0 = res$b[[1]],
                         b1 = res$b[[2]],
                         tau2r = res$tau2),
          name = "Meta-regression with categorical predictor")

# Meta-regression with numerical moderator --------------------------------

set.seed(seed)
k <- 1e3 # number of studies
b0 <- 0.3 # the intercept i.e., average yi when x is 0
b1 <- 0.1 # the beta1 i.e., the increase in yi for an increase in 1 year
tau2r <- 0.1 # the residual heterogeneity after including x1
n <- 1e3 # number of participants per study, per group
x1 <- runif(k, 20, 40) # random mean-age for each study
x10 <- x1 - mean(x1) # centering the age

sim <- make_data(k = k, nc = n, nt = n, age0 = x10)

sim$deltai <- rnorm(k, 0, sqrt(tau2r))

sim$es <- b0 + sim$deltai + b1*sim$age0

sim <- sim_studies(sim$es, sim$nc, sim$nt, data = sim)

res <- rma(yi, vi, mods = ~age0, method = "REML", data = sim)

check_sim(fixed = list(b0 = b0,
                       b1 = b1,
                       tau2r = tau2r), 
          results = list(b0 = res$b[[1]],
                         b1 = res$b[[2]],
                         tau2r = res$tau2),
          name = "Meta-regression with numerical predictor")


# Meta-regression with numerical predictor, fixing R2 ---------------------

set.seed(seed)
k <- 1e3 # the number of studies
r2 <- 0.2 # the desired r2 value
tau2 <- 0.3 # the overall tau2
b0 <- 0.3 # the intercept i.e., average yi when x1 is 0
b1_2 <- tau2 * r2 # the beta1^2 i.e., the increase in yi for an increase in 1 year
b1 <- sqrt(b1_2) # b1_2 is squared, back to the original scale
tau2r <- tau2 - b1_2 # the residual heterogeneity after including x1
n <- 1e3 # number of participants per study, per group
x1 <- runif(k, 20, 40) # random mean-age for each study

sim <- make_data(k = k, nt = n, nc = n, age = x1)

sim$deltai <- rnorm(k, 0, sqrt(tau2r))
sim$age0 <- scale(sim$age, center = TRUE, scale = TRUE) # standardize the moderator
sim$es <- b0 + sim$deltai + b1*sim$age0

sim <- sim_studies(sim$es, sim$nt, sim$nc, data = sim)
res <- rma(yi, vi, mods = ~age0, method = "REML", data = sim)

check_sim(fixed = list(b0 = b0,
                       b1 = b1,
                       tau2r = tau2r,
                       r2 = r2*100), 
          results = list(b0 = res$b[[1]],
                         b1 = res$b[[2]],
                         tau2r = res$tau2,
                         r2 = res$R2),
          name = "Meta-regression fixing R2")

# For the R2 we tried to repeat the simulation and check the
# distribution of simulated values

nsim <- 1e3
k <- 100 # the number of studies
r2 <- 0.2 # the desired r2 value
tau2 <- 0.3 # the overall tau2
b0 <- 0.3 # the intercept i.e., average yi when x1 is 0
b1_2 <- tau2 * r2 # the beta1^2 i.e., the increase in yi for an increase in 1 year
b1 <- sqrt(b1_2) # b1_2 is squared, back to the original scale
tau2r <- tau2 - b1_2 # the residual heterogeneity after including x1
n <- 100 # number of participants per study, per group

r2s <- rep(0, nsim)

for(i in 1:nsim){
  x1 <- runif(k, 20, 40) # random mean-age for each study
  sim <- make_data(k = k, nt = n, nc = n, age = x1)
  sim$deltai <- rnorm(k, 0, sqrt(tau2r))
  sim$age0 <- scale(sim$age, center = TRUE, scale = TRUE) # standardize the moderator
  sim$es <- b0 + sim$deltai + b1*sim$age0
  sim <- sim_studies(sim$es, sim$nt, sim$nc, data = sim)
  res <- rma(yi, vi, mods = ~age0, method = "REML", data = sim)
  r2s[i] <- res$R2
}

title <- sprintf("Median $R^2 = %.3f$ (%s simulations)", median(r2s), nsim)
hist(r2s, main = latex2exp::TeX(title))

# Power Analysis ----------------------------------------------------------

# for the power analysis we simulate the type-1 error rate that should be
# ~ 0.05

set.seed(seed)
nsim <- 1000 # number of simulations per condition
k <- c(20, 30, 100) # number of studies
delta <- 0 # the average effect size
tau2 <- 0.2 # the heterogeneity
navg <- 20 # average sample size per study
nmin <- 10 # minimum sample size per study
alpha <- 0.05 # the alpha level

# creating all combinations
sim_grid <- tidyr::expand_grid(k, mu = delta, tau2, navg, nmin, nsim)

# apply the simulation to all combinations
res <- purrr::pmap(sim_grid, do_sim)

sim_grid$type1error <- unlist(res)
sim_grid

# 3-level meta-analysis ---------------------------------------------------

set.seed(seed)
i <- 500 # Number of studies (level 3)
j <- 5 # Number of effect sizes within each study (level 2)
tau2 <- 0.3 # heterogeneity between studies
omega2 <- 0.1 # heterogeneity within studies
icc <- tau2 / (tau2 + omega2) # real ICC
n <- 1e3 # sample size for each group, within each study
mu <- 0.3 # real average effect size

delta_i <- rnorm(i, 0, sqrt(tau2))
zeta_ij <- rnorm(i*j, 0, sqrt(omega2))

sim <- tidyr::expand_grid(
  study = 1:i,
  effect = 1:j,
  nt = n,
  nc = n,
  mu = mu
)

sim$delta_i <- delta_i[sim$study]
sim$zeta_ij <- zeta_ij

sim$es <- sim$mu + sim$delta_i + sim$zeta_ij
sim <- sim_studies(es = sim$es, nt = sim$nt, nc = sim$nc, data = sim)
res <- rma.mv(yi, vi, random = ~1|study/effect, data = sim, sparse = TRUE)

check_sim(fixed = list(mu = mu,
                       tau2 = tau2,
                       omega2 = omega2), 
          results = list(mu = res$b[[1]],
                         tau2 = res$sigma2[1],
                         omega2 = res$sigma2[2]),
          name = "3level model")

# Multivariate meta-analysis ----------------------------------------------

set.seed(seed)

k <- 500 # number of studies
p <- 3 # number of outcomes per study
rho <- 0.7 # population level correlation between outcomes
r <- 0.5 # correlation between sampling errors
n <- 1e3 # number of participants per group per study

mus <- c(0.1, 0, 1) # population level real effect sizes
taus2 <- c(0.3, 0.3, 0.3) # population level real tau2s

# creating the dataframe structure
sim <- tidyr::expand_grid(
  id = 1:k,
  p = 1:p,
  nt = n,
  nc = n
)

# adding the real effects and the random-effect adjustments
sim$mu <- rep(mus, k)
sim$delta_i <- sim_T(k, p, taus2, rho)

# average effect + random-effect adjustment
sim$es <- sim$mu + sim$delta_i

# splitting to list of studies for improving the simulation speed
siml <- split(sim, sim$id)

# simulating the effect sizes for each study
resl <- lapply(siml, function(x) sim_study_m(x$es, x$nt, x$nc, r))
res <- dplyr::bind_rows(resl) # combine everything
sim <- cbind(sim, res) # append to the sim dataframe

# create the block variance-covariance matrix
# cluster is the study (i.e. each "block")
# obs identify each outcome (i.e., each "block" number of rows/columns)
# r is the fixed sampling errors correlation
V <- vcalc(vi = vi, cluster = id, obs = p, rho = r, data = sim)

# outcome to character with a meaningful prefix
sim$p <- paste0("outcome", sim$p)

# fitting the model
res <- rma.mv(yi, 
              V, 
              mods = ~ 0 + p, 
              random = ~p|id, 
              data = sim, 
              struct = "CS", 
              sparse = TRUE)

fixed <- list(mu1 = mus[1],
             mu2 = mus[2],
             mu3 = mus[3],
             tau2 = taus2[1], # CS thus the same tau2 for each outcome
             rho = rho)

results <- list(mu1 = res$b[[1]],
                mu2 = res$b[[2]],
                mu3 = res$b[[3]],
                tau2 = res$tau2,
                rho = res$rho)

check_sim(fixed, results, name = "Multivariate model")
