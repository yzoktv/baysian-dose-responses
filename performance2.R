# ==============================
# BAYESIAN DOSE-RESPONSE PIPELINE
# ==============================
library(tidyverse)
library(simsurv)
library(rstan)
library(survival)
library(purrr)
library(splines)

# ~~~~~ USER CONFIGURATION ~~~~~
set.seed(123)
n <- 200
n_sim <- 50
scenarios <- c("linear", "sigmoid", "upward", "downward")
time_points <- c(2, 5, 8)
J <- 4
dose_levels <- c(2, 5, 8)

# ---- Data Generation ----
generate_data <- function(n, scenario) {
  f <- switch(scenario,
              linear = function(d) 0.1*d,
              upward = function(d) 0.01 * d + 0.01 * d^2,
              downward = function(d) 0.5 * d * exp(-0.2 * d),
              sigmoid = function(d) 1 / (1 + exp(5 - d)) - 0.01)
  
  dose <- numeric(0)
  while (length(dose) < n) {
    new_samples <- rgamma(n, shape = 2, rate = 0.4)
    dose <- c(dose, new_samples[new_samples < 20])
  }
  dose <- dose[1:n]
  df <- data.frame(dose = dose, lp = f(dose))
  sim_event <- simsurv(dist = "weibull",
                       lambdas = 0.1,
                       gammas = 1.5,
                       x = df,
                       betas = c(lp = 1),
                       maxt = 100)
  censor_time <- rexp(n, rate = 0.1)
  obs_time <- pmin(sim_event$eventtime, censor_time)
  status <- as.integer(sim_event$eventtime <= censor_time)
  return(data.frame(dose = dose, time = obs_time, status = status))
}

# ---- Model Fitting ----
# ---- Model Fitting ----
fit_models <- function(df, J = 4) {
  cox_linear <- coxph(Surv(time, status) ~ dose, data = df)
  cox_spline <- coxph(Surv(time, status) ~ ns(dose, df = 2), data = df)
  cox_pspline <- coxph(Surv(time, status) ~ pspline(dose, df = 2), data = df)
  
  knots <- quantile(df$dose, probs = c(0.1, 0.3, 0.5, 0.7, 0.9))
  B <- bs(df$dose, knots = knots, degree = 2, intercept = TRUE)
  
  # Construct time_knots of length J + 1
  if (sum(df$status == 1) > 0) {
    tq <- quantile(df$time[df$status == 1], probs = seq(0, 1, length.out = J))
    time_knots <- unique(c(0, tq))  # ensure starts at 0
    if (length(time_knots) < (J + 1)) {
      time_knots <- c(time_knots, max(df$time))  # pad if needed
    }
  } else {
    time_knots <- seq(0, max(df$time), length.out = J + 1)
  }
  
  stan_data <- list(
    N = nrow(df),
    H = ncol(B),
    B = B,
    status = df$status,
    time = df$time,
    J = J,
    time_knots = time_knots,
    dose_knots = seq(min(df$dose), max(df$dose), length.out = ncol(B))
  )
  
  bayes_fit <- sampling(stan_model("bayes_cox.stan"), data = stan_data,
                        chains = 3, iter = 5000, warmup = 2000, refresh = 0,
                        control = list(adapt_delta = 0.99, max_treedepth = 15))
  
  list(
    cox_linear = cox_linear,
    cox_spline = cox_spline,
    cox_pspline = cox_pspline,
    bayes = bayes_fit,
    B = B,
    knots = knots,
    time_knots = time_knots
  )
}
# Evaluate performance
evaluate_performance <- function(fit, df, scenario) {
  f_true <- switch(scenario,
                   linear = function(d) 0.1*d,
                   upward = function(d) 0.01 * d + 0.01 * d^2,
                   downward = function(d) 0.5 * d * exp(-0.2 * d),
                   sigmoid = function(d) 1 / (1 + exp(5 - d)) - 0.01
  )
  
  dose_grid <- seq(min(df$dose), max(df$dose), length.out = 100)
  B_grid <- bs(dose_grid, knots = fit$knots, degree = 2, intercept = TRUE)
  
  beta_post <- as.matrix(rstan::extract(fit$bayes, "beta")$beta)
  if (is.null(beta_post)) {
    return(data.frame(
      dose = dose_grid,
      truth = f_true(dose_grid),
      bias = NA,
      rmse = NA,
      coverage = NA,
      model = "Bayesian",
      scenario = scenario
    ))
  }
  
  f_hat <- B_grid %*% t(beta_post)
  
  truth <- f_true(dose_grid)
  truth <- truth - mean(truth)
  
  f_hat_centered <- sweep(f_hat, 2, colMeans(f_hat))  # subtract mean of each draw
  mean_f_hat <- rowMeans(f_hat_centered)
  bias <- mean_f_hat - truth
  sd <- apply(f_hat_centered, 1, sd)
  rmse <- sqrt(rowMeans((f_hat_centered - truth)^2))
  lower <- apply(f_hat_centered, 1, quantile, probs = 0.025)
  upper <- apply(f_hat_centered, 1, quantile, probs = 0.975)
  coverage <- (truth >= lower) & (truth <= upper)
  
  # Add frequentist model evaluations
  # Cox-linear predictions
  pred_linear <- predict(fit$cox_linear, newdata = data.frame(dose = dose_grid), type = "lp") # log-relative
  f_hat_linear <- pred_linear - mean(pred_linear)  # Center for comparison
  
  # Cox-spline predictions
  pred_spline <- predict(fit$cox_spline, newdata = data.frame(dose = dose_grid), type = "lp")
  f_hat_spline <- pred_spline - mean(pred_spline)
  
  # Cox-spline predictions
  pred_pspline <- predict(fit$cox_pspline, newdata = data.frame(dose = dose_grid), type = "lp")
  f_hat_pspline <- pred_pspline - mean(pred_pspline)
  
  # Combine metrics
  metrics <- bind_rows(
    # Bayesian metrics
    data.frame(dose = dose_grid, truth = truth, mean = mean_f_hat, bias = bias, rmse = rmse, 
               coverage = coverage, model = "Bayesian", scenario = scenario),
    # Cox-linear metrics
    data.frame(dose = dose_grid, truth = truth, mean = f_hat_linear, bias = f_hat_linear - truth,
               rmse = sqrt(mean((f_hat_linear - truth)^2)), 
               coverage = NA, model = "Cox-Linear", scenario = scenario),
    # Cox-spline metrics
    data.frame(dose = dose_grid, truth = truth, mean = f_hat_spline, bias = f_hat_spline - truth,
               rmse = sqrt(mean((f_hat_spline - truth)^2)), 
               coverage = NA, model = "Cox-Spline", scenario = scenario),
    # Cox-pspline metrics
    data.frame(dose = dose_grid, truth = truth, mean = f_hat_pspline, bias = f_hat_pspline - truth,
               rmse = sqrt(mean((f_hat_pspline - truth)^2)), 
               coverage = NA, model = "Cox-Pspline", scenario = scenario)
  )
  
  return(metrics)
}

# ---- Evaluation of f(d) ----
evaluate_Ftd <- function(fit, df, scenario, dose_levels = c(5, 10, 15), time_grid = NULL) {
  if (is.null(time_grid)) {
    time_grid <- seq(0.1, max(df$time), length.out = 100)
  }
  
  time_knots <- fit$time_knots
  lambda_draws <- rstan::extract(fit$bayes, "lambda")$lambda
  beta_draws <- rstan::extract(fit$bayes, "beta")$beta
  
  stopifnot(length(time_knots) == ncol(lambda_draws) + 1)
  
  # True F(t|d)
  true_Ftd <- function(d, t) {
    lp <- switch(scenario,
                 linear = 0.1 * d,
                 upward = 0.01 * d + 0.01 * d^2,
                 downward = 0.5 * d * exp(-0.2 * d),
                 sigmoid = 1 / (1 + exp(5 - d)) - 0.01)
    lambda <- 0.1 * exp(lp)
    1 - exp(-(lambda * t)^1.5)
  }
  
  # Bayesian F(t|d)
  bayes_Ftd <- function(d, t) {
    B_pred <- bs(d, knots = fit$knots, degree = 2, intercept = TRUE)
    valid_idx <- which(time_knots < t)
    j_t <- if (length(valid_idx) > 0) max(valid_idx) else NA_integer_
    if (is.na(j_t) || j_t < 1 || j_t > ncol(lambda_draws)) return(c(mean = NA, lower = NA, upper = NA))
    
    cumhaz_vals <- sapply(1:nrow(beta_draws), function(i) {
      linpred <- as.vector(B_pred %*% beta_draws[i, ])
      cumhaz <- sum(sapply(1:j_t, function(j) {
        t_start <- time_knots[j]
        t_end <- if (j == j_t) t else time_knots[j + 1]
        lambda_draws[i, j] * (t_end - t_start)
      }))
      exp(linpred) * cumhaz
    })
    
    cumhaz_vals <- cumhaz_vals[is.finite(cumhaz_vals)]
    if (length(cumhaz_vals) == 0) return(c(mean = NA, lower = NA, upper = NA))
    Ftd_samples <- 1 - exp(-cumhaz_vals)
    c(mean = mean(Ftd_samples), lower = quantile(Ftd_samples, 0.025), upper = quantile(Ftd_samples, 0.975))
  }
  
  # Frequentist F(t|d)
  get_Ftd_cox <- function(fit_model, d, t) {
    bh <- basehaz(fit_model, centered = FALSE)
    lp <- predict(fit_model, newdata = data.frame(dose = d), type = "lp")
    H0_t <- bh$hazard[which.max(bh$time >= t)]
    if (length(H0_t) == 0 || is.na(H0_t)) return(NA_real_)
    1 - exp(-H0_t * exp(lp))
  }
  
  expand_grid(dose = dose_levels, time = time_grid) %>%
    rowwise() %>%
    mutate(
      truth      = true_Ftd(dose, time),
      bayes      = list(bayes_Ftd(dose, time)),
      cox_linear = get_Ftd_cox(fit$cox_linear, dose, time),
      cox_spline = get_Ftd_cox(fit$cox_spline, dose, time),
      cox_pspline = get_Ftd_cox(fit$cox_pspline, dose, time),
      scenario   = scenario
    ) %>%
    unnest_wider(bayes)
}

# ---- Main Pipeline ----
results_nested <- set_names(scenarios) %>%
  map(function(scenario) {
    map(1:n_sim, function(sim) {
      cat("\nRunning", scenario, "- Simulation", sim, "/", n_sim)
      df <- generate_data(n, scenario)
      fits <- fit_models(df)
      list(
        fd = evaluate_performance(fits, df, scenario),
        ftd = evaluate_Ftd(fits, df, scenario)
      )
    }) %>% transpose()
  })


dose_response_results <- map_dfr(results_nested, ~ bind_rows(.x$fd), .id = "scenario")
Ftd_results <- map_dfr(results_nested, ~ bind_rows(.x$ftd), .id = "scenario")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~ RESULTS VISUALIZATION 1 ~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
performance_table <- dose_response_results %>%
  group_by(scenario, model) %>%
  summarise(
    avg_bias = mean(bias),
    avg_rmse = mean(rmse),
    coverage_rate = mean(coverage),
    .groups = "drop"
  )

# Plot: RMSE
ggplot(dose_response_results, aes(x = dose, y = rmse, color = model)) +
  geom_smooth(se = FALSE) +
  facet_wrap(~scenario) +
  labs(title = "RMSE Comparison Across Scenarios", x = "Dose", y = "RMSE") +
  theme_bw()

# Plot: mean
ggplot(dose_response_results, aes(x = dose, y = mean, color = model)) +
  geom_smooth(se = FALSE) +
  facet_wrap(~scenario) +
  labs(title = "Mean Comparison Across Scenarios", x = "Dose", y = "Mean") +
  theme_bw()

# Plot: bias
ggplot(dose_response_results, aes(x = dose, y = bias, color = model)) +
  geom_smooth(se = FALSE) +
  facet_wrap(~scenario) +
  labs(title = "Bias Comparison Across Scenarios", x = "Dose", y = "Bias") +
  theme_bw()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~ Comparison 2 ~~~~~~~
# ---- Visualization ----
Ftd_summary <- Ftd_results %>%
  filter(!is.na(mean)) %>%
  mutate(
    bias = mean - truth,
    squared_error = (mean - truth)^2,
    coverage = (truth >= `lower.2.5%`) & (truth <= `upper.97.5%`)
  ) %>%
  group_by(scenario, dose) %>%
  summarise(
    avg_bias = mean(bias),
    rmse = sqrt(mean(squared_error)),
    coverage_rate = mean(coverage),
    .groups = "drop"
  )

print(Ftd_summary)

Ftd_long <- Ftd_results %>%
  pivot_longer(cols = c(mean, truth, cox_linear, cox_spline, cox_pspline),
               names_to = "method", values_to = "Ftd") %>%
  mutate(method = recode(method,
                         mean = "Bayesian",
                         truth = "Truth",
                         cox_linear = "Cox-Linear",
                         cox_spline = "Cox-Spline",
                         cox_pspline = "Cox-Pspline"))

ggplot(Ftd_long, aes(x = time, y = Ftd, color = method, linetype = method)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = Ftd_results,
              aes(x = time, ymin = `lower.2.5%`, ymax = `upper.97.5%`),
              fill = "blue", alpha = 0.2, inherit.aes = FALSE) +
  facet_grid(rows = vars(scenario), cols = vars(dose), labeller = "label_value") +
  labs(title = "F(t | d) Curve Comparison", x = "Time", y = "F(t | d)",
       color = "Method", linetype = "Method") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_viridis_d()

