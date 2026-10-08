fake_surv_fit = function(coefficients, scale){
  structure(list(coefficients = coefficients, scale = scale), class = "survreg")
}

test_that("pi convergence check fails when coefficients decrease", {
  mu = list(fake_surv_fit(c(1, 2), 1), fake_surv_fit(c(3, 4), 1))
  pi_old = list(coefficients = c(0, 0))
  pi_new = list(coefficients = c(-5, -5))

  expect_false(
    model_coefficient_checks.surv(mu, pi_new, mu, pi_old, 1e-5, ncomp = 2)
  )
  expect_true(
    model_coefficient_checks.surv(mu, pi_old, mu, pi_old, 1e-5, ncomp = 2)
  )
})

test_that("get_scale returns sigma for survreg and mgcv fits", {
  set.seed(1)
  df = data.frame(t = runif(200, 0, 10))
  df$y = 1 + 0.2 * df$t + rnorm(200, sd = 2)
  df$left = df$y
  df$right = df$y

  surv_fit = survival::survreg(
    survival::Surv(left, right, type = "interval2") ~ t,
    data = df, dist = "gaussian")
  gam_fit = mgcv::gam(cbind(left, right) ~ t, family = mgcv::cnorm(),
                      data = df)

  expect_equal(get_scale(surv_fit), surv_fit$scale)
  expect_equal(get_scale(gam_fit), gam_fit$family$getTheta(TRUE))
  expect_equal(get_scale(gam_fit), get_scale(surv_fit), tolerance = 0.05)
})

test_that("set_scale_log puts bounds and tested range on the log2 scale", {
  log_data = sim_two_component_data(scale = "log")
  mic_data = sim_two_component_data(scale = "MIC")

  converted = set_scale_log(mic_data, NULL)

  expect_equal(attr(converted, "scale"), "log")
  expect_equal(converted$left_bound, log_data$left_bound)
  expect_equal(converted$right_bound, log_data$right_bound)
  expect_equal(converted$low_con, log_data$low_con)
  expect_equal(converted$high_con, log_data$high_con)

  # already-log data is left alone, and converting twice is a no-op
  expect_equal(set_scale_log(converted, NULL)$left_bound, converted$left_bound)
})

test_that("add_scale normalizes scale names and rejects bad ones", {
  df = tibble::tibble(left_bound = 1, right_bound = 2)

  expect_equal(attr(add_scale(df, "fold"), "scale"), "log")
  expect_equal(attr(add_scale(df, "concentration"), "scale"), "MIC")
  expect_error(add_scale(df, NULL), "scale attribute")
  expect_error(add_scale(df, "linear"), "Invalid value of scale")
})

test_that("invalid arguments raise errors instead of passing silently", {
  expect_error(
    fit_all_mu_models(tibble::tibble(), ncomp = 2, mu_formula = list(),
                      approach = "nonsense", fixed_side = NULL,
                      maxiter_survreg = 30)
  )
})

test_that("maxiter_survreg reaches survreg", {
  possible_data = sim_two_component_data() %>%
    dplyr::mutate(c = "1", `P(C=c|y,t)` = 1)
  formula = survival::Surv(time = left_bound, time2 = right_bound,
                           type = "interval2") ~ t

  fit = fit_mu_model_safe_formatted(possible_data, "1", formula,
                                    maxiter_survreg = 1)

  expect_lte(fit$iter[1], 1)
})

test_that("rsev draws from the standard smallest extreme value distribution", {
  set.seed(1)
  x = rsev(1e5)
  expect_equal(mean(x), digamma(1), tolerance = 0.02)
  expect_equal(var(x), pi^2 / 6, tolerance = 0.02)
})

test_that("EM_algorithm gives the same fit for MIC-scale and log-scale input", {
  fit_log = EM_algorithm(sim_two_component_data(scale = "log"),
                         mu_formula = survival::Surv(time = left_bound,
                                                     time2 = right_bound,
                                                     type = "interval2") ~ t,
                         max_it = 300, verbose = 0, scale = "log")
  fit_mic = EM_algorithm(sim_two_component_data(scale = "MIC"),
                         mu_formula = survival::Surv(time = left_bound,
                                                     time2 = right_bound,
                                                     type = "interval2") ~ t,
                         max_it = 300, verbose = 0, scale = "MIC")

  expect_equal(fit_log$converge, "YES")
  expect_equal(fit_mic$converge, "YES")
  expect_equal(coef(fit_mic$mu_model[[1]]), coef(fit_log$mu_model[[1]]),
               tolerance = 1e-6)
  expect_equal(coef(fit_mic$mu_model[[2]]), coef(fit_log$mu_model[[2]]),
               tolerance = 1e-6)
})

test_that("the default logit pi model is fitted, not rejected", {
  possible_data = tibble::tibble(
    t = rep(seq(0, 10, length.out = 50), each = 2),
    c = rep(c("1", "2"), 50),
    `P(C=c|y,t)` = rep(c(0.7, 0.3), 50)
  )
  expect_s3_class(
    fit_mgcv_pi_model(c == "2" ~ s(t), "logit", possible_data), "gam")
  expect_s3_class(
    fit_mgcv_pi_model(c == "2" ~ t, "logit_simple", possible_data), "glm")
  expect_error(fit_mgcv_pi_model(c == "2" ~ t, "probit", possible_data),
               "link function")
})
