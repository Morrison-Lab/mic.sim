mean_fn = function(t, c){
  dplyr::case_when(c == "1" ~ -2 + 0.1 * t, c == "2" ~ 3, TRUE ~ NaN)
}

test_that("simulate_mics accepts mean_function and the deprecated E[X|T,C]", {
  set.seed(1)
  new = simulate_mics(n = 50, mean_function = mean_fn) %>% suppressMessages()
  set.seed(1)
  expect_warning(
    old <- simulate_mics(n = 50, `E[X|T,C]` = mean_fn) %>% suppressMessages(),
    "deprecated")
  expect_equal(old, new)
})

test_that("simulate_mics rejects both names at once and unknown arguments", {
  expect_error(simulate_mics(n = 5, mean_function = mean_fn, `E[X|T,C]` = mean_fn),
               "only one")
  expect_error(simulate_mics(n = 5, not_an_argument = 1), "Unused argument")
})
