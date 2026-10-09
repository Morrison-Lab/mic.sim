sim_two_component_data = function(n = 300, seed = 1, scale = "log"){
  set.seed(seed)
  simulate_mics(
    n = n,
    t_dist = function(n){runif(n, min = 0, max = 10)},
    pi = function(t) {
      z <- 0.3 + 0.02 * t
      tibble::tibble("1" = 1 - z, "2" = z)
    },
    mean_function = function(t, c){
      dplyr::case_when(c == "1" ~ -2 + 0.1 * t,
                       c == "2" ~ 3,
                       TRUE ~ NaN)
    },
    sd_vector = c("1" = 1, "2" = 1),
    low_con = if (scale == "MIC") 2^-5 else -5,
    high_con = if (scale == "MIC") 2^6 else 6,
    scale = scale
  ) %>% suppressMessages()
}
