#' Prepare Simulated Data for the EM Algorithm
#'
#' Intermediate function to prepare simulated data for use in EM algorithm:
#' selects the interval bounds, time, covariates, and tested concentration
#' limits from the output of simulate_mics(), converts the bounds to the log2
#' scale if needed, and adds an observation id.
#'
#' @param data.sim Data frame of simulated data, usually the output of simulate_mics()
#' @param left_bound_name String, name of the column in data.sim containing the left bound of each censoring interval
#' @param right_bound_name String, name of the column in data.sim containing the right bound of each censoring interval
#' @param time String, name of the column in data.sim containing the time variable
#' @param covariate_names Character vector or NULL, names of covariate columns in data.sim to keep
#' @param scale String or NULL, scale of the bounds in data.sim: "MIC" (bounds are converted to log2) or "log" (bounds are kept as is). If NULL, the "scale" attribute of data.sim is used.
#' @param keep_truth Logical, if TRUE the true (uncensored) simulated value and true component of each observation are kept in columns observed_value and comp
#' @param observed_value_name String, name of the column in data.sim containing the true simulated value, used if keep_truth is TRUE
#' @param comp_name String, name of the column in data.sim containing the true component, used if keep_truth is TRUE
#' @param low_con_name String, name of the column in data.sim containing the lowest tested concentration
#' @param high_con_name String, name of the column in data.sim containing the highest tested concentration
#'
#' @return A tibble with columns obs_id, the covariate_names columns, the time
#'   column, left_bound and right_bound (on the log2 scale), low_con, and
#'   high_con, plus observed_value and comp if keep_truth is TRUE.
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select all_of mutate n relocate
#'
#' @examples
#' prep_sim_data_for_em(simulate_mics(n = 50))
prep_sim_data_for_em <- function(
    data.sim = simulate_mics(),
    left_bound_name = "left_bound",
    right_bound_name = "right_bound",
    time = "t",
    covariate_names = NULL,
    scale = NULL,
    keep_truth = FALSE,
    observed_value_name = "observed_value",
    comp_name = "comp",
    low_con_name = "low_con",
    high_con_name = "high_con"
) {

if(keep_truth){
truth <- data.sim %>%
  rename(observed_value = match(paste0(observed_value_name), names(data.sim))) %>%
  rename(comp = match(paste0(comp_name), names(data.sim))) %>%
  select("observed_value", "comp")

}

if(is.null(scale) && attr(data.sim, "scale") == "MIC")  {
    df <- data.sim %>%
      select(all_of(c(covariate_names, time)), left_bound = all_of(left_bound_name), right_bound = all_of(right_bound_name), low_con = all_of(low_con_name), high_con = all_of(high_con_name)) %>%
      mutate(obs_id = 1:n(),
             left_bound = log2(left_bound),
             right_bound = log2(right_bound)) %>%
      relocate(obs_id, .before = everything())
    if(keep_truth){
      df <- cbind(df, truth) %>% tibble()
    }
  }
  else if (is.null(scale) && attr(data.sim, "scale") == "log"){
    df <- data.sim %>%
      select(all_of(c(covariate_names, time)), left_bound = all_of(left_bound_name), right_bound = all_of(right_bound_name), low_con = all_of(low_con_name), high_con = all_of(high_con_name)) %>%
      mutate(obs_id = 1:n()) %>%
      relocate(obs_id, .before = everything())
    if(keep_truth){
      df <- cbind(df, truth) %>% tibble()
    }
  }
  else if(scale == "MIC"){
    df <- data.sim %>%
    select(all_of(c(covariate_names, time)), left_bound = all_of(left_bound_name), right_bound = all_of(right_bound_name), low_con = all_of(low_con_name), high_con = all_of(high_con_name)) %>%
    mutate(obs_id = 1:n(),
           left_bound = log2(left_bound),
           right_bound = log2(right_bound)) %>%
    relocate(obs_id, .before = everything())
    if(keep_truth){
      df <- cbind(df, truth) %>% tibble()
    }
}
else if (scale == "log"){
  df <- data.sim %>%
      select(all_of(c(covariate_names, time)), left_bound = all_of(left_bound_name), right_bound = all_of(right_bound_name), low_con = all_of(low_con_name), high_con = all_of(high_con_name)) %>%
      mutate(obs_id = 1:n()) %>%
      relocate(obs_id, .before = everything())
  if(keep_truth){
    df <- cbind(df, truth) %>% tibble()
  }
  }
  else{warning("Either (A.) set scale variable to MIC or log or (B.) data.sim should have a scale attribute of either MIC or log", call. = FALSE)}

return(df)

  }
