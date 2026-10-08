#' Logistic Regression for Proportion Above an MIC Split Over Time
#'
#' Dichotomizes MICs at an ECOFF, a breakpoint, or a visual split and fits a
#' logistic generalized additive model (binomial family, smooth term s(t)) for
#' the proportion of isolates above the split over time. Used by plot_fm() to
#' compare the estimated component weights with these simpler summaries.
#'
#' @param data either possible_data or a set of mic data ready you want to run import_mics on
#' @param split_by String, which value to dichotomize the MICs at: "ecoff" (default), "visual_split", "s_breakpoint" (or "S"), or "r_breakpoint" (or "R"). The corresponding argument must be supplied.
#' @param data_type use either "possible_data" if passing in a possible_data object or "import" if you want to import mics and run logistic regression
#' @param drug NULL if using possible_data, if importing data should be the name of the column of the mics or a vector of the mic column and the sign column
#' @param date_col string, what is the name of the column in the data that corresponds to time of sampling
#' @param date_type string, either "decimal", "date", or "year" use decimal if using t from possible data, date or year if importing mic data and the date column is a date or just a year respectively
#' @param first_year NULL if date_type is "decimal", otherwise a numeric year or decimal year value if using "year" or "date" for date_type respectively
#' @param s_breakpoint string, the breakpoint on the MIC scale for what constitutes a susceptible isolate, e.g. <=8 (ug/mL, do not include units)
#' @param r_breakpoint string, the breakpoint on the MIC scale for what constitutes a resistant isolate, e.g. >=128 (ug/mL, do not include units)
#' @param ecoff string or numeric, see plot_fm()
#' @param visual_split string or numeric, see plot_fm()
#' @param k NULL (default) or a value passed as argument k to the mgcv::gam() call that fits the logistic regression
#' @return A fitted \code{gam} object (from mgcv) for the logistic regression
#'   of the dichotomized outcome (1 = resistant or non-wild-type, i.e. above the
#'   split; 0 = otherwise) on s(t).
#' @export
#'
#' @examples
#' \donttest{
#' data = simulate_mics()
#' output = fit_EM(model = "pspline", approach = "full", pre_set_degrees = c(4, 4),
#'                 visible_data = data, max_it = 300, verbose = 0)
#' log_reg(output$possible_data, split_by = "ecoff", data_type = "possible_data",
#'         drug = NULL, date_col = "t", date_type = "decimal", first_year = NULL,
#'         ecoff = 2)
#' }
log_reg <- function(data, split_by = "ecoff", data_type, drug, date_col, date_type, first_year, s_breakpoint, r_breakpoint, ecoff, visual_split, k = NULL){
  #assume ecoff is by default "an MIC that is less than or equal to ecoff is WT"
  if(date_type == "decimal"){
    df_temp = data %>% rename(t = date_col)
  }else if(date_type == "date"){
    df_temp = data %>% rename(date = date_col) %>%
      mutate(t = decimal_date(date) - first_year) %>%
      suppressWarnings()
  } else if(date_type == "year"){
    df_temp = data %>% rename(date = date_col) %>%
      mutate(t = as.numeric(date) - first_year) %>%
      suppressWarnings()
  }else{
    stop("pick decimal or year", call. = FALSE)
  }

  if(split_by == "S" | split_by == "s_breakpoint"){
    if(!is.null(s_breakpoint)){
      divider = case_when(grepl(pattern = "(\u2264)|(<=)|(=<)", x = s_breakpoint) ~ parse_number(as.character(s_breakpoint)),
                          grepl(pattern = "(<)", x = s_breakpoint) & !grepl(pattern = "(\u2264)|(<=)|(=<)", x = s_breakpoint) ~ parse_number(as.character(s_breakpoint)) - 0.00001,
                          TRUE ~ parse_number(as.character(s_breakpoint))
      )
    }else{
      stop("split_by set to s_breakpoint and no s_breakpoint provided", call. = FALSE)
    }

  }else if(split_by == "R" | split_by == "r_breakpoint"){
    if(!is.null(r_breakpoint)){
      divider = case_when(grepl(pattern = "(\u2265)|(>=)|(=>)", x = r_breakpoint) ~ parse_number(as.character(r_breakpoint)),
                          grepl(pattern = "(>)", x = r_breakpoint) & !grepl(pattern = "(\u2265)|(>=)|(=>)", x = r_breakpoint) ~ parse_number(as.character(r_breakpoint)) + 0.00001,
                          TRUE ~ parse_number(as.character(r_breakpoint))
      )
    }else{
      stop("split_by set to r_breakpoint and no r_breakpoint provided", call. = FALSE)

    }
  }else if(split_by == "visual_split"){

    if(!is.null(visual_split)){
      divider = case_when(#grepl(pattern = "(\u2265)|(>=)|(=>)", x = visual_split) ~ parse_number(as.character(visual_split)),
        grepl(pattern = "(<)", x =visual_split) & !grepl(pattern = "(\u2264)|(<=)|(=<)", x = visual_split) ~ parse_number(as.character(visual_split)) - 0.00001,
        TRUE ~ parse_number(as.character(visual_split))
      )
    }else{
      stop("split_by set to visual_split and no visual_split provided, either change split_by to s_breakpoint or r_breakpoint or ecoff or else provide a value for visual_split", call. = FALSE)

    }
  }else{

    if(!is.null(ecoff)){
      divider = case_when(#grepl(pattern = "(\u2265)|(>=)|(=>)", x = ecoff) ~ parse_number(as.character(ecoff)),
        grepl(pattern = "(<)", x =ecoff) & !grepl(pattern = "(\u2264)|(<=)|(=<)", x = ecoff) ~ parse_number(as.character(ecoff)) - 0.00001,
        TRUE ~ parse_number(as.character(ecoff))
      )
    }else{
      stop("split_by set to ecoff and no ecoff provided, either change split_by to s_breakpoint or r_breakpoint or visual_split or else provide a value for ecoff", call. = FALSE)

    }

  }





  if(data_type == "import"){
    df = import_mics((data %>% select(all_of(drug))) %>%
                       pull(drug)) %>%
      tibble(., df_temp) %>%
      filter(!is.na(mic_column)) %>%
      mutate(
        cens = case_when(
          left_bound == 0 ~ "LC",
          right_bound == Inf ~ "RC",
          TRUE ~ "int"
        ),
        mic = case_when(
          cens == "int" ~ parse_number(as.character(mic_column)),
          cens == "LC" ~ parse_number(as.character(mic_column)),
          cens == "RC" ~ parse_number(as.character(mic_column)) * 2,
          TRUE ~ NaN
        )
      )

  } else if(data_type == "possible_data"){
    log_divider = ceiling(log2(divider))
    df = df_temp %>% filter(c == 1) %>%
      mutate(
        cens = case_when(
          left_bound == -Inf ~ "LC",
          right_bound == Inf ~ "RC",
          TRUE ~ "int"
        ),
        mic = case_when(
          cens == "int" ~ right_bound,
          cens == "LC" ~ right_bound, #previous - 0.01
          cens == "RC" ~ left_bound + 1, #previous + 0.01
          TRUE ~ NaN
        )
      )


  }else{
    stop("choose either import or possible_data for data_type", call. = FALSE)
  }
  if(split_by == "s_breakpoint" | split_by == "S"){
    split_df = df %>%
      mutate(
        sir = case_when(
          mic <= log_divider ~ "S",
          TRUE ~ "R"
        )
      ) %>%
      mutate(dichot_res = case_when(
        sir %in% c("R") ~ 1,
        TRUE ~ 0
      ))
    if(is.null(k)){
      split_df %>%
      mgcv::gam(formula = dichot_res ~ s(t), method = "REML", family = binomial(link = "logit")) %>% return()
    }else{
      split_df %>%
      mgcv::gam(formula = dichot_res ~ s(t), k = k, method = "REML", family = binomial(link = "logit")) %>% return()
      }

  }else if(split_by == "r_breakpoint" | split_by == "R"){
    split_df = df %>%
      mutate(
        sir = case_when(
          mic >= log_divider ~ "R",
          TRUE ~ "S"
        )
      ) %>%
      mutate(dichot_res = case_when(
        sir %in% c("R") ~ 1,
        TRUE ~ 0
      ))

    if(is.null(k)){
      split_df %>%
        mgcv::gam(formula = dichot_res ~ s(t), method = "REML", family = binomial(link = "logit")) %>% return()
    }else{
      split_df %>%
        mgcv::gam(formula = dichot_res ~ s(t), k = k, method = "REML", family = binomial(link = "logit")) %>% return()
    }
  }else{
    #currently handling visual_split and ecoff the same
    split_df = df %>%
      mutate(
        wtnwt = case_when(
          mic <= log_divider ~ "WT",
          TRUE ~ "NWT"
        )
      ) %>%
      mutate(dichot_res = case_when(
        wtnwt == "NWT" ~ 1,
        TRUE ~ 0
      ))
      if(is.null(k)){
        split_df %>%
          mgcv::gam(formula = dichot_res ~ s(t), method = "REML", family = binomial(link = "logit")) %>% return()
      }else{
        split_df %>%
          mgcv::gam(formula = dichot_res ~ s(t), k = k, method = "REML", family = binomial(link = "logit")) %>% return()
      }
  }


}
