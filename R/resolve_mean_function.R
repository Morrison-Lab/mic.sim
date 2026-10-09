#' Accept the deprecated argument name `E[X|T,C]`
#'
#' simulate_mics() and simulation_run() used to call their component-mean
#' argument `E[X|T,C]`. That name can't be documented in Rd (Rd splits
#' argument names on commas), so it is now `mean_function`. Calls that still
#' pass `E[X|T,C]` by name land in `...`; this helper picks it up, warns, and
#' rejects any other unused arguments.
#'
#' @param mean_function The value of `mean_function` in the caller.
#' @param mean_function_missing Whether the caller's `mean_function` was missing.
#' @param dots `list(...)` from the caller.
#'
#' @return The component-mean function to use.
#' @keywords internal
#' @noRd
resolve_mean_function = function(mean_function, mean_function_missing, dots){
  old_name = "E[X|T,C]"
  if (old_name %in% names(dots)) {
    if (!mean_function_missing) {
      stop("Supply only one of `mean_function` and `E[X|T,C]`.", call. = FALSE)
    }
    warning("The argument `E[X|T,C]` is deprecated; use `mean_function` instead.",
            call. = FALSE)
    mean_function = dots[[old_name]]
    dots = dots[names(dots) != old_name]
  }
  if (length(dots) > 0) {
    unused = names(dots)
    if (is.null(unused)) unused = rep("", length(dots))
    unused[unused == ""] = "<unnamed>"
    stop("Unused argument(s): ", paste(unused, collapse = ", "), call. = FALSE)
  }
  mean_function
}
