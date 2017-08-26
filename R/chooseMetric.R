## chooseMetric.R
## Function to, given a mertic name, check against the available metrics in the package to use and instantiate the
## appropriate metric, if available (stop otherwise). Returns a checked external pointer to the Metric objects memory.
chooseMetric <- function(metric_name, config = list(), ...){
  ## Choose the metric
  metric_idx <- pmatch(metric_name, .supported_metrics)
  if (is.na(metric_idx)) {
    available_metrics <- paste0(.supported_metrics, collapse = ", ")
    stop(paste0("Unknown metric specified. Please use one of: < ", available_metrics, " >"))
  }
  metric_par <- append(config, list(...))
  metric_ptr <- chooseMetric_int(metric_name, metric_par);
  if (class(metric_ptr) == "externalptr" && deparse(metric_ptr) != "<pointer: 0x0>"){
    return(metric_ptr);
  } else { stop("Unable to use the supplied metric with the parameters given.") }
}