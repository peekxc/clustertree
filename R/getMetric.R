getMetric_int <- function(metric_name, ...){
  config <- as.list(...)
  # TODO: Test for the parameter settings
  metric_type <- chooseMetric(metric_name, config);
  return(metric_type)
}