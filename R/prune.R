#' @title Prune an empirical cluster tree
#' @name prune
#' @description More details coming soon...
#' @param eps A tuning parameter (> 0) representing how aggressively to prune.
#' @references See KC and SD (2014).
#' @export
prune <- function(C_n, eps, delta){
  n <- length(C_n$hc$height) + 1
  d <- ifelse(is.null(C_n[["d"]]), 1, C_n[["d"]])
  k <- C_n$k
  beta_n <- sqrt((4/n) * ((d * log(2*n)) + log(8 / delta)))
  tmp <- d*log(n) + log(1/delta)
  C_o <- n * ((beta_n^2 + beta_n * sqrt(C_n$k/n))/(tmp + sqrt(C_n$k * tmp)))
  C_del <- 2 * C_o * log(2/delta)

  # If the significance level is too high, then all attempts at pruning
  # C_del has to be <= k - sqrt(k * d * log(n))
  const <- (C_del/n) * sqrt(k * d * log(n))
  if ((k / n) - const < 0){
    message <- "Pruning did not result any changes to the cluster tree. Try lowering the universal constant.\n"
    message <- paste0(message, sprintf("Cdelta = %f, eps = %f, delta = %f"))
    warning(message)
    return(C_n)
  }

  # Given a lambda, generate the appropriate r
  r_ <- function(lambda) { ((k/n + const)/(clustertree:::vol_nSphere(d)*lambda))^(1/d) }

  # Upper bound of the universal constant
  prune_lvls <- sapply(C_n$hc$height, function(r){ (1/clustertree:::vol_nSphere(d) * r^(1/d)) * comp })

  # Get the lambda values for which to
  valid_r <- which((prune_lvls - eps) > 0)
  r_prime <- r_(prune_lvls - eps)[valid_r]
  r <- C_n$hc$height[valid_r]

}





# plot(get_del, from = 0, to = 1)



