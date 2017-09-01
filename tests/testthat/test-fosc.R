library("clustertree")
sl <- clustertree::clustertree(iris[, 1:4], k = 5L)
sl_si <- clustertree:::simplified_hclust(sl$hc, 2)

sf <- function(cl) { if(length(cl) == 0) return(0) else return(max(cl)) }
what <- clustertree:::fosc(sl_si, sf)