#function for mode
getmode <- function(v) {
  uniqv <- unique(na.omit(v)) #ZW added na.omit
  uniqv[which.max(tabulate(match(v, uniqv)))]
}