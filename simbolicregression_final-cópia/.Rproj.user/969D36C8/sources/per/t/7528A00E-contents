getDistance <- function(Y, Ye = NULL) {
  if (is.null(Ye)) {
    Y  <- as.matrix(Y)
    Ymean <- cbind(matrix(rep(mean(Y[, 1]), nrow(Y))), 0)
    DM <- DistMallows(1, Y, Ymean)
    return(DM)
  } else {
    Y  <- as.matrix(Y)
    Ye <- as.matrix(Ye)
    DM <- DistMallows(1, Y, Ye)
    return(DM)
  }
}
