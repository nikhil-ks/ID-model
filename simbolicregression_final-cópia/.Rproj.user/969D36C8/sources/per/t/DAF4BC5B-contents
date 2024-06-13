ToCentreRange <- function(X) {

  simbolic <- matrix(data = 0, nrow = nrow(X), ncol = 2)
  simbolic[,1] <- (X[,1] + X[,2])/2
  simbolic[,2] <- (X[,2] - X[,1])/2
  return(simbolic)
}

TOCENTRERANGE <- function(X) {

  var <- ncol(X)-1
  A <- NULL

  for(i in seq(from = 1, to = var, by = 2)) {
    D <- ToCentreRange(X[,c(i,i+1)])

    A <- cbind(A,D)
  }
  return(A)
}
