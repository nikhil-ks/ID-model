SDM <- function(x,media) {

  sum <- 0
  n <- nrow(x)

  for(i in 1:n) {
    DM <- DistMallows2(x[i,1],x[i,2],x[i,1],media,0,media)
    sum <- sum + DM
  }
  return(sum)
}

SDM1 <- function(x,media) {

  sum <- 0
  n <- nrow(x)

  for(i in 1:n) {
    DM <- DistMallows2(x[i,1],x[i,2],x[i,1],media[i,1],media[i,2],media[i,1])
    sum <- sum + DM
  }
  return(sum)
}
