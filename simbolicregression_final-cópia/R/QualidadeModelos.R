RLQualidade <- function(Model, Y, Ye) {

  Y  <- as.matrix(Y)
  Ye <- as.matrix(Ye)

  Ym    <- matrix(rep(mean(Y[, 1]), nrow(Y)))
  Ymean <- cbind(Ym, 0)

  DMe     <- DistMallows(Model, Ye, Ymean)
  DM      <- DistMallows(Model, Y, Ymean)
  DMrmsem <- DistMallows(Model, Y, Ye)

  CD <- DMe / DM

  #transformação centro e raios em extremos
  Yeu <- Ye[, 1] + Ye[, 2]
  Yel <- Ye[, 1] - Ye[, 2]
  Yru <- Y[, 1] + Y[, 2]
  Yrl <- Y[, 1] - Y[, 2]

  A1 <- Yru - Yeu
  A2 <- Yrl - Yel

  B1 <- (mean(A1^2))
  B2 <- (mean(A2^2))

  RMSEM <- sqrt(DMrmsem / nrow(Ye))
  RMSEU <- B1^0.5
  RMSEL <- B2^0.5

  return(list(c("CD"= CD, "RMSEM" = RMSEM, "RMSEL" = RMSEL, "RMSEU" = RMSEU)))

}

#------------------------------------------------------------------------------------#
RLQualidade2 <- function(Y, Ye) {

  Y  <- as.matrix(Y)
  Cy <- Y[,1]
  Ry <- Y[,2]
  My <- Y[,3]

  Ye <- as.matrix(Ye)
  Cye <- Ye[,1]
  Rye <- Ye[,2]
  Mye <- Ye[,3]

  n  <- nrow(Ye)
  DistMallFinal <- 0

  for (i in 1:n) {
    DM <- DistMallows2(Cye[i],Rye[i],Mye[i],Cy[i],Ry[i],My[i])

    DistMallFinal <- DistMallFinal + DM
  }

  #transformação centro e raios em extremos
  Yeu <- Ye[, 1] + Ye[, 2]
  Yel <- Ye[, 1] - Ye[, 2]
  Yru <- Y[, 1] + Y[, 2]
  Yrl <- Y[, 1] - Y[, 2]

  A1 <- Yru - Yeu
  A2 <- Yrl - Yel

  B1 <- (mean(A1^2))
  B2 <- (mean(A2^2))

  RMSEM <- sqrt(DistMallFinal/ nrow(Ye))
  RMSEU <- B1^0.5
  RMSEL <- B2^0.5

  return(list(c("RMSEM" = RMSEM, "RMSEL" = RMSEL, "RMSEU" = RMSEU)))

}
