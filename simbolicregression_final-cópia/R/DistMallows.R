DistMallows <- function(Model, A, B) {

  Centros <- cbind(as.matrix(A[,1]), as.matrix(B[,1]))
  Raios   <- cbind(as.matrix(A[,2]), as.matrix(B[,2]))

  # C <- matrix(rep(0, x = nrow(Centros)))
  C <- matrix(rep(0, nrow(Centros)))
  R <- matrix(rep(0, nrow(Raios)))

  for(i in 1:nrow(Centros)) {
    C[i, ] <- (Centros[i, 1] - Centros[i, 2])^2
  }

  for(i in 1:nrow(Raios)) {
    if(Model == 1) {
      R[i, ] <- (1 / 3) * (Raios[i, 1] - Raios[i, 2])^2
    } else if (Model == 2) {
      R[i, ] <- (1 / 6) * (Raios[i, 1] - Raios[i, 2])^2
    } else {
      print("escolher modelo 1-Uniforme, 2-Triangular Simetrico, 3-Triangular Geral")
    }
  }

  M  <- C + R
  DM <- sum(M[, 1])
  return(DM)
}

