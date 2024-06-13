RLParametrosModelo1 <- function(Y, X) {

  #  Modelo1 -- Distribuição Uniforme
  #  Y -- Variável Resposta
  #         Y -> yc,yr
  #  X -- p Variáveis Explicativas
  #         X -> xc1,xr1,xc2,xr2,...,xcp,xrp

  ################################################################################################################################

  Y <- as.matrix(Y)
  X <- as.matrix(X)

  Yc <- as.matrix(Y[, 1])                                      #matriz com os centros da variável resposta
  Yr <- as.matrix(Y[, 2])                                      #matriz com os raios da variável resposta
  Xc <- as.matrix(X[, seq(from = 1, to = ncol(X), by = 2)])    #matriz com os centros das variaveis explicativas
  Xr <- as.matrix(X[, seq(from = 2, to = ncol(X), by = 2)])    #matriz com os raios das variáveis explicativas
  p  <- ncol(Xc)                                               #número de variáveis

  ################################################################################################################################

  NHhessianaBlocoI <- function(Y,X) {

    H   <- matrix(0, p * 2, p * 2)
    Xcc <- crossprod(Xc)
    Xrr <- crossprod(Xr)

    for(k in 1:(p * 2)) {
      for(l in 1:(p * 2)) {
        if(k %% 2 != 0) {
          H <- H
        } else if(l %% 2 != 0) {
          H <- H
        } else {
          H[k, l]     <-  2 * Xcc[k/2, l/2] + 2/3 * Xrr[k/2, l/2]
          H[k-1, l-1] <-  2 * Xcc[k/2, l/2] + 2/3 * Xrr[k/2, l/2]
          H[k-1, l]   <-  -2 * Xcc[k/2, l/2] + 2/3 * Xrr[k/2, l/2]
          H[k, l-1]   <-  -2 * Xcc[k/2, l/2] + 2/3 * Xrr[k/2, l/2]
        }
      }
    }

    H <- round(H, 2)

    return(H)
  }

  H1 <- tryCatch({
    NHhessianaBlocoI(Y, X)
  }, error = function(e) {
    return(NA)
  })

  if(is.na(H1[1]))
    return(NA)


  ################################################################################################################################

  NHhessianaBlocoLinha <- function(Y, X) {

    HA <- matrix(0, 1, p * 2)

    for(l in 1:(p * 2)) {
      if(l %% 2 != 0) {
        HA <- HA
      } else {
        HA[, l]   <- -2 * sum(Xc[, l/2])
        HA[, l-1] <- 2 * sum(Xc[, l/2])
      }
    }

    HA <- list("down" = HA, "side" = t(HA), "constante" = nrow(Xc) * 2)

    return(HA)
  }

  HA <- tryCatch({
    NHhessianaBlocoLinha(Y, X)
  }, error = function(e) {
    return(NA)
  })

  if(is.na(HA[1]))
    return(NA)

  ################################################################################################################################

  NHhessiana <- function(Y, X) {

    NH <- cbind(rbind(H1, HA$down), rbind(HA$side, HA$constante))
    return(NH)
  }

  #H <- NHhessiana(Y, X)
  H <- tryCatch({
    NHhessiana(Y, X)
  }, error = function(e) {
    return(NA)
  })

  if(is.na(H[1]))
    return(NA)

  ################################################################################################################################

  NHfuncao <- function(Y, X) {

    Hf  <- matrix(0, p * 2, 1)                             #matriz hessianafuncaoI preenchida com zeros
    Cf  <- matrix(0, 1, 1)
    YXc <- matrix(0, p, 1)
    YXr <- matrix(0, p, 1)

    for(i in 1:p) {
      YXc[i, 1] <- crossprod(Yc[, 1], Xc[, i])
      YXr[i, 1] <- crossprod(Yr[, 1], Xr[, i])
    }

    for(l in 1:p * 2) {
      if(l %% 2 != 0) {
        Hf <- Hf
      } else {
        Hf[l, ]   <- 2 * YXc[l/2, ] - 2/3 * YXr[l/2, ]
        Hf[l-1, ] <- -2 * YXc[l/2, ] - 2/3 * YXr[l/2, ]
      }
    }

    Cf[1, 1] <- -2 * sum(Yc)
    w1 <- rbind(Hf, Cf)

    return(w1)
  }

  #F <-  NHfuncao(Y, X) * (-1)
  F <- tryCatch({
    NHfuncao(Y, X) * (-1)
  },
  error = function(e) {
    return(NA)
  })

  if(is.na(F[1]))
    return(NA)
  ################################################################################################################################

  library("quadprog")

  A  <- matrix(0, p*2, p*2)

  diag(A) <- 1

  A  <- rbind(A, matrix(0, 1, p*2))

  b  <- matrix(0, p*2, 1)

  #RL <- solve.QP(H, F, A, b)
  RL <- tryCatch({
    solve.QP(H, F, A, b)
  },
  error = function(e) {
    return(NA)
    #stop("Error in solve.QP(H, F, A, b) : matrix D in quadratic function is not positive definite! Try reducing number of clusters/iterations.")
  })

  if(is.na(RL[1]))
    return(NA)

  return(RL$solution)

}
