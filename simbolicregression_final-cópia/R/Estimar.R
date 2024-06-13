RLEstimar <- function(Xe, RL) {
  #  Modelo 1 e 2
  #  Xe - Variaveis explicativas para estimar
  #  RL - parametros resultantes dos modelos: a1, b1, ... ap,bp,const

  Xe <- as.matrix(Xe)

  parametros <- matrix(RL[-length(RL)], (length(RL) - 1) / 2, 2, byrow = TRUE)

  parametros_raios <- as.matrix(parametros[, 1] + parametros[, 2])

  parametros_centros <- rbind(as.matrix(parametros[, 1] - parametros[, 2]), RL[length(RL)])

  ###############################################################################################################

  Xc <- as.matrix(Xe[, seq(from = 1, to = ncol(Xe), by = 2)])    #matriz com os centros das variaveis explicativas
  Xr <- as.matrix(Xe[, seq(from = 2, to = ncol(Xe), by = 2)])    #matriz com os raios das variáveis explicativas

  ###############################################################################################################

  Yec <- cbind(Xc, 1) %*% parametros_centros
  Yer <- cbind(Xr)    %*% parametros_raios

  return(cbind(Yec, Yer))
}

#------------------------------------------------------------------------------------------------------------------------------------

Modelo3Previsao <- function(Xe,RL) {
  #  Modelo 3
  #  Xe - Variaveis explicativas para estimar
  #  RL - parametros resultantes dos modelos: alfa, beta e gama

  alfa <- RL[1]
  beta <- RL[2]
  gama <- RL[3]

  X  <- as.matrix(Xe)
  n  <- nrow(X)                                                #número de observações
  Ye <- matrix(0, n, 3)

  XC <- as.matrix(X[, 1])    #matriz com os centros das variaveis explicativas
  XR <- as.matrix(X[, 2])    #matriz com os raios das variáveis explicativas
  XM <- as.matrix(X[, 3])    #matriz com as modas das variáveis explicativas

  for(i in 1:n) {

    Xc <- 0
    Xr <- 0
    Xm <- 0

    Xc <- XC[i, ]
    Xr <- XR[i, ]
    Xm <- XM[i, ]
    Xa <- Xc - Xr
    Xb <- Xc + Xr

    Ymest <- 0
    Yaest <- 0
    Ybest <- 0
    Yrest <- 0
    Ycest <- 0


    if(round((Xm - Xa),4) > round((Xb - Xm),4)) {
      #t <- (alfa^(2 / 3)) / ((alfa^(2 / 3)) + beta^(2 / 3))
      if((alfa > ((beta * (Xb - Xm)^(3 / 2)) / (Xm - Xa)^(3 / 2))) & (alfa < ((beta * (Xm - Xa)^(3 / 2)) / (Xb - Xm)^(3 / 2)))) {
        Ye[i,3] <- (alfa-beta)*Xa + sqrt((Xb-Xa)*(Xm-Xa))*((alfa^(4/3)-beta^(4/3))/sqrt(alfa^(2/3)+beta^(2/3))) + gama
        Yaest   <- alfa * Xa - beta * Xb + gama
        Ybest   <- alfa * Xb - beta * Xa + gama
        Ye[i,1] <- (Ybest + Yaest)/2
        Ye[i,2] <- (Ybest - Yaest)/2
      } else if(alfa < ((beta * (Xb - Xm) ^ (3 / 2)) / (Xm - Xa)^(3 / 2))) {
        Ye[i,3] <- alfa * Xa - beta * Xm + alfa * sqrt((Xb-Xm) * (Xm-Xa)) + gama
        Yaest   <- alfa * Xa - beta * Xb + gama
        Ybest   <- alfa * Xb - beta * Xa + gama
        Ye[i,1] <- (Ybest + Yaest)/2
        Ye[i,2] <- (Ybest - Yaest)/2
      } else if(alfa > ((beta * (Xm - Xa)^(3 / 2)) / (Xb - Xm)^(3 / 2))) {
        Ye[i,3] <- alfa * Xm - beta * Xa + beta * sqrt((Xb-Xm) * (Xm-Xa)) + gama
        Yaest   <- alfa * Xa - beta * Xb + gama
        Ybest   <- alfa * Xb - beta * Xa + gama
        Ye[i,1] <- (Ybest + Yaest)/2
        Ye[i,2] <- (Ybest - Yaest)/2
      }

    } else if(round((Xm - Xa),4) < round((Xb - Xm),4)) {
     # t <- (beta^(2 / 3)) / ((alfa^(2 / 3)) + beta^(2 / 3))
      if((alfa > ((beta * (Xm - Xa)^(3 / 2)) / (Xb - Xm)^(3 / 2))) & (alfa < ((beta * (Xb - Xm)^(3 / 2)) / (Xm - Xa)^(3 / 2)))) {
        Ye[i,3] <- (alfa-beta)*Xb + sqrt((Xb-Xa)*(Xb-Xm))*((beta^(4/3)-alfa^(4/3))/sqrt(alfa^(2/3)+beta^(2/3))) + gama
        Yaest   <- alfa * Xa - beta * Xb + gama
        Ybest   <- alfa * Xb - beta * Xa + gama
        Ye[i,1] <- (Ybest + Yaest)/2
        Ye[i,2] <- (Ybest - Yaest)/2
      } else if(alfa < ((beta * (Xm - Xa)^(3 / 2)) / (Xb - Xm)^(3 / 2))) {
        Ye[i,3] <- alfa * Xb - beta * Xm + alfa * sqrt((Xb-Xm) * (Xm-Xa)) + gama
        Yaest   <- alfa * Xa - beta * Xb + gama
        Ybest   <- alfa * Xb - beta * Xa + gama
        Ye[i,1] <- (Ybest + Yaest)/2
        Ye[i,2] <- (Ybest - Yaest)/2
      } else if(alfa > ((beta * (Xb - Xm)^(3 / 2)) / (Xm - Xa)^(3 / 2))) {
        Ye[i,3] <- alfa * Xm - beta * Xb + alfa * sqrt((Xb-Xm) * (Xm-Xa)) + gama
        Yaest   <- alfa * Xa - beta * Xb + gama
        Ybest   <- alfa * Xb - beta * Xa + gama
        Ye[i,1] <- (Ybest + Yaest)/2
        Ye[i,2] <- (Ybest - Yaest)/2
      }
    } else if(round((Xm - Xa),4) == round((Xb - Xm),4)) {
      Yaest   <- alfa * Xa - beta * Xb + gama
      Ybest   <- alfa * Xb - beta * Xa + gama
      Ye[i,1] <- (Ybest + Yaest)/2
      Ye[i,2] <- (Ybest - Yaest)/2
      Ye[i,3] <- (Ybest + Yaest)/2
    } # if

  } # for

  return(Ye)

}
