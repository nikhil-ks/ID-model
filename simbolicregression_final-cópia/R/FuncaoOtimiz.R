OptFuction <- function(initialRL,Y,X) {

  Y <- as.matrix(Y)
  X <- as.matrix(X)

  YC <- as.matrix(Y[, 1])                                      #matriz com os centros da variável resposta
  YR <- as.matrix(Y[, 2])                                      #matriz com os raios da variável resposta
  YM <- as.matrix(Y[, 3])                                      #matriz com as modas da variável resposta
  XC <- as.matrix(X[, 1])    #matriz com os centros das variaveis explicativas
  XR <- as.matrix(X[, 2])    #matriz com os raios das variáveis explicativas
  XM <- as.matrix(X[, 3])    #matriz com as modas das variáveis explicativas                                             #número de variáveis
  p  <- ncol(XC)
  n  <- nrow(X)                                                #número de observações

  alfa  <- initialRL[1]
  beta  <- initialRL[2]
  gama  <- initialRL[3]

  try(if(p > 1) stop("model only supports one variable"))

  #---------------------------------------------------------------------------------------------------------------------#

  sum <- 0

  for(i in 1:n) {

    Yc <- 0
    Yr <- 0
    Ym <- 0
    Xc <- 0
    Xr <- 0
    Xm <- 0

    Yc <- YC[i, ]
    Yr <- YR[i, ]
    Ym <- YM[i, ]

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
        Ymest <- (alfa-beta)*Xa + sqrt((Xb-Xa)*(Xm-Xa))*((alfa^(4/3)-beta^(4/3))/sqrt(alfa^(2/3)+beta^(2/3))) + gama
        Yaest <- alfa * Xa - beta * Xb + gama
        Ybest <- alfa * Xb - beta * Xa + gama
        Ycest <- (Ybest + Yaest)/2
        Yrest <- (Ybest - Yaest)/2
      } else if(alfa <= ((beta * (Xb - Xm) ^ (3 / 2)) / (Xm - Xa)^(3 / 2))) {
        Ymest <- alfa * Xa - beta * Xm + alfa * sqrt((Xb-Xm) * (Xm-Xa)) + gama
        Yaest <- alfa * Xa - beta * Xb + gama
        Ybest <- alfa * Xb - beta * Xa + gama
        Ycest <- (Ybest + Yaest)/2
        Yrest <- (Ybest - Yaest)/2
      } else if(alfa >= ((beta * (Xm - Xa)^(3 / 2)) / (Xb - Xm)^(3 / 2))) {
        Ymest <- alfa * Xm - beta * Xa + beta * sqrt((Xb-Xm) * (Xm-Xa)) + gama
        Yaest <- alfa * Xa - beta * Xb + gama
        Ybest <- alfa * Xb - beta * Xa + gama
        Ycest <- (Ybest + Yaest)/2
        Yrest <- (Ybest - Yaest)/2
      }

    } else if(round((Xm - Xa),4) < round((Xb - Xm),4)) {
      #t <- (beta^(2 / 3)) / ((alfa^(2 / 3)) + beta^(2 / 3))
      if((alfa > ((beta * (Xm - Xa)^(3 / 2)) / (Xb - Xm)^(3 / 2))) & (alfa < ((beta * (Xb - Xm)^(3 / 2)) / (Xm - Xa)^(3 / 2)))) {
        Ymest <- (alfa-beta)*Xb + sqrt((Xb-Xa)*(Xb-Xm))*((beta^(4/3)-alfa^(4/3))/sqrt(alfa^(2/3)+beta^(2/3))) + gama
        Yaest <- alfa * Xa - beta * Xb + gama
        Ybest <- alfa * Xb - beta * Xa + gama
        Ycest <- (Ybest + Yaest)/2
        Yrest <- (Ybest - Yaest)/2
      } else if(alfa <= ((beta * (Xm - Xa)^(3 / 2)) / (Xb - Xm)^(3 / 2))) {
        Ymest <- alfa * Xb - beta * Xm + alfa * sqrt((Xb-Xm) * (Xm-Xa)) + gama
        Yaest <- alfa * Xa - beta * Xb + gama
        Ybest <- alfa * Xb - beta * Xa + gama
        Ycest <- (Ybest + Yaest)/2
        Yrest <- (Ybest - Yaest)/2
      } else if(alfa >= ((beta * (Xb - Xm)^(3 / 2)) / (Xm - Xa)^(3 / 2))) {
        Ymest <- alfa * Xm - beta * Xb + alfa * sqrt((Xb-Xm) * (Xm-Xa)) + gama
        Yaest <- alfa * Xa - beta * Xb + gama
        Ybest <- alfa * Xb - beta * Xa + gama
        Ycest <- (Ybest + Yaest)/2
        Yrest <- (Ybest - Yaest)/2
      }
    } else if(round((Xm - Xa),4) == round((Xb - Xm),4)) {
      Yaest <- alfa * Xa - beta * Xb + gama
      Ybest <- alfa * Xb - beta * Xa + gama
      Ycest <- (Ybest + Yaest)/2
      Yrest <- (Ybest - Yaest)/2
      Ymest <- Ycest
    } # if

  DM <- DistMallows2(Yc,Yr,Ym,Ycest,Yrest,Ymest)
  sum <- sum + DM
  } # for

  return(sum)

} # function
