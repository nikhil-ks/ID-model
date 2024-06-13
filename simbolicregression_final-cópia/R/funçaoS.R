S <- function(initialRL,Y,X) {

  Y <- as.matrix(Y)
  X <- as.matrix(X)

  YC <- as.matrix(Y[, 1])                                      #matriz com os centros da variável resposta
  YR <- as.matrix(Y[, 2])                                      #matriz com os raios da variável resposta
  YM <- as.matrix(Y[, 3])                                      #matriz com as modas da variável resposta
  XC <- as.matrix(X[, seq(from = 1, to = ncol(X), by = 3)])    #matriz com os centros das variaveis explicativas
  XR <- as.matrix(X[, seq(from = 2, to = ncol(X), by = 3)])    #matriz com os raios das variáveis explicativas
  XM <- as.matrix(X[, seq(from = 3, to = ncol(X), by = 3)])    #matriz com as modas das variáveis explicativas
  p  <- ncol(XC)                                               #número de variáveis
  n  <- nrow(X)                                                #número de observações

  alfa  <- initialRL[1]
  beta  <- initialRL[2]

  try(if(p > 1) stop("*****model only supports one variable*****"))

##################################################################################################################

  sum   <- 0
  DM    <- 0
  Ycest <- 0
  Yrest <- 0
  Ymest <- 0
  t     <- 0



  for(i in 1:n) {
    Yc <- YC[i, ]
    Yr <- YR[i, ]
    Ym <- YM[i, ]

    Xc <- XC[i, ]
    Xr <- XR[i, ]
    Xm <- XM[i, ]
    Xa <- Xc - Xr
    Xb <- Xc + Xr

    if((Xm-Xa) > (Xb-Xm)) {
      t <- (alfa^(2/3))/((alfa^(2/3)) + beta^(2/3))
      if((t>((Xb-Xm)/(Xb-Xa)))&(t<((Xm-Xa)/(Xb-Xa)))){
        Ymest <- (alfa-beta)*Xa + sqrt((Xb-Xa)*(Xm-Xa))*((alfa^(4/3)-beta^(4/3))/sqrt(alfa^(2/3)+beta^(2/3)))
        Yaest <- alfa * Xa - beta * Xb
        Ybest <- alfa * Xb - beta * Xa
      } else if(alfa < ((beta*(Xb-Xm)^(3/2))/(Xm-Xa)^(3/2))) {
        Ymest <- alfa * Xa - beta * Xm + alfa * sqrt((Xb-Xm) * (Xm-Xa))
        Yaest <- alfa * Xa - beta * Xb
        Ybest <- alfa * Xb - beta * Xa
      } else if(alfa > ((beta*(Xm-Xa)^(3/2))/(Xb-Xm)^(3/2))) {
        Ymest <- alfa * Xm - beta * Xa + beta * sqrt((Xb-Xm) * (Xm-Xa))
        Yaest <- alfa * Xa - beta * Xb
        Ybest <- alfa * Xb - beta * Xa
      }
    }
    else {
      if((Xm-Xa) < (Xb-Xm)) {
        t <- (beta^(2/3))/((alfa^(2/3)) + beta^(2/3))
        if((t>((Xm-Xa)/(Xb-Xa)))&(t<((Xb-Xm)/(Xb-Xa)))){
          Ymest <- (alfa-beta)*Xb + sqrt((Xb-Xa)*(Xb-Xm))*((beta^(4/3)-alfa^(4/3))/sqrt(alfa^(2/3)+beta^(2/3)))
          Yaest <- alfa * Xa - beta * Xb
          Ybest <- alfa * Xb - beta * Xa
        } else if(alfa < ((beta*(Xm-Xa)^(3/2))/(Xb-Xm)^(3/2))) {
          Ymest <- alfa * Xb - beta * Xm + alfa * sqrt((Xb-Xm) * (Xm-Xa))
          Yaest <- alfa * Xa - beta * Xb
          Ybest <- alfa * Xb - beta * Xa
        } else if(alfa > ((beta*(Xb-Xm)^(3/2))/(Xm-Xa)^(3/2))) {
          Ymest <- alfa * Xm - beta * Xb + alfa * sqrt((Xb-Xm) * (Xm-Xa))
          Yaest <- alfa * Xa - beta * Xb
          Ybest <- alfa * Xb - beta * Xa
        } else if((Xm-Xa) < (Xb-Xm)) {
          Yaest <- alfa * Xa - beta * Xb
          Ybest <- alfa * Xb - beta * Xa
          Ymest <- (Ybest + Yaest)/2
        }
      }
    }

    Ycest <- (Ybest + Yaest)/2
    Yrest <- (Ybest - Yaest)/2

    DM <- DistMallows2(Yc,Yr,Ym,Ycest,Yrest,Ymest)
    sum <- sum + DM # CONFIRMAR
  }

  return(sum)
}
