getMode <- function(x){
  xdens <- density(x)
  modex <- xdens$x[which.max(xdens$y)]
  return(modex)
}

ToSimbolic <- function(X,m) {
 # agregador <- unique(X[,1])
  agregador <- X[,1]

  if(m == 1 | m == 2) {

    simbolic <- matrix(data = 0, nrow = length(agregador), ncol = 2)

    for(i in 1:length(agregador)) {
      data <- 0
      data <- X[X[,1] == agregador[i],]
      simbolic[i,1] <- (max(data[,2]) + min(data[,2]))/2
      simbolic[i,2] <- (max(data[,2]) - min(data[,2]))/2
    }

    return(simbolic)

  } else if(m == 3) {

    simbolic <- matrix(data = 0, nrow = length(agregador), ncol =  3)

    for(i in 1:length(agregador)) {
      data <- 0
      data <- X[X[,1] == agregador[i],]
      simbolic[i,1] <- (max(data[,2]) + min(data[,2]))/2
      simbolic[i,2] <- (max(data[,2]) - min(data[,2]))/2
      if(nrow(data)==1) {simbolic[i,3] <- (max(data[,2]) + min(data[,2]))/2}
      else{simbolic[i,3] <- getMode(data[,2])}
    }
    return(simbolic)
  } else {

    return("Escolher 1 para utilizar o modelo assumindo uma distribuição uniforme, 2 para utilizar o modelo assumindo uma distribuição triangular simétrica ou 3 para utilizar o modelo assumindo uma distribuição triangular geral")
  }
}

TOSIMBOLIC <- function(X,m) {

  var <- ncol(X)-1
  A <- NULL

  for(i in 1:var) {
    D <- ToSimbolic(X[,c(1,i+1)],m)

    A <- cbind(A,D)
  }
  return(A)
}
