NewModelo <- function(Y, X, Xe = X) {

  RL <- RLParametrosModelo1(Y, X)
  Ye <- RLEstimar(Xe, RL)
  DM <- DistMallows(1, Y, Ye)
  #return(ModeloFinal)

}
