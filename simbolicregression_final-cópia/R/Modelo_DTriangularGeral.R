RLParametrosModelo3 <- function(Y, X) {

  #  Modelo3 -- Distribuição TriangularGeral
  #  Y -- Variável Resposta
  #         Y -> yc,yr,ym
  #  X -- p Variáveis Explicativas
  #         X -> xc1,xr1,xm1,xc2,xr2,,xm2,...,xcp,xrp,xmp

  ################################################################################################################################

  initialRL <- c(10,10,10)
  RLvar <- c(55,55,55)

  SolOpt <- RepLOptim(parmean = initialRL,parsd = RLvar, fr = OptFuction,lower = c(0,0,-Inf), Y = Y, X = X)

  return(SolOpt)
}

