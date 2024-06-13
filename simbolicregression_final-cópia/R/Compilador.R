Modelo <- function(Model, Y, X, Xe = X) {

  if(Model == 1) {

    RL <- RLParametrosModelo1(Y, X)
    Ye <- RLEstimar(Xe, RL)
    Metricas <- RLQualidade(Model, Y, Ye)
    ModeloFinal <-list("Parametros do Modelo" = RL,"Head Intervalos Estimados" = head(Ye),  "Medidas de Qualidade do Modelo" = Metricas)

    return(ModeloFinal)

  } else if(Model == 2) {

    RL <- RLParametrosModelo2(Y, X)
    Ye <- RLEstimar(Xe, RL)
    Metricas <- RLQualidade(Model, Y, Ye)
    ModeloFinal <- list("Parametros do Modelo" = RL,"Head Intervalos Estimados" = head(Ye),  "Medidas de Qualidade do Modelo" = Metricas)

    return(ModeloFinal)

  } else if(Model == 3) {

    RL <- RLParametrosModelo3(Y,X)$par
    Ye <- Modelo3Previsao(Xe,RL)
    Metricas <- RLQualidade2(Y,Ye)
    ModeloFinal <-list("Parametros do Modelo" = RL,"Head Intervalos Estimados" = Ye,  "Medidas de Qualidade do Modelo" = Metricas)

    return(ModeloFinal)

  }

}

