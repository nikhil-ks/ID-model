getParameters <- function(clusters, n) {

  parameters <- list()

  #for each cluster, calculate parameters
  for (i in 1:n) {
    if (length(clusters) < n || is.null(clusters[[i]]))
      return(NA)
    clusterData = clusters[[i]]
    X = clusterData[,c(-ncol(clusterData) + 1, -ncol(clusterData))]
    Y = clusterData[,c(ncol(clusterData) - 1, ncol(clusterData))]
    parameters[i] = list(RLParametrosModelo1(Y, X))
    if (is.na(parameters[i]))
      return(NA)
  }
  return(parameters)
}
