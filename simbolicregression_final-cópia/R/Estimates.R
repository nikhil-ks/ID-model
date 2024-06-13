getEstimates <- function(clusters, n, parameters) {

  estimates <- matrix(list(), nrow = n, ncol = n)

  #estimates for each cluster using obtained parameters
  for (i in 1:n) {
    clusterData = clusters[[i]]
    Xe = clusterData[,c(-ncol(clusterData) + 1, -ncol(clusterData))]
    #for each paramter obtained
    for (j in 1:length(parameters)) {
      estimates[[i,j]] <- list(RLEstimar(Xe, parameters[[j]]))
    }
  }
  return(estimates)
}
