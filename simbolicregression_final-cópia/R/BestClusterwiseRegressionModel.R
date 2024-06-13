bestClusterRegressionModel <- function(dataset, nInitializations = 10, nClusters = 2, rowLabel = F) {
  bestModel <- list()
  for (i in 1:nInitializations) {
    tempModel <- clusterwiseRegression(dataset = dataset, nClusters = nClusters, rowLabel = rowLabel, nInitializations = nInitializations)
    if (i == 1) {
      if(!is.na(tempModel[1]))
        bestModel <- tempModel
      else
        stop("Unable to find an acceptable initial partition, check that data and number of clusters is reasonable")
    }
    else if (!is.na(tempModel[1]) && tempModel[["Quality Measures"]][[nClusters + 1]][1] < bestModel[["Quality Measures"]][[nClusters + 1]][1])
      bestModel <- tempModel
    else if (!is.na(tempModel[1]) && tempModel[["Quality Measures"]][[nClusters + 1]][1] == bestModel[["Quality Measures"]][[nClusters + 1]][1])
      if (tempModel[["Quality Measures"]][[nClusters + 5]][1] > bestModel[["Quality Measures"]][[nClusters + 5]][1])
        bestModel <- tempModel
    else if (is.na(tempModel[1])) {
      warning("returning best model due to lack of solutions")
      return(bestModel)
    }
  }
  return(bestModel)
}
