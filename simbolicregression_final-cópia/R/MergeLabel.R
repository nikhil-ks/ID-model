mergeLabel <- function(clusters, rowLabels, n) {

  if (is.null(rowLabels)) {
    return(clusters)
  }

  mergedClusters <- list()

  #for each cluster, cbind with identifier of observations
  for (i in 1:n) {
    mergedClusters[[i]] <- cbind(label = rowLabels[rownames(clusters[[i]]),], clusters[[i]])
  }
  return(mergedClusters)
}
