clusterwiseRegression <- function(dataset, nClusters = 2, nIterations = 40, rowLabel = F, nInitializations = 10) {

  if (((ncol(dataset) - 1) %% 2 == 0 && rowLabel == T) || (ncol(dataset) %% 2 == 0 && rowLabel == F)) {

    #saving row labels with row identifiers
    if (rowLabel == T) {
      rowLabels <- as.data.frame(dataset[,1])
    } else {
      rowLabels <- NULL
    }
    #convert into centres and half ranges
    datasetCR <- as.data.frame(TOCENTRERANGE(dataset[,-1]))

    #save some measures for later use (to calculate omega' etc)
    Yi <- datasetCR[,c(ncol(datasetCR) - 1, ncol(datasetCR))]
    dYYmean <- getDistance(Yi)
    RLid <- RLParametrosModelo1(Yi, datasetCR[,c(-ncol(datasetCR) + 1, -ncol(datasetCR))])
    Yid <- RLEstimar(datasetCR[,c(-ncol(datasetCR) + 1, -ncol(datasetCR))], RLid) #as.data.frame
    dYYid <- getDistance(Yi, Yid)

    #intial parition of data into nCLusters
    clusters <- split(datasetCR, sample(1:nClusters, nrow(datasetCR), replace = T))

    #beginning nIterations
    for (i in 1:nIterations) {

      tempParameters = getParameters(clusters, nClusters)
      if (is.na(tempParameters[1])) {
        if (i == 1) {
          tired = 0 #force find an acceptable partition
          while(is.na(tempParameters[1])) {
            clusters <- split(datasetCR, sample(1:nClusters, nrow(datasetCR), replace = T))
            tempParameters = getParameters(clusters, nClusters)
            tired = tired + 1
            if (tired == nInitializations) #unable to find an acceptable partition 7 times in a row, give up
              return(NA)
              #stop("Unable to find an acceptable initial parition, check that data and number of clusters is reasonable")
          }
        }
        else {
          qualityMeasures <- getQualityMeasures(prevClusters, dYYmean, dYYid, Yi, as.data.frame(Yid), estimates, nClusters, nrow(datasetCR))
          mergedClusters <- mergeLabel(prevClusters, rowLabels, nClusters)
          return(list("Parameters" = parameters, "Clusters" = mergedClusters, "Quality Measures" = qualityMeasures, "Iteration" = i - 1))
        }
      }
      parameters <- tempParameters
      estimates <- getEstimates(clusters, nClusters, parameters)
      prevClusters <- clusters  #fallback cluster in case in the next iteration parameters cannot be solved
      tempClusters <- vector("list", nClusters)

      #for each cluster
      for (j in 1:nClusters) {
        clusterData = clusters[[j]]
        X = clusterData[,c(-ncol(clusterData) + 1, -ncol(clusterData))]
        Y = clusterData[,c(ncol(clusterData) - 1, ncol(clusterData))]
        #for each row in cluster calculate mallows distance between observed and estimated
        for (k in 1:nrow(clusterData)) {
          observed = Y[k,]
          mallowsDistances <- list()
          #compute mallows distance
          for (d in 1:ncol(estimates)) {
            estimate = estimates[[j,d]][[1]][k,] #or  as.matrix(Y[k,])
            mallowsDistances[d] <- DistMallows(1, as.matrix(observed), t(as.matrix(estimate)))
          }
          #reassign clusters
          clusterIndex <- which.min(mallowsDistances)
          tempClusters[[clusterIndex]] <- rbind(tempClusters[[clusterIndex]], cbind(X[k,], Y[k,]))
        }
      }

      if(isTRUE(all.equal(tempClusters, clusters))) {
        clusters <- tempClusters
        qualityMeasures <- getQualityMeasures(clusters, dYYmean, dYYid, Yi, as.data.frame(Yid), estimates, nClusters, nrow(datasetCR))
        mergedClusters <- mergeLabel(clusters, rowLabels, nClusters)
        return(list("Parameters" = parameters, "Clusters" = mergedClusters, "Quality Measures" = qualityMeasures, "Iteration" = i))
      } else {
          clusters <- tempClusters
      }

      if (i == nIterations) {
        parameters <- getParameters(clusters, nClusters)
        if (is.na(parameters[1])) {
          parameters <- getParameters(prevClusters, nClusters)
        }
        estimates <- getEstimates(clusters, nClusters, parameters)
        qualityMeasures <- getQualityMeasures(clusters, dYYmean, dYYid, Yi, as.data.frame(Yid), estimates, nClusters , nrow(datasetCR))
        mergedClusters <- mergeLabel(clusters, rowLabels, nClusters)
        warning("Last iteration reached")
        return(list("Parameters" = parameters, "Clusters" = mergedClusters, "Quality Measures" = qualityMeasures))
      }
    }
  } else {
      stop("Number of columns of dataset without rowlabels/identifiers should be divisible by 2")
  }

}
