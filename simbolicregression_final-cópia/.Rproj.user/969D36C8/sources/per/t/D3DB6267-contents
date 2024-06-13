getQualityMeasures <- function(clusters, dYYmean, dYYid, Yi, Yid, estimates, nClusters, nRows) {

  qualityMeasures <- list()
  weightedCD <- 0
  sumOmega <- 0
  sumResiduals <- 0
  totalSilhouetteCoef <- 0
  silhouetteCoefs <- vector("list", nClusters)
  clusterCoefs <- list()

  #calculate quality measues for each cluster
  for (i in 1:nClusters) {
    qualityMeasures[i] <- RLQualidade(1, clusters[[i]][,c(ncol(clusters[[i]]) - 1, ncol(clusters[[i]]))], estimates[[i,i]][[1]])
    weightedCD <- weightedCD + ((nrow(clusters[[i]]) / nRows) * qualityMeasures[[i]][1])
    sumOmega <- sumOmega + getDistance(estimates[[i,i]][[1]], cbind(matrix(rep(mean(clusters[[i]][,ncol(clusters[[i]]) - 1]), nrow(clusters[[i]]))), 0))
    sumResiduals <- sumResiduals + getDistance(estimates[[i,i]][[1]], Yid[rownames(clusters[[i]]),])
    totalClusterSilhouetteCoef <- 0
    # sum silhouette coefficient
    for (j in 1:nrow(clusters[[i]])) {
      ai = DistMallows(1, as.matrix(Yi[rownames(clusters[[i]])[j],]), t(as.matrix(estimates[[i,i]][[1]][j,])))
      mallowsDistances <- vector()
      for (k in 1:ncol(estimates)) {
        if (k != i)
          mallowsDistances[k] <- DistMallows(1, as.matrix(Yi[rownames(clusters[[i]])[j],]), t(as.matrix(estimates[[i,k]][[1]][j,])))
      }
      bi = min(mallowsDistances, na.rm = T)
      totalSilhouetteCoef <- totalSilhouetteCoef + ((bi - ai) / max(ai, bi))
      totalClusterSilhouetteCoef <- totalClusterSilhouetteCoef + ((bi - ai) / max(ai, bi))
      silhouetteCoefs[[i]] <- rbind(silhouetteCoefs[[i]], cbind(rownames(clusters[[i]])[j], (bi - ai) / max(ai, bi)))
    }
    clusterCoefs[i] <- totalClusterSilhouetteCoef / nrow(clusters[[i]])
  }

  DMrmsem <- 0
  globalRMSEL <- 0
  globalRMSEU <- 0
  goodnessOfFit <- 0

  #calculate global RMSEs
  for (i in 1:nClusters) {
    DMrmsem <- DMrmsem + (nrow(clusters[[i]])*(qualityMeasures[[i]][2] ^ 2))
    globalRMSEL <- globalRMSEL + (nrow(clusters[[i]])*(qualityMeasures[[i]][3] ^ 2))
    globalRMSEU <- globalRMSEU + (nrow(clusters[[i]])*(qualityMeasures[[i]][4] ^ 2))
  }

  qualityMeasures[length(qualityMeasures) + 1][1] <- list(c("Goodness of Fit" = DMrmsem))
  qualityMeasures[length(qualityMeasures) + 1][1] <- list(c("Weighted CD" = weightedCD))
  qualityMeasures[length(qualityMeasures) + 1][1] <- list(c("CD'" = (sumOmega / dYYmean)))
  #Ratio between residual inter-class variance and total residual variance
  qualityMeasures[length(qualityMeasures) + 1][1] <- list(c("C" = (sumResiduals / dYYid)))
  qualityMeasures[length(qualityMeasures) + 1][1] <- list(c("Global Silhouette Coefficient" = totalSilhouetteCoef / nRows))
  qualityMeasures[length(qualityMeasures) + 1][1] <- list(c("Class Silhouette Coefficients" = clusterCoefs))
  qualityMeasures[length(qualityMeasures) + 1][1] <- list(c("Global RMSEM" = sqrt(DMrmsem/nRows), "Global RMSEL" = sqrt(globalRMSEL/nRows), "Global RMSEU" = sqrt(globalRMSEU/nRows)))
  qualityMeasures[length(qualityMeasures) + 1][1] <- list(c("Silhouette Coefficients" = silhouetteCoefs))

  return(qualityMeasures)
}
