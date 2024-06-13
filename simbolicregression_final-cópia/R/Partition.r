partition <- function(dataset) {
  #partition of data in two clusters(cold season and hot season)
  coolSeason <- dataset[which(dataset$Month %in% c(1,2,3,4,12)),]
  hotSeason <- dataset[which(dataset$Month %in% c(5,6,7,8,9,10,11)),]

  #converting partition into centre-halfranges
  coolSeason <- TOCENTRERANGE(coolSeason[,c(4:ncol(coolSeason))])
  hotSeason <- TOCENTRERANGE(hotSeason[,c(4:ncol(hotSeason))])

  for(i in 1:3) { #or repeat

    #creating X and Y for initial partition
    coolSeasonX <- coolSeason[,c(-ncol(coolSeason) + 1, -ncol(coolSeason))] #t1
    coolSeasonY <- coolSeason[,c(ncol(coolSeason) - 1, ncol(coolSeason))]

    hotSeasonX <- hotSeason[,c(-ncol(hotSeason) + 1, -ncol(hotSeason))] #t2
    hotSeasonY <- hotSeason[,c(ncol(hotSeason) - 1, ncol(hotSeason))]

    #finding model for each cluster
    RL1 <- RLParametrosModelo1(coolSeasonY, coolSeasonX)
    RL2 <- RLParametrosModelo1(hotSeasonY, hotSeasonX)

    #estimations using models
    Ye1 <- RLEstimar(coolSeasonX, RL1)
    Ye2 <- RLEstimar(coolSeasonX, RL2)
    Ye3 <- RLEstimar(hotSeasonX, RL1)
    Ye4 <- RLEstimar(hotSeasonX, RL2)

    #temp clusters
    c1 <- NULL
    c2 <- NULL

    #reassignment
    for(i in 1:nrow(coolSeasonY)) {
      d1 <- DistMallows(1, t(as.matrix(coolSeasonY[i,])), t(as.matrix(Ye1[i,])))
      d2 <- DistMallows(1, t(as.matrix(coolSeasonY[i,])), t(as.matrix(Ye2[i,])))
      if(d1 > d2) {
        c2 <- rbind(c2, cbind(t(coolSeasonX[i,]), t(coolSeasonY[i,])))
      }
      else {
        c1 <- rbind(c1, cbind(t(coolSeasonX[i,]), t(coolSeasonY[i,])))
      }
    }

    for(i in 1:nrow(hotSeasonY)) {
      d1 <- DistMallows(1, t(as.matrix(hotSeasonY[i,])), t(as.matrix(Ye3[i,])))
      d2 <- DistMallows(1, t(as.matrix(hotSeasonY[i,])), t(as.matrix(Ye4[i,])))
      if(d1 > d2) {
        c2 <- rbind(c2, cbind(t(hotSeasonX[i,]), t(hotSeasonY[i,])))
      }
      else {
        c1 <- rbind(c1, cbind(t(hotSeasonX[i,]), t(hotSeasonY[i,])))
      }
    }

    #closing condition
    if (all.equal(c1, coolSeason) == TRUE && all.equal(c2, hotSeason) == TRUE) {
      coolSeason <- c1
      hotSeason <- c2
      break
    } else {
      coolSeason <- c1
      hotSeason <- c2
    }

  }
  return(list(c("coolSeasonParam"= RL1, "hotSeasonParam"= RL2,"c1"=nrow(coolSeason), "c2"=nrow(hotSeason))))
}
