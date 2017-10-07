#' @include dissever-package.R

# Horrible hack to keep CRAN happy and suppress NOTES about parts of the code that use non-standard evaluation.
# See:
# http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
# https://github.com/smbache/magrittr/issues/29
utils::globalVariables(c( "cell", "diss", ".", "matches", "i"))

# simple wrapper around raster::as.data.frame that handles categorical data columns correctly
.as_data_frame_factors <- function(x, ...) {
  res <- as.data.frame(x, ...)
  if (any(is.factor(x))) {
    # Get names of original stack (not affected by bug)
    nm <- names(x)
    # Test if coordinates have been returned
    args <- list(...)
    if("xy" %in% names(args)) { xy <- args[['xy']] } else { xy <- FALSE }
    # If coordinates asked to be returned, we need to take this into account
    if (xy) { names(res)[-1*1:2] <- nm } else { names(res) <- nm }
  }
  res
}

.join_interpol <- function(coarse_df, fine_df, attr, by = 'cell') {
  # Nearest-neighbour interpolation as an inner SQL join
  left_join(fine_df, coarse_df , by = by) %>% select(matches(attr))
}

.create_lut_fine <- function(coarse, fine) { extract(coarse, coordinates(fine)) }

.default_control_init <- caret::trainControl( method = 'cv', number = 5 )

# fits one model to the entire training set
.default_control_iter <- caret::trainControl(method = 'none')

.create_tune_grid <- function(model, tune_length) {
  params <- modelLookup(model)$parameter
  grid <- expand.grid(lapply(1:length(params), function(x) 1:tune_length))
  names(grid) <- as.character(params)
  grid
}

# Computes the carret regression model between some coarse data and the stack of covariates
.update_model <- function(vars, y, method = 'rf', control, tune_grid, data_type="numeric", latLong=NULL) {
  # Pick the parameters of the model using error on the first run.
  # Then use the optimised parameters in the iteration loop to save on computing time.
  # Basically we just need to change the trainControl object to do that.
  y_aux = y
  if ( data_type == 'categorical' ) { y_aux = factor( y_aux ) }
  if ( method == 'gwr' ) { 
    fit <- gw( as.formula(paste("x~",paste(names(vars), collapse="+"))) , data= data.frame( vars , x=y_aux ) )
  } else if ( method == 'mlp' ) {
    data <- mx.symbol.Variable("data")
    fc1 <- mx.symbol.FullyConnected(data, num_hidden=10)
    act1 <- mx.symbol.Activation(fc1, act_type="relu")
    fc2 <- mx.symbol.FullyConnected(act1, num_hidden=1)
    fit <- mx.model.FeedForward.create(mx.symbol.LinearRegressionOutput(fc2), X=vars, y=y_aux, num.round=50, array.batch.size=20, learning.rate=2e-6, momentum=0.9, eval.metric=mx.metric.rmse)
  } else if ( method == 'lme' ) {
    fit <- data.frame( vars , out=y_aux , lat=latLong$lat , long=latLong$long , dummy=rep.int( 1 , length(y_aux) ) )
    fit <- lme( fixed=out ~ . - dummy - lat - long , data=fit , random = ~ 1 | dummy, correlation = corGaus(form = ~ lat+long | dummy ) )
  } else fit <- train( x = vars, y = y_aux, method = method, trControl = control, tuneGrid  = tune_grid )
  fit
}

.has_parallel_backend <- function() getDoParWorkers() > 1

.get_split_idx <- function(n, p) sort(1:n %% p)

# Generates prediction intervals using bootstraping
.bootstrap_ci <- function(fit, fine_df, level = 0.9, n = 50L, data_type="numeric", latLong=NULL ) {
  df <- fit$trainingData
  reg_method <- fit$method
  boot_samples <- boot(df, function(data, idx, method = reg_method) {
    bootstrap_df <- data[idx, ]
    # if ( data_type == "categorical" ) { bootstrap_df <- factor(bootstrap_df) }
    if ( method == 'gwrm' ) bootstrap_fit <- gw(.outcome ~ ., data = bootstrap_df )
    else if ( method == 'mlp' ) {
      data <- mx.symbol.Variable("data")
      fc1 <- mx.symbol.FullyConnected(data, num_hidden=10)
      act1 <- mx.symbol.Activation(fc1, act_type="relu")
      fc2 <- mx.symbol.FullyConnected(act1, num_hidden=1)
      # TODO : fix for MLP
      # bootstrap_fit <- mx.model.FeedForward.create(mx.symbol.LinearRegressionOutput(fc2), X=vars, y=y_aux, num.round=50, array.batch.size=20, learning.rate=2e-6, momentum=0.9, eval.metric=mx.metric.rmse)      
      # bootstrap_fit <- train(.outcome ~ ., data = bootstrap_df, method = method, trControl = trainControl(method = "none"), tuneGrid = fit$bestTune)
    } else if ( method == 'lme' ) bootstrap_fit <- lme( fixed=.outcome ~ . - dummy - lat - long , data=data.frame( bootstrap_df , lat=latLong$lat[idx], lat=latLong$long[idx], dummy=rep.int( 1 , nrow(bootstrap_df) ) ) , random = ~ 1 | dummy, correlation = corGaus(form = ~ lat+long | dummy ) )
    else bootstrap_fit <- train(.outcome ~ ., data = bootstrap_df, method = method, trControl = trainControl(method = "none"), tuneGrid = fit$bestTune)
    # generate predictions
    if ( method == 'lme' ) predict(bootstrap_fit, data.frame(fine_df,latLong,dummy=rep.int(1,nrow(fine_df))) )
    else { 
      res <- predict(bootstrap_fit, fine_df)
      if ( !is.null( nrow(res) ) ) res <- res[,1]
      as.numeric( res )
    }
  }, n)
  # If level is a number < 1, or instead if level is a function
  if (is.numeric(level)) {
    ci <- c((1 - level) / 2, 1 - (1 - level) / 2)
    if ( data_type == "categorical" ) {
      res <- data.frame(
        lower = aaply(boot_samples$t, 2, quantile, probs = ci[1]),
        mean = aaply(boot_samples$t, 2, median),
        upper = aaply(boot_samples$t, 2, quantile, probs = ci[2])
      )
    } else {
      res <- data.frame(
        lower = aaply(boot_samples$t, 2, quantile, probs = ci[1]),
        mean = aaply(boot_samples$t, 2, mean),
        upper = aaply(boot_samples$t, 2, quantile, probs = ci[2])
      )
    }
  } else if (is.function(level)) {
    if ( data_type == "categorical" ) {
      res <- data.frame( mean = aaply(boot_samples$t, 2, median), uncert = aaply(boot_samples$t, 2, level) )
    } else {
      res <- data.frame( mean = aaply(boot_samples$t, 2, mean), uncert = aaply(boot_samples$t, 2, level ) )
    }
  } else {
    stop('Incorrect value for the "level" option.')
  }
  as.numeric( res )
}

.generate_ci <- function(object, covariates, level = 0.9, n = 50L, latLong=NULL) {
  fine_df <- na.exclude(.as_data_frame_factors(covariates, xy = TRUE))
  b <- .bootstrap_ci(object$fit, fine_df, level = 0.9, n = 50L, latLong=latLong)
  res <- rasterFromXYZ( data.frame( fine_df[, 1:2], b ), res = res(covariates), crs = projection(covariates) )
  res
}

.predict_map <- function(fit, data, split = NULL, boot = NULL, level = 0.9, data_type="numeric", latLong=NULL) {
  if (.has_parallel_backend()) {
    n_workers <- length(unique(split))
    if (n_workers < 1) stop('Wrong split vector')
    res <- foreach( i = 0:(n_workers - 1), .combine = c, .packages = 'caret' ) %dopar% {
      if (is.null(boot)) {
        if (! is.null(latLong) ) res <- predict(object=fit , newdata = data.frame(data,latLong)[split == i, ])
        else res <- predict(object=fit , newdata = data[split == i, ])
        if ( !is.null( nrow(res) ) ) res <- res[,1]
        as.numeric( res )
      } else {
        .bootstrap_ci(fit = fit, fine_df = data[split == i, ], level = level, n = boot, data_type=data_type, latLong=latLong)
      }
    }
  } else {
    if (is.null(boot)) {
      if (! is.null(latLong) ) res <- predict( object=fit , newdata=data.frame(data,latLong) )
      else res <- predict( object=fit , newdata=data )
      if ( !is.null( nrow(res) ) ) res <- res[,1]
    } else res <- .bootstrap_ci(fit = fit, fine_df = data, level = level, n = boot, data_type=data_type, latLong=latLong)
  }
  as.numeric( res )
}

.gw.dist<- function( dp.locat , rp.locat ) {
  n.rp<-length(rp.locat[,1])
  n.dp<-length(dp.locat[,1])
  dist.res <- big.matrix(nrow=n.dp , ncol=n.rp)
  for (i in 1:n.rp) dist.res[,i]<-spDistsN1(dp.locat, matrix(rp.locat[i,],nrow=1),longlat=TRUE)
  dist.res
}




.pycnoSeq <- function( x, pops, celldim, r=0.2, converge=3, verbose=TRUE ) {
  if (!is(celldim,"SpatialGrid")) {
    
    bbx <- slot(x,'bbox')
    offset <- bbx[,1]
    extent <- bbx[,2] - offset
    shape <- ceiling(extent / celldim)
    gr <- SpatialGrid(GridTopology(offset,c(celldim,celldim),shape)) 
  } else { gr <- celldim }
  
  px <- CRS(proj4string(x))
  proj4string(gr) <- px
  
  gr <- SpatialPixelsDataFrame(coordinates(gr),data.frame(zone=SpatialPoints(coordinates(gr),proj4string=px) %over% as(x,"SpatialPolygons")))
  gr <- as(gr,"SpatialGridDataFrame")
  gr.dim <- slot(getGridTopology(gr), "cells.dim" )
  
  zones <- gr[[1]]
  dim(zones) <- gr.dim
  attr(zones,'na') <- is.na(zones)
  
  zones[is.na(zones)] <- max(zones,na.rm=T) + 1  
  zone.list <- sort(unique(array(zones)))
  pops <- c(pops,0)
  x <- zones * 0
  
  foreach (item = zone.list, .inorder=FALSE, .export = c("x","zones")) %do% {
    zone.set <- (zones == item)
    x[zone.set] <- pops[item] / sum(zone.set)
  }
  
  stopper <- max(x, na.rm = TRUE) * 10^(-converge)
  
  repeat {
    old.x <- x
    mval <- mean(x)
    s1d <- function(s) unclass(stats::filter(s,c(0.5,0,0.5)))
    pad <- rbind(mval,cbind(mval,x,mval),mval)
    pad <- (t(apply(pad,1,s1d)) + apply(pad,2,s1d))/2
    sm <- (pad[2:(nrow(x)+1),2:(ncol(x)+1)])
    x <- x*r + (1-r)*sm
    foreach (item = zone.list, .inorder=FALSE, .export = c("x","zones","pops")) %do% {
      zone.set <- (zones == item)
      correct <- (pops[item] - sum(x[zone.set])) / sum(zone.set)
      x[zone.set] <- x[zone.set] + correct    
    } 
    x[ x < 0 ] <- 0
    foreach (item = zone.list, .inorder=FALSE, .export = c("x","zones","pops")) %do% {
      zone.set <- (zones == item)
      correct <- pops[item] / sum(x[zone.set])
      x[zone.set] <- x[zone.set] * correct 
    }
    if (verbose) {
      flush.console()
      cat(sprintf("Maximum Change: %12.5f - will stop at %12.5f\n", max(abs(old.x - x), na.rm = TRUE), stopper))
    }
    if (max(abs(old.x - x), na.rm = TRUE) <= stopper) break
    print(max(abs(old.x - x), na.rm = TRUE))
  }
  if (!is.null(attr(x,'na'))) x[attr(x,'na')] <- NA
  result <- SpatialPixelsDataFrame( coordinates(gr), data.frame(dens=array(x)) )
  result <- as( result , "SpatialGridDataFrame" )
  proj4string(result) <- px
  return(result)
}

#######################################################################################
#                              Pycno worker functions                                 #
#######################################################################################

.overlap <- function(x, partitions, lockname=NULL, dbpath=NULL){
  
  con <- dbConnect(RSQLite::SQLite(), dbname = dbpath)
  
  sectionShape <- c(globalShape[1], subSection)
  #sectionShape <- globalShape
  
  backCoords <- paste(workerId,"Coords", ".bin", sep = "")
  descCoords <- paste(workerId, "Coords", ".desc", sep = "")
  bigCoords <- filebacked.big.matrix(nrow = sectionShape[[1]]*2, ncol = sectionShape[[2]], type = options()$bigmemory.default.type, separated = FALSE, backingfile = backCoords, descriptorfile = descCoords)
  options(bigmemory.typecast.warning=FALSE)
  
  if(missing(partitions)){
    sizeMb <- 100
    nMaxcellsToLoad <- (sizeMb*1024*1024)/80
    nColsToLoad <- floor(nMaxcellsToLoad/sectionShape[[1]])
    partitions <- ceiling(sectionShape[[2]]/nColsToLoad)
  }
  
  intervaly <- floor(sectionShape/partitions)[[2]]
  subSubSections <- rep(intervaly, partitions)
  
  i <- 1
  while(sum(subSubSections) < sectionShape[2]){
    subSubSections[i] <- subSubSections[i] + 1
    i <- i + 1
  }
  
  px <- CRS(proj4string(poligonDataFrame))
  
  i <- 1
  cursorCoords <- subSection
  while(i <= partitions){
    
    start.time  <- Sys.time()
    
    if( partitions == 1 ){
      #print("1")
      gridElement <- SpatialGrid(GridTopology(offset,c(celldim,celldim),sectionShape))
      rowSumStart <- 0
      ncols=(sectionShape)[[2]]
    } else if (i == 1){
      #print("2")
      gridElement <- SpatialGrid(GridTopology(offset,c(celldim,celldim),c(sectionShape[1], subSubSections[partitions - (i-1)])))
      ncols=subSubSections[partitions - (i-1)]
      rowSumStart <- sum(head(subSubSections,length(subSubSections)-i))
    } else if (i > 1){
      #print("3")
      offset[2] <- offset[2] + (celldim * subSubSections[partitions-(i-2)])
      #rowSum <- rowSum + subSubSections[i-1] - 1
      rowSumStart <- sum(head(subSubSections,length(subSubSections)-i))
      gridElement <- SpatialGrid(GridTopology(offset,c(celldim,celldim),c(sectionShape[1], subSubSections[partitions - (i-1)])))
      ncols=subSubSections[partitions - (i-1)]
    }
    
    #startPosRow <- 1
    #startPosCol <- 1
    #dimCol <- ncol(x)
    #dimRow <- nrow(x)
    
    proj4string(gridElement) <- px
    
    bigCoords[, (cursorCoords-((ncols) - 1)):cursorCoords] <- c(t(coordinates(gridElement)))
    cursorCoords <- cursorCoords - ncols
    
    tryCatch(idZones <- SpatialPixelsDataFrame(coordinates(gridElement),data.frame(zone=SpatialPoints(coordinates(gridElement),proj4string=px) %over% as(poligonDataFrame,"SpatialPolygons"))), 
             
             error = function(c) {
               if(grepl("cannot allocate", c$message) ){
                 g <- 10
                 print(c$message)
                 #partitions <<- partitions + 10
                 try(rm(gridElement))
                 try(rm(idZones))
                 gc()
                 if(i > 1){
                   offset[2] <- offset[2] - (celldim * subSubSections[partitions-(i-2)])
                 }
                 next
                 #idZones <<- .overlap(x)
               }
               print("out") 
               #c$message <- paste0(c$message, "in" )
               #stop(c)
             })
    end.time <- Sys.time()
    time.taken <- difftime(end.time , start.time, units = "secs")
    print(paste(workerId, "- Grid overlaped - partition", i , "of", partitions, "taking",  round(time.taken, digits = 0), "Secs /", "Size of", round(object.size(idZones)/1024/1024, digits = 0), "Mb", "at", Sys.time(), sep = " "))
    rm(gridElement)
    #try(rm(idZones))
    gr.dim <- slot(getGridTopology(idZones), "cells.dim" )
    zones <- idZones[[1]]
    dim(zones) <- gr.dim
    attr(zones,'na') <- is.na(zones)
    zones[is.na(zones)] <- naId
    rm(idZones)
    
    start.time  <- Sys.time()
    for(item in zone.list){
      indices <- which((zones == item), arr.ind = TRUE)
      if(length(indices) == 0){next}
      indices[,2] <- indices[,2] + (colSpanStart - 1) + rowSumStart
      valToInsert <- c(rep(item, nrow(indices)), c(indices))
      valToInsert <- matrix(valToInsert, ncol = 3, nrow = nrow(indices))
      colnames(valToInsert) <- c("idZona", "x", "y")
      
      if(!is.null(lockname)) {
        file.lock = lock(lockname)
      }
      # The lines below are the "critical section"
      dbWriteTable(con, "indices", as.data.frame(valToInsert), append = TRUE)
      
      if(!is.null(lockname)) {
        unlock(file.lock)
      }
      
      #result[[toString(item)]] <- rbind(indices, result[[toString(item)]])
      #result[[toString(item)]] <- indices
      rm(indices)
    }
    end.time <- Sys.time()
    time.taken <- difftime(end.time , start.time, units = "secs")
    print(paste(workerId, "- Indices computed for partition", i , "of", partitions, "- taking",  round(time.taken, digits = 0), "Secs -", "at", Sys.time(), sep = " "))
    rm(zones)
    i <- i + 1
  }
  #return(result)
  dbDisconnect(con)
  return(TRUE)
}

.worker5 <- function(x, index = NULL, dbpath=NULL){
  descx <- paste("x", toString(index), ".desc", sep = "")
  #descZones <- paste("zones", toString(index), ".desc", sep = "")
  x <- attach.big.matrix(descx)
  
  sizeMb <- 200
  maxLengthInds <- (sizeMb*1024*1024)/8
  con <- dbConnect(RSQLite::SQLite(), dbname = dbpath)
  
  if(Sys.info()["sysname"][[1]] == "Windows"){
    #http://www.sqlite.org/faq.html#q5
    res <- dbSendQuery(con, "PRAGMA busy_timeout=10;")
    dbClearResult(res)
  }
  
  countsList <- vector(mode="list", length=length(zoneList))
  names(countsList) <- zoneList
  print(paste0(workerId," - zoneList - ",length(zoneList)))
  print(paste0(workerId," - countsList - ",length(countsList)))
  
  countNrowsQuery <- dbSendQuery(con, "SELECT Count(idZona) FROM indices WHERE idZona = :x")
  print(paste0(workerId," - DONE"))
  for (zone in zoneList) {
    #print(paste(workerId, "- Count Zone =", zone , "at", Sys.time(), sep = " "))
    dbBind(countNrowsQuery, param = list(x = zone))
    lengthInds <- dbFetch(countNrowsQuery)[[1]]
    # if(lengthInds==0){
    #   countsList[[as.character(zone)]] <- NULL
    #   next
    # }
    countsList[[as.character(zone)]] <- lengthInds
    #dbClearResult(countNrowsQuery)
  }
  
  
  print(paste0(workerId," - zoneList - ",length(zoneList)))
  print(paste0(workerId," - countsList - ",length(countsList)))
  dbClearResult(countNrowsQuery)
  assign("countsList",countsList,.GlobalEnv)
  #library(profvis)
  #profvis({
  getIndicesQuery <- dbSendQuery(con, "SELECT x, y FROM indices WHERE idZona = :x")
  for (zone in zoneList) {
    lengthInds <- countsList[[as.character(zone)]]
    if(lengthInds==0){
      next
    }
    popCell <- pops[zone] / lengthInds
    boolPartition <- lengthInds >= maxLengthInds
    dbBind(getIndicesQuery, param = list(x = zone))
    #rest <- dbFetch(getIndicesQuery)
    
    if(!boolPartition){
      indexes <- as.matrix(dbFetch(getIndicesQuery))
      x[indexes] <- popCell
    }else {
      while(!dbHasCompleted(getIndicesQuery)){
        indexes <- as.matrix(dbFetch(getIndicesQuery, n = maxLengthInds))
        x[indexes] <- popCell
      }
    }
    
  }
  print(paste0(workerId," - zoneList - ",length(zoneList)))
  print(paste0(workerId," - countsList - ",length(countsList)))
  dbClearResult(getIndicesQuery)
  #})
  dbDisconnect(con)
  return(TRUE)
}

.maxBigMatrix <- function(x, descriptor){
  bm <- attach.big.matrix(descriptor)
  
  maxValues <- c()
  for(i in 1:nrow(workerSections)){
    startCol <- workerSections[i,1]
    endCol <- workerSections[i,2]
    maxValues <- cbind(maxValues , max(bm[,startCol:endCol], na.rm = TRUE))
  }
  
  return(max(maxValues))
}

.parallelAvg <- function(x, index = NULL){
  
  descx <- paste("x", toString(index-1), ".desc", sep = "")
  #descZones <- paste("zones", toString(index), ".desc", sep = "")
  x <- attach.big.matrix(descx)
  newDescx <- paste("x", toString(index), ".desc", sep = "")
  newx <- attach.big.matrix(newDescx)
  
  
  #write.table(chunkMatrix, "C:/Users/Carlos/Desktop/chunkMatrix", sep="\t")
  
  for(i in 1:nrow(workerSections)){
    startCol <- workerSections[i,1]
    endCol <- workerSections[i,2]
    
    if(startCol>1){tempStartCol <- startCol - 1} else {tempStartCol <- startCol}
    if(endCol<dim(x)[2]){tempEndCol <- endCol + 1} else {tempEndCol <- endCol}
    
    chunkMatrix <- x[,tempStartCol:tempEndCol]
    
    ####################################
    #Solution BaseLine
    
    # mval <- mean(chunkMatrix, na.rm = TRUE)
    # s1d <- function(s) unclass(stats::filter(s,c(0.5,0,0.5)))
    # pad <- rbind(mval,cbind(mval,chunkMatrix,mval),mval)
    # pad <- (t(apply(pad,1,s1d)) + apply(pad,2,s1d))/2
    # sm <- (pad[2:(nrow(chunkMatrix)+1),2:(ncol(chunkMatrix)+1)])
    
    ####################################
    #Solution 1
    # start.time <- Sys.time()
    mval <- mean(chunkMatrix, na.rm = TRUE)
    s1d <- function(s) {
      a <- complete.cases(s)
      if(sum(a, na.rm=TRUE)<3){return(s)}
      res <- unclass(stats::filter(s[a],c(0.5,0,0.5)))
      s[a] <- res
      return(s)
    }
    
    pad <- rbind(mval,cbind(mval,chunkMatrix,mval),mval) # Adiciona mval á volta da matrix inicial
    pad <- (t(apply(pad,1,s1d)) + apply(pad,2,s1d))/2 #appaly -> aplica s1d a cada linha da matrix pad
    sm <- (pad[2:(nrow(chunkMatrix)+1),2:(ncol(chunkMatrix)+1)])
    
    # end.time <- Sys.time()
    # time.taken <- end.time - start.time
    # time.taken
    
    ####################################
    #Solution 2
    # start.time <- Sys.time()
    # addresses <- expand.grid(x = 1:nrow(chunkMatrix), y = 1:ncol(chunkMatrix))
    # m2<-cbind(NA,rbind(NA,chunkMatrix,NA),NA)
    # naLocation <- is.na(chunkMatrix)
    # ret<-c()
    # for(i in 1:-1)
    #   for(j in 1:-1)
    #     if(i!=0 || j !=0)
    #       ret<-rbind(ret,m2[addresses$x+i+1+nrow(m2)*(addresses$y+j)])
    # 
    # sm <- matrix(colMeans(ret, na.rm = TRUE) , ncol = ncol(chunkMatrix), nrow = nrow(chunkMatrix))
    # write.table(ret, file="C:/Users/Carlos/Desktop/ret.dat",sep="\t", quote=FALSE)
    # sm[naLocation] <- NA
    # end.time <- Sys.time()
    # time.taken <- end.time - start.time
    # time.taken
    ##############################
    
    chunkMatrix <- chunkMatrix*r + (1-r)*sm
    #Note: The situation where (startCol==1 && endCol==dim(x)[2]) is not necessey as long as no_cores > 1
    if(startCol>1 && endCol<dim(x)[2]){
      newx[,startCol:endCol] <- chunkMatrix[,2:(dim(chunkMatrix)[2]-1)]
    }else if(startCol==1){
      newx[,startCol:endCol] <- chunkMatrix[,1:(dim(chunkMatrix)[2]-1)]
    }else if(endCol==dim(x)[2]){
      newx[,startCol:endCol] <- chunkMatrix[,2:dim(chunkMatrix)[2]]}
    rm(chunkMatrix)
    
  }
  return(TRUE)
}

.worker7 <- function(x, index = NULL, dbpath=NULL){
  
  descx <- paste("x", toString(index), ".desc", sep = "")
  #descZones <- paste("zones", toString(index), ".desc", sep = "")
  x <- attach.big.matrix(descx)
  
  sizeMb <- 200
  maxLengthInds <- (sizeMb*1024*1024)/8
  con <- dbConnect(RSQLite::SQLite(), dbname = dbpath)
  if(Sys.info()["sysname"][[1]] == "Windows"){
    #http://www.sqlite.org/faq.html#q5
    res <- dbSendQuery(con, "PRAGMA busy_timeout=10;")
    dbClearResult(res)
  }
  getIndicesQuery <- dbSendQuery(con, "SELECT x, y FROM indices WHERE idZona = :x")
  for (zone in zoneList) {
    sumZone <- 0
    #print(paste(workerId, "- Write Zone =", zone , "at", Sys.time(), sep = " "))
    
    lengthInds <- countsList[[as.character(zone)]]
    if(lengthInds == 0){
      next
    }
    boolPartition <- lengthInds >= maxLengthInds
    dbBind(getIndicesQuery, param = list(x = zone))
    
    if(!boolPartition){
      indexes <- as.matrix(dbFetch(getIndicesQuery))
      sumZone <- sum(x[indexes])
    }else {
      while(!dbHasCompleted(getIndicesQuery)){
        indexes <- as.matrix(dbFetch(getIndicesQuery, n = maxLengthInds))
        sumZone <- sumZone + sum(x[indexes])
      }
    }
    correct <- (pops[zone] - sumZone) / lengthInds
    
    if(!boolPartition){
      x[indexes] <- x[indexes] + correct
    }else {
      dbBind(getIndicesQuery, param = list(x = zone))
      while(!dbHasCompleted(getIndicesQuery)){
        indexes <- as.matrix(dbFetch(getIndicesQuery, n = maxLengthInds))
        x[indexes] <- x[indexes] + correct
      }
    }
  }
  dbClearResult(getIndicesQuery)
  dbDisconnect(con)
}

.zeroReplace <- function (x, index = NULL){
  
  desc <- paste("x", toString(index), ".desc", sep = "")
  bm <- attach.big.matrix(desc)
  for(i in 1:nrow(workerSections)){
    startCol <- workerSections[i,1]
    endCol <- workerSections[i,2]
    matrixChunk <- bm[,startCol:endCol]
    matrixChunk[ matrixChunk < 0 ] <- 0
    bm[,startCol:endCol] <- matrixChunk
  }
}

.worker8 <- function(x, index = NULL, dbpath=NULL){
  descx <- paste("x", toString(index), ".desc", sep = "")
  #descZones <- paste("zones", toString(index), ".desc", sep = "")
  x <- attach.big.matrix(descx)
  
  sizeMb <- 200
  maxLengthInds <- (sizeMb*1024*1024)/8
  con <- dbConnect(RSQLite::SQLite(), dbname = dbpath)
  if(Sys.info()["sysname"][[1]] == "Windows"){
    #http://www.sqlite.org/faq.html#q5
    res <- dbSendQuery(con, "PRAGMA busy_timeout=10;")
    dbClearResult(res)
  }
  getIndicesQuery <- dbSendQuery(con, "SELECT x, y FROM indices WHERE idZona = :x")
  for (zone in zoneList) {
    sumZone <- 0
    #print(paste(workerId, "- Write Zone =", zone , "at", Sys.time(), sep = " "))
    
    lengthInds <- countsList[[as.character(zone)]]
    if(lengthInds == 0){
      next
    }
    boolPartition <- lengthInds >= maxLengthInds
    dbBind(getIndicesQuery, param = list(x = zone))
    
    if(!boolPartition){
      indexes <- as.matrix(dbFetch(getIndicesQuery))
      sumZone <- sum(x[indexes])
    }else {
      while(!dbHasCompleted(getIndicesQuery)){
        indexes <- as.matrix(dbFetch(getIndicesQuery, n = maxLengthInds))
        sumZone <- sumZone + sum(x[indexes])
      }
    }
    correct <- pops[zone]/sumZone
    
    if(!boolPartition){
      x[indexes] <- x[indexes] * correct
    }else {
      dbBind(getIndicesQuery, param = list(x = zone))
      while(!dbHasCompleted(getIndicesQuery)){
        indexes <- as.matrix(dbFetch(getIndicesQuery, n = maxLengthInds))
        x[indexes] <- x[indexes] * correct
      }
    }
    rm(indexes)
  }
  dbClearResult(getIndicesQuery)
  dbDisconnect(con)
}

.maxAbsDiff <- function(x, index = NULL){
  
  descx <- paste("x", toString(index-1), ".desc", sep = "")
  x <- attach.big.matrix(descx)
  
  newDescx <- paste("x", toString(index), ".desc", sep = "")
  newx <- attach.big.matrix(newDescx)
  
  values <- c()
  for(i in 1:nrow(workerSections)){
    values <- c(values, workerSections[i,1])
    temp <- workerSections[i,1] + floor(length(workerSections[i,1]:workerSections[i,2])/2)
    values <- c(values, temp)
    values <- c(values, temp + 1)
    values <- c(values, workerSections[i,2])
  }
  halfWorkerSections <- matrix(values, ncol = 2, nrow = nrow(workerSections)*2, byrow = TRUE)
  
  maxValues <- c()
  for(i in 1:nrow(halfWorkerSections)){
    startCol <- halfWorkerSections[i,1]
    endCol <- halfWorkerSections[i,2]
    newMatrixChunk <- newx[,startCol:endCol]
    oldMatrixChunk <- x[,startCol:endCol]
    maxValues <- c(maxValues, max(abs(oldMatrixChunk - newMatrixChunk), na.rm = TRUE))
  }
  
  return(max(maxValues))
  
}

.rasterResult <- function (x, index = NULL){
  rasterOptions(todisk = TRUE, tmpdir = getwd())
  desc <- paste("x", toString(index), ".desc", sep = "")
  bm <- attach.big.matrix(desc)
  
  desc <- paste(workerId, "Coords", ".desc", sep = "")
  coord1 <- attach.big.matrix(desc)
  
  sizeMb <- 200
  nMaxcellsToLoad <- (sizeMb*1024*1024)/28
  nColsToLoad <- floor(nMaxcellsToLoad/(dim(bm)[1]))
  
  sections <- matrix(nrow = 0 ,ncol = 2)
  startCol <- 1
  endCol <- nColsToLoad
  sections <- rbind(sections, c(startCol,endCol))
  #print(paste("startCol-", startCol, "endCol-", endCol))
  k <- 0
  while(endCol < dim(coord1)[2]){
    startCol <- endCol + 1
    endCol <- endCol + nColsToLoad
    if(endCol > dim(coord1)[2]){
      endCol <- dim(coord1)[2]
    }
    sections <- rbind(sections, c(startCol,endCol))
    k <- k + 1
  }
  
  if(endCol > dim(coord1)[2]){sections[1,2] <- dim(coord1)[2]}
  
  #total <- rbind(t(matrix(coord1[,],nrow=2)), t(matrix(coord2[,],nrow=2)), t(matrix(coord3[,],nrow=2)), t(matrix(coord4[,],nrow=2)))
  rasterList <- list()
  for(i in 1:nrow(sections)){
    coordStart <- sections[i,][1]
    coordEnd <- sections[i,][2]
    bmStart <- colSpanStart + coordStart - 1
    bmEnd <- colSpanStart + coordEnd - 1
    rasterList[[i]] <- raster(SpatialPixelsDataFrame(t(matrix(coord1[,coordStart:coordEnd],nrow=2)) , data.frame(dens=c(array(bm[,bmStart:bmEnd])))))
  }
  if(length(rasterList)>1){
    rasters1.mosaicargs <- rasterList
    rasters1.mosaicargs$fun <- mean
    result <- do.call(mosaic, rasters1.mosaicargs)
  } else {result <- rasterList[[1]]}
  
  return(result)
  
}

#######################################################################################
#                              Pycno main function                                    #
#######################################################################################

# Pycnophylactic interpolation, adapted from the pycno package by Chris Brunsdon.
# Given a SpatialPolygonsDataFrame and a set of populations for each polygon, compute a population density estimate based on Tobler's pycnophylactic interpolation algorithm. The result is a SpatialGridDataFrame.
.pycno <- function( x, pops, celldim, r=0.2, converge=3, verbose=TRUE ) {
  no_cores <- length(cl)
  bbx <- slot(x,'bbox')
  offset <- bbx[,1]
  extent <- bbx[,2] - offset
  globalShape <- ceiling(extent / celldim)
  intervaly <- floor(globalShape/no_cores)[[2]]
  subSections <- rep(intervaly, no_cores)
  
  i <- 1
  while(sum(subSections) < globalShape[[2]]){
    subSections[i] <- subSections[i] + 1
    i <- i + 1
  }
  
  clusterExport(cl, "globalShape", envir = environment())
  clusterExport(cl, "celldim", envir = environment())
  for(i in 1:no_cores)  {
    if (!is(celldim,"SpatialGrid")) {
      if( i == 1){
        offset <- bbx[,1]
      }else{
        offset[2] <- offset[2] + (celldim * subSections[i-1])
      }			
      clusterExport(cl[no_cores - (i-1)], "offset", envir = environment())
      subSection <- subSections[i]
      clusterExport(cl[no_cores - (i-1)], "subSection", envir = environment())
      #gridList[[i]] <- SpatialGrid(GridTopology(offset,c(celldim,celldim),c(globalShape[1], subSections[i])))
      
    } else {
      #TODO:
      gr <- celldim
      gridList[[i]] <- gr
    }
  }
  px <- CRS(proj4string(x))
  
  clusterExport(cl, "px", envir = environment())
  poligonDataFrame <- x
  clusterExport(cl, "poligonDataFrame", envir = environment())
  colSpanStart <- 1
  for(i in 1:no_cores){
    workerId <- i
    clusterExport(cl[i], "workerId", envir = environment())
    #gridElement <- gridList[[no_cores - (i-1)]]
    colSpanEnd <- subSections[no_cores - (i-1)] + colSpanStart - 1
    #clusterExport(cl[i], "gridElement", envir = environment())
    clusterExport(cl[i], "colSpanStart", envir = environment())
    #clusterExport(cl[i], "colSpanEnd", envir = environment())
    colSpanStart <- colSpanEnd + 1
  }
  
  clusterEvalQ(cl, library(rgdal))
  clusterEvalQ(cl, library(bigmemory))
  clusterEvalQ(cl, library(DBI))
  clusterEvalQ(cl, library(flock))
  clusterEvalQ(cl, library(RSQLite))
  
  zone.list <- sapply(slot(x, "polygons"), function(x) as.numeric(slot(x, "ID")))
  zone.list <- tail(zone.list, length(zone.list)-1) #removes first ID zero
  zone.list <- c(zone.list, max(zone.list) + 1, max(zone.list) + 2)
  naId <- max(zone.list)
  
  clusterExport(cl, "naId", envir = environment())
  clusterExport(cl, "zone.list", envir = environment())
  
  lockname = tempfile(tmpdir = getwd())
  dbpath <- tempfile(tmpdir = getwd())
  con <- dbConnect(RSQLite::SQLite(), dbname=dbpath)
  df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df) <- c("idZona", "x", "y")
  dbWriteTable(con, "indices", df)
  qry1 <- dbSendQuery(con,"PRAGMA synchronous=OFF")
  dbClearResult(qry1)
  qry2 <- dbSendQuery(con,"PRAGMA count_changes=OFF")
  dbClearResult(qry2)
  qry3 <- dbSendQuery(con,"PRAGMA journal_mode=MEMORY")
  dbClearResult(qry3)
  qry4 <- dbSendQuery(con,"PRAGMA temp_store=MEMORY")
  dbClearResult(qry4)
  #qry5 <- dbSendQuery(con,"PRAGMA read_uncommitted=TRUE")
  #dbClearResult(qry5)
  
  
  print(paste("Start overlap at", Sys.time(),sep = " "))
  result <- parLapply(cl, 1:no_cores, .overlap, lockname=lockname, dbpath=dbpath)
  print(paste("End overlap at", Sys.time(),sep = " "))
  qry <- dbSendQuery(con,"CREATE INDEX index_zona ON indices (idZona)")
  dbClearResult(qry)
  #dbDisconnect(con)
  
  pops <- c(pops,0)
  j <- 1
  backx <- paste("x", toString(j), ".bin", sep = "")
  descx <- paste("x", toString(j), ".desc", sep = "")
  
  x <- big.matrix(globalShape[[1]], globalShape[[2]], type = "double",  init = 0, separated = FALSE, backingfile = backx, descriptorfile = descx)
  
  ################################
  # TODO: Give equal size zoneIndexes subsets no all workers
  # aux2 <- list()
  # sizeResultsMerged <- object.size(resultsMerged)
  # i <- 1
  # end <- length(resultsMerged)
  # 
  # while( object.size(aux2) < sizeResultsMerged/no_cores){
  #   if(i == end + 1){break()}
  #   aux2[names(resultsMerged[i])] <- resultsMerged[i]
  #   i <- i + 1
  # }
  #################################
  
  splitsZoneList <- split(zone.list, sort(zone.list%%no_cores))
  clusterApply(cl, splitsZoneList, function(m) {zoneList <<- m; NULL})
  #clusterApply(cl, resultsMerged, function(m) {zonesIndex <<- m; NULL})
  clusterExport(cl, "pops", envir = environment())
  
  # for(i in 1:no_cores){
  #   zonesIndex <- resultsMerged[[i]]
  #   clusterExport(cl[i], "zonesIndex", envir = environment())
  #   clusterExport(cl[i], "pops", envir = environment())
  # }
  
  print(paste("Start worker5 at", Sys.time(),sep = " "))
  noresult <-parLapply(cl, 1:no_cores, .worker5, index = j, dbpath=dbpath)
  print(paste("End worker5 at", Sys.time(),sep = " "))
  #defines the region on wich worker nodes will search for a max value
  
  
  #Start Partitioning Strategy
  sizeMb <- 200
  nMaxcellsToLoad <- (sizeMb*1024*1024)/8
  nColsToLoad <- floor(nMaxcellsToLoad/(dim(x)[1]))
  nRowsToLoad <- dim(x)[1]
  
  sections <- matrix(nrow = 0 ,ncol = 2)
  startCol <- 1
  endCol <- nColsToLoad
  sections <- rbind(sections, c(startCol,endCol))
  #print(paste("startCol-", startCol, "endCol-", endCol))
  k <- 0
  while(endCol < dim(x)[2]){
    startCol <- endCol + 1
    endCol <- endCol + nColsToLoad
    if(endCol > dim(x)[2]){
      endCol <- dim(x)[2]
    }
    sections <- rbind(sections, c(startCol,endCol))
    k <- k + 1
  }
  
  if(nrow(sections)<no_cores){
    sections <- matrix(nrow = 0 ,ncol = 2)
    startCol <- 1
    endCol <- subSections[1]
    sections <- rbind(sections, c(startCol,endCol))
    for(i in 2:no_cores){
      startCol <- endCol + 1
      endCol <- endCol + subSections[i]
      sections <- rbind(sections, c(startCol,endCol))
      #print(paste("startCol-", startCol, "endCol-", endCol))
    }
  }
  
  colsToLoad <- rep( floor(nrow(sections)/no_cores), no_cores)
  i <- 1
  while(sum(colsToLoad) < nrow(sections)){
    colsToLoad[i] <- colsToLoad[i] + 1
    i <- i + 1
  }
  
  tempRow <- 1
  for(i in 1:no_cores){
    workerSections <- sections[tempRow:(tempRow + colsToLoad[i] - 1),]
    tempRow <- tempRow + colsToLoad[i]
    workerSections<- matrix(c(workerSections), ncol = 2)
    clusterExport(cl[i] , "workerSections", envir = environment())
  }
  
  #End Partitioning Strategy
  print(paste("Start maxBigMatrix at", Sys.time(),sep = " "))
  maxValues <-parLapply(cl, 1:no_cores, .maxBigMatrix, descriptor = descx)
  print(paste("End maxBigMatrix at", Sys.time(),sep = " "))
  #TODO: tornar operação de encontrar maximo paralela, por forma a poder carregar em memoria segmentos que no total prefaçam tamanho < RAM
  #x <- as.matrix(x)
  stopper <- max(unlist(maxValues)) * 10^(-converge)
  
  
  clusterExport(cl , "r", envir = environment())
  
  repeat {
    j <- j + 1
    old.x <- x
    back <- paste("x", toString(j), ".bin", sep = "")
    desc <- paste("x", toString(j), ".desc", sep = "")
    x <- big.matrix(globalShape[[1]], globalShape[[2]], type = "double",  init = NA, separated = FALSE, backingfile = back, descriptorfile = desc)
    
    print(paste("Start parallelAvg at", Sys.time(),sep = ""))
    parLapply(cl, 1:no_cores, .parallelAvg, index = j)
    print(paste("End parallelAvg at", Sys.time(),sep = ""))
    
    # mval <- mean(x)
    # s1d <- function(s) unclass(stats::filter(s,c(0.5,0,0.5)))
    # pad <- rbind(mval,cbind(mval,x,mval),mval)
    # pad <- (t(apply(pad,1,s1d)) + apply(pad,2,s1d))/2
    # sm <- (pad[2:(nrow(x)+1),2:(ncol(x)+1)])
    #x <- x*r + (1-r)*sm
    
    print(paste("Start worker7 at", Sys.time(),sep = ""))
    parLapply(cl, 1:no_cores, .worker7, index = j, dbpath=dbpath)
    print(paste("End worker7 at", Sys.time(),sep = ""))
    
    print(paste("Start zeroReplace at", Sys.time(),sep = ""))
    parLapply(cl, 1:no_cores, .zeroReplace, index = j)
    print(paste("End zeroReplace at", Sys.time(),sep = ""))
    
    print(paste("Start worker8 at", Sys.time(),sep = ""))
    parLapply(cl, 1:no_cores, .worker8, index = j, dbpath=dbpath)
    print(paste("End worker8 at", Sys.time(),sep = ""))
    
    print(paste("Start maxAbsDiff at", Sys.time(),sep = ""))
    maxAbsDiffs <-parLapply(cl, 1:no_cores, .maxAbsDiff, index = j)
    print(paste("End maxAbsDiff at", Sys.time(),sep = ""))
    
    stopperTest <- max(unlist(maxAbsDiffs))
    
    if (verbose) {
      flush.console()
      cat(sprintf("Maximum Change: %12.5f - will stop at %12.5f\n", stopperTest, stopper))
    }
    
    if (stopperTest <= stopper) break
  }
  
  #if (!is.null(attr(x,'na'))) x[attr(x,'na')] <- NA
  clusterEvalQ(cl, library(raster))
  print(paste("Start rasterResult at", Sys.time(),sep = " "))
  rasterListFinal <-parLapply(cl, 1:no_cores, .rasterResult, index = j)
  print(paste("End rasterResult at", Sys.time(),sep = " "))
  
  rasters1.mosaicargs <- rasterListFinal
  rasters1.mosaicargs$fun <- mean
  rasterFinal <- do.call(mosaic, rasters1.mosaicargs)
  
  #coordList <- parLapply(cl, 1:no_cores, .workerGetCoords)
  
  #coordListMerged <- rbindlist(coordList)
  
  #result <- SpatialPixelsDataFrame( coordListMerged, data.frame(dens=array(x)) )
  #result <- as( result , "SpatialGridDataFrame" )
  proj4string(rasterFinal) <- px
  dbDisconnect(con)
  print("END")
  return(rasterFinal)
}

#######################################################################################
#                       Methods regarding secondary memory                            #
#######################################################################################

.mass_preserv_ffdf <- function(diss_result){
  
  summary_metric <- function( data, type ) {
    if (type == 'count') { return( sum(data) ) }
    if (type == 'categorical') { return( median(data) ) }
    return( mean(data) )
  }
  
  splitby <- as.character.ff(diss_result$cell)
  
  grp_qty <- ffdfdply(x=diss_result[c("cell","diss", nm_coarse)], 
                      split=splitby, 
                      FUN = function(data){
                        ## This happens in RAM - containing **several** split elements so here we can use data.table which works fine for in RAM computing
                        require(data.table)
                        data <- as.data.table(data)
                        result <- data[, list(diss = summary_metric(diss,data_type), count = length(diss), layer=first(layer)), by = list(cell)]
                        as.data.frame(result)
                      })
  
  
  
  #Mass preservation
  #temp_inds2 <- ffwhich(grp_qty, diss==0)
  #if(!is.null(temp_inds)){grp_qty$diss[temp_inds2] <- grp_qty$layer[temp_inds2]}
  
  #calculate mass preserving corrections1
  corrections <- ff(vmode="double", length = nrow(grp_qty))
  cursor <- 1
  parts <- chunk(grp_qty)
  for(part in parts){
    temp <- ((grp_qty[part,,drop=FALSE]["layer"]-grp_qty[part,,drop=FALSE]["diss"]) / grp_qty[part,,drop=FALSE]["count"])[,]
    write.ff(corrections, value=as.vector(temp), i=cursor, add = FALSE)
    cursor <- cursor + length(temp)
  }
  
  grp_qty$corrections <- corrections
  
  #Aplay corrections1 to 
  diss_result <- merge(diss_result,grp_qty[c("cell","corrections")], by="cell", all.x=TRUE, trace = TRUE)
  parts <- chunk(diss_result)
  for(part in parts){
    diss_result$diss[part] <- diss_result$diss[part] + diss_result$corrections[part]
  }
  
  diss_result$corrections <- NULL
  
  temp_inds <- ffwhich(diss_result, diss<0)
  if(!is.null(temp_inds)){diss_result$diss[temp_inds] <- ff(0, length = length(temp_inds))}
  
  grp_qty2 <- ffdfdply(x=diss_result[c("cell","diss", nm_coarse)], 
                       split=splitby, 
                       FUN = function(data){
                         ## This happens in RAM - containing **several** split elements so here we can use data.table which works fine for in RAM computing
                         require(data.table)
                         data <- as.data.table(data)
                         result <- data[, list(diss = summary_metric(diss,data_type), count = length(diss), layer=first(layer)), by = list(cell)]
                         as.data.frame(result)
                       })
  
  #calculate mass preserving corrections2
  corrections2 <- ff(vmode="double", length = nrow(grp_qty))
  cursor <- 1
  parts <- chunk(grp_qty2)
  for(part in parts){
    temp <- (grp_qty2[part,,drop=FALSE]["layer"]/grp_qty2[part,,drop=FALSE]["diss"])[,]
    write.ff(corrections2, value=as.vector(temp), i=cursor, add = FALSE)
    cursor <- cursor + length(temp)
  }
  grp_qty2$corrections2 <- corrections2
  
  
  #Aplay corrections2 to results
  diss_result <- merge(diss_result,grp_qty2[c("cell","corrections2")], by="cell", all.x=TRUE, trace = TRUE)
  parts <- chunk(diss_result)
  for(part in parts){
    diss_result$diss[part] <- diss_result$diss[part] * diss_result$corrections[part]
  }
  
  diss_result$corrections2 <- NULL
  
  return(diss_result)
}

.ffdf_to_raster <- function(diss_result, map, res, crs){
  
  parts <- chunk(diss_result)
  rasterOptions(todisk = TRUE)
  
  rasterList <- foreach( i = 1:length(parts), .packages = c('ff','ffbase', 'raster')) %dopar% {
    rasterOptions(todisk = TRUE)
    #writeRaster(rasterFromXYZ( data.frame( diss_result[parts[[i]],,drop=FALSE][,c('x', 'y')], diss = map[parts[[i]]]), res = res, crs = crs), filename = paste("raster", parts[[i]][[1]], "to", parts[[i]][[2]], sep = ""), overwrite=TRUE)
    rasterFromXYZ( data.frame( diss_result[parts[[i]],,drop=FALSE][,c('x', 'y')], diss = map[parts[[i]]]), res = res, crs = crs)
    #teste <- rasterFromXYZ( data.frame( as.ram(diss_result[parts[[i]],,drop=FALSE][c('x', 'y')]), diss = map ,res = res, crs = crs))
  }
  
  # rows <- 0
  # maxcols <- 0
  # rasterList <- list()
  # for ( i in 1:length(parts)){
  #   rasterList[[i]] <- raster(paste("raster", parts[[i]][[1]], "to", parts[[i]][[2]], sep = ""))
  #   rows <- rows + dim(rasterList[[i]])[1]
  #   maxcols <- max(maxcols, dim(rasterList[[i]])[2])
  # }
  
  if(length(rasterList)>1){
    #rasters1.mosaicargs <- rasterList
    rasterList$fun <- mean
    rasterList$filename <- "final3.grd"
    rasterList$overwrite <- TRUE
    result <- do.call(mosaic, rasterList)
  } else {result <- rasterList[[1]]}
  
  removeTmpFiles()
  
  # tempRast <- rasterList[[1]]
  # if(length(rasterList)>1){
  #   for(i in 2:length(rasterList)){
  #     tempRast <- merge(tempRast,  rasterList[[i]], filename=paste("final", i ,".grd", sep = ""), overwrite=TRUE)
  #     #tempRast <- merge(tempRast,  rasterList[[i]], filename="final4.grd", overwrite=TRUE)
  #     if (i > 2){
  #       file.remove(paste("final", i-1 ,".grd", sep = ""),paste("final", i-1 ,".gri", sep = ""))
  #     }
  #     
  #   }
  # }
  
  # tempRast <- raster(system.file("external/test.grd", package="raster"))
  # extent(tempRast) <- c(min(diss_result$x), max(diss_result$x), min(diss_result$y), max(diss_result$y))
  # dim(tempRast) <- c(rows, maxcols)
  # res(tempRast) <- res(rasterList[[1]])
  # crs(tempRast) <- crs(rasterList[[1]])
  # tempRast <- writeStart(teste, "final.grd", overwrite=TRUE)
  # #rowCursor <- nrow(tempRast)
  # rowCursor <- 1
  # if(length(rasterList)>1){
  #   for(i in 1:length(rasterList)){
  #     #tempRast <- merge(tempRast,  rasterList[[i]], filename="final.grd", overwrite=TRUE)
  #     v <- getValuesBlock(rasterList[[i]], row=1, nrows=nrow(rasterList[[i]]) )
  #     tempRast <- writeValues(tempRast, v, rowCursor)
  #     rowCursor <- rowCursor + nrow(rasterList[[i]])
  #   }
  # }
  # tempRast <- writeStop(tempRast)
  
  return(result)
}

.raster_to_ff_data_frame_xy <- function(rstr, ...) {
  bs <- blockSize(rstr)
  vars <- list()
  for(i in 1:length(names(rstr))){
    vars[[names(rstr)[i]]] <- ff(vmode="double", length = rstr@ncols*rstr@nrows)
  }
  #data <- ff(vmode="double", length = rstr@ncols*rstr@nrows)
  y <- ff(vmode="double", length = rstr@ncols*rstr@nrows)
  x <- ff(vmode="double", length = rstr@ncols*rstr@nrows)
  cursor <- 1
  for (i in 1:bs$n) {
    #v <- raster(matrix(getValues(rstr, row=bs$row[i], nrows=bs$nrows[i]) , ncol = dim(rstr)[2] , byrow = TRUE))
    tempData <- getValues(rstr, row=bs$row[i], nrows=bs$nrows[i])
    if(is(tempData,"matrix")){
      for(j in 1:length(colnames(tempData))){
        write.ff(vars[[colnames(tempData)[j]]], value=tempData[,colnames(tempData)[j]], i=cursor, add = FALSE)
      }
    } else if (is(tempData,"numeric")){
      write.ff(vars[[names(rstr)]], value=tempData, i=cursor, add = FALSE)
    }
    rm(tempData)
    tempY <- rep(yFromRow(rstr, bs$row[i]:(bs$row[i] + bs$nrows[i] - 1)), each = dim(rstr)[2])
    write.ff(y, value=tempY, i=cursor, add = FALSE)
    rm(tempY)
    tempX <- rep(xFromCol(rstr, 1:dim(rstr)[2]), times = bs$nrows[i])
    write.ff(x, value=tempX, i=cursor, add = FALSE)
    cursor <- cursor + length(tempX)
    rm(tempX)
  }
  raster_ffdf <- ffdf(x=x, y=y)
  for(i in 1:length(names(rstr))){
    raster_ffdf[[names(rstr)[i]]] <- vars[[names(rstr)[i]]]
  }
  if (any(is.factor(x))) {
    #TODO: ?
  }
  return(raster_ffdf)
}

.raster_values_to_ff <- function(rstr, ...) {
  bs <- blockSize(rstr)
  cell <- ff(vmode="double", length = rstr@ncols*rstr@nrows)
  cursor <- 1
  for (i in 1:bs$n) {
    #v <- raster(matrix(getValues(rstr, row=bs$row[i], nrows=bs$nrows[i]) , ncol = dim(rstr)[2] , byrow = TRUE))
    tempCells <- as.integer(getValues(rstr, row=bs$row[i], nrows=bs$nrows[i]))
    write.ff(cell, value=tempCells, i=cursor, add = FALSE)
    cursor <- cursor + length(tempCells)
    rm(tempCells)
  }
  if (any(is.factor(rstr))) {
    #TODO: ?
  }
  return(cell)
}

.create_lut_fine_ff <- function(coarse, fine) {
  bs <- blockSize(fine)
  coarseVals <- ff(vmode="integer", length = fine@ncols*fine@nrows)
  cursor <- 1
  for (i in 1:bs$n) {
    print(paste("Witing to disk block", i, "of", bs$n, sep = " "))
    temp_df <- data.frame(x=rep(xFromCol(fine, 1:dim(fine)[2]), times = bs$nrows[i]), 
                          y=rep(yFromRow(fine, bs$row[i]:(bs$row[i] + bs$nrows[i] - 1)), each = dim(fine)[2]))
    
    tempVals <- as.integer(extract(coarse, temp_df))
    rm(temp_df)
    write.ff(coarseVals, value=tempVals, i=cursor, add = FALSE)
    cursor <- cursor + length(tempVals)
    rm(tempVals)
  }
  return(coarseVals)
}

.create_ff_indices <- function(rowsLenght) {
  coarseVals <- ff(vmode="integer", length = rowsLenght)
  for (i in chunk(coarseVals)){
    message("Writing chunk - ", i[[1]],":",i[[2]])
    coarseVals[i] <- i[[1]]:i[[2]]
  }
  return(coarseVals)
}

.bigsample_in_chunks <- function(ff_inds, size){
  chunks <- chunk(ff_inds)
  id_spl <- ff(vmode="integer", length = size)
  
  interval <- ceiling(size/length(chunks))
  sections <- rep(interval, length(chunks))
  
  correction <- sum(sections) - size
  sections[length(sections)] <- sections[length(sections)] - correction
  # i <- 1
  # while(sum(sections) < size){
  #   sections[i] <- sections[i] + 1
  #   i <- i + 1
  # }
  
  i <- 1
  cursor <- 1
  for (chunk in chunks){
    message("Writing chunk - ", chunk[[1]],":",chunk[[2]])
    temp <- bigsample(ff_inds[chunk], size = sections[i])
    write.ff(id_spl, value=temp, i=cursor, add = FALSE)
    cursor <- cursor + length(temp)
    i <- i + 1
  }
  return(id_spl)
}

.na_exclude_ffdf <- function(fine_df){
  nms <- names(fine_df)
  nms <- nms[!nms %in% c("x","y")]
  for (nm in nms){
    inds <- !is.na(fine_df[[nm]])
    fine_df <- fine_df[inds]
    rm(inds)
  }
  return(fine_df)
}

.worker <- function(x){
  print(paste ("Split", x, "start", Sys.time(), sep = " "))
  part <- rasterize(poligonDataFrame[parts[[x]],], raster( resolution=minres * 1.01, ext=extent(poligonDataFrame) ), coarse_var_name, fun='first')
  print(paste ("Split", x, "end", Sys.time(), sep = " "))
  return(part)
}

.parallel_rasterize <- function(coarse, coarse_var_name, minres){
  no_cores <- length(cl)
  # Number of polygons features in SPDF
  features <- 1:nrow(coarse[,])
  # Split features in n parts
  n <- no_cores
  parts <- split(features, cut(features, n))
  clusterEvalQ(cl, library(raster))
  #TODO: coarse is exported to the cluster inside pycno method, export in dissever method instead
  #clusterExport(cl, "coarse", envir=environment())
  if(!all(unlist(clusterEvalQ(cl, exists("poligonDataFrame"))))){
    poligonDataFrame <- coarse
    clusterExport(cl, "poligonDataFrame", envir=environment())
  }
  clusterExport(cl, "parts", envir=environment())
  clusterExport(cl, "minres", envir=environment())
  clusterExport(cl, "coarse_var_name", envir=environment())
  rParts <- parLapply(cl, 1:no_cores, .worker)
  ids_coarse <- do.call(merge, rParts)
  rm(rParts)
  return(ids_coarse)
}

.join_interpol_ffdf <- function(coarse_df, fine_df, attr, by) {
  # Nearest-neighbour interpolation as an inner SQL join
  left_join_ffdf <- merge(fine_df, coarse_df, by = by, all.x=TRUE, all.y=FALSE, trace = TRUE) #Left join
  cols_to_remove<- names(left_join_ffdf)[!grepl(paste0(attr, collapse = "|"), names(left_join_ffdf))]
  return(left_join_ffdf[setdiff(colnames(left_join_ffdf), cols_to_remove)])
}

.cbind_ffdf2 <- function(d1, d2){
  D1names <- colnames(d1)
  D2names <- colnames(d2)
  mergeCall <- do.call("ffdf", c(physical(d1), physical(d2)))
  colnames(mergeCall) <- c(D1names, D2names)
  mergeCall
}

#NOTA: alteraçao dos agrumentos de chamada face ao metofo nao ffdf
.update_model_ffdf <- function(vars, y, method = 'rf', control, tune_grid, data_type="numeric", long=NULL , lat=NULL) {
  # Pick the parameters of the model using error on the first run.
  # Then use the optimised parameters in the iteration loop to save on computing time.
  # Basically we just need to change the trainControl object to do that.
  if ( method == 'lm' ){
    form <- as.formula(paste("x~",paste(names(vars), collapse="+")))
    temp_data <- .cbind_ffdf2(vars, ffdf(x=y))
    for (i in chunk(temp_data)){
      if (i[1]==1){
        message("first chunk is: ", i[[1]],":",i[[2]])
        fit <- biglm(form, data=temp_data[i,,drop=FALSE])
      }else{
        message("next chunk is: ", i[[1]],":",i[[2]])
        fit <- update(fit, temp_data[i,,drop=FALSE])
      }
    }
    return(fit)
  } else {
    #TODO: temp fix, use sample strategy?
    vars = as.ram(fine_df[id_spl, nm_covariates])
    latLong = data.frame( long=as.ram(long), lat = as.ram(lat))
    y <- as.ram(y)
  }
  
  y_aux = y
  if ( data_type == 'categorical' ) { y_aux = factor( y_aux ) }
  if ( method == 'gwr' ) { 
    fit <- gw( as.formula(paste("x~",paste(names(vars), collapse="+"))) , data= data.frame( vars , x=y_aux ) )
  } else if ( method == 'mlp' ) {
    data <- mx.symbol.Variable("data")
    fc1 <- mx.symbol.FullyConnected(data, num_hidden=10)
    act1 <- mx.symbol.Activation(fc1, act_type="relu")
    fc2 <- mx.symbol.FullyConnected(act1, num_hidden=1)
    fit <- mx.model.FeedForward.create(mx.symbol.LinearRegressionOutput(fc2), X=vars, y=y_aux, num.round=50, array.batch.size=20, learning.rate=2e-6, momentum=0.9, eval.metric=mx.metric.rmse)
  } else if ( method == 'lme' ) {
    fit <- data.frame( vars , out=y_aux , lat=latLong$lat , long=latLong$long , dummy=rep.int( 1 , length(y_aux) ) )
    fit <- lme( fixed=out ~ . - dummy - lat - long , data=fit , random = ~ 1 | dummy, correlation = corGaus(form = ~ lat+long | dummy ) )
  } else {
    fit <- train( x = vars, y = y_aux, method = method, trControl = control, tuneGrid  = tune_grid )}
  fit
}

.predict_map_ffdf <- function(fit, data, boot = NULL, level = 0.9, data_type="numeric", long=NULL , lat=NULL) {
  
  #parts <- chunk(data, BATCHBYTES = 300000)
  parts <- chunk(data) #Uses default batch size
  res <- foreach( i = 1:length(parts), .packages = c('caret', 'ff','ffbase')) %dopar% {
    if(inherits(fit,"biglm")){
      library(biglm)
    }
    if (is.null(boot)) {
      if (!is.null(long) && !is.null(long) ) {
        latLong = data.frame( long=long[parts[[i]]], lat = lat[parts[[i]]])
        worker_res <- predict(object=fit , newdata = data.frame(data[parts[[i]],,drop=FALSE],latLong))
      } else {
        worker_res <- predict(object=fit , newdata = data[parts[[i]],,drop=FALSE])
      }
      if ( !is.null( nrow(worker_res) ) ) worker_res <- worker_res[,1]
      ff(as.numeric( worker_res ),finalizer = "close")
    } else {
      #TODO:
      #.bootstrap_ci(fit = fit, fine_df = data[parts[[i]],,drop=FALSE], level = level, n = boot, data_type=data_type, latLong=latLong)
    }
  }
  
  merge_res <- res[[1]]
  if(length(res) > 1){
    for(i in 2:length(res)){
      merge_res <- c(merge_res, res[[i]])
    }
  }
  return(merge_res)
}

#######################################################################################
#                              Dissever main method                                   #
#######################################################################################

.dissever <- function(
  coarse,
  fine,
  coarse_var_names = NULL,
  method = "lm",
  p = NULL,
  sample_method = 'random',
  nmax = NULL,
  thresh = 0.01,
  min_iter = 5,
  max_iter = 100,
  boot = NULL,
  level = 0.9,
  tune_length = 3,
  tune_grid = .create_tune_grid(model = method, tune_length = tune_length),
  train_control_init = .default_control_init,
  train_control_iter = .default_control_iter,
  data_type = "count",
  add_pycno = 0,
  no_cores = 4,
  verbose = FALSE
) {
  
  if(is.null(no_cores)){no_cores = detectCores()}
  if(no_cores==1){no_cores = 2}
  data_type <<- data_type
  cl <<- makePSOCKcluster(no_cores, methods = TRUE, outfile="cluster.log")
  setDefaultCluster(cl=cl)
  assign("cl",cl,.GlobalEnv)
  registerDoParallel(cl)
  
  input_polygons = class(coarse) == "SpatialPolygonsDataFrame"
  if ( data_type != "count" && add_pycno > 0 ) {
    stop('Initialization based on pycnophylactic interpolation should only be used with count data')
  }
  if ( !input_polygons && class(coarse) != "RasterLayer" ) {
    stop('The course data should be provided as a SpatialPolygonsDataFrame or as a RasterLayer')
  }
  if ( !input_polygons && !is.null(coarse_var_names) ) {
    stop('The parameter coarse_var_names should only be used when providing a SpatialPolygonsDataFrame as the coarse data')
  }
  # Horrible hack, avoid division by 0 in pycnophylactic interpolation
  if(input_polygons && data_type == "count" && nrow(coarse[which(coarse[[coarse_var_names[2]]] == 0),]) > 0) {
    coarse[[coarse_var_names[2]]] = coarse[[coarse_var_names[2]]] + 0.00001
  }
  if ( input_polygons ) {
    if ( is.null(coarse_var_names) ) { coarse_var_names <- names( coarse ) }
    if ( length(coarse_var_names) > 2 ) {
      stop('The parameter coarse_var_names should be used to provide the names for attributes corresponding to the IDs of polygons and the quantity to be downscaled')
    }
    minres <- min(res(fine))
    if ( add_pycno > 0 ) { pycnolayer <- .pycno( coarse, coarse[[coarse_var_names[2]]], min(minres), converge=add_pycno, verbose=FALSE ) }
    else if ( data_type == "count" ) { pycnolayer <- .pycno( coarse, coarse[[coarse_var_names[2]]], min(minres), converge=0, verbose=FALSE ) }    
    #ids_coarse <- rasterize(coarse, raster( resolution=minres * 1.01, ext=extent(coarse) ), coarse_var_names[1], fun='first')
    ids_coarse <- .parallel_rasterize(coarse, coarse_var_names[1],minres)
    names(ids_coarse) <- 'cell'
    #coarse <- rasterize(coarse, raster( resolution=minres * 1.01, ext=extent(coarse) ), coarse_var_names[2], fun='first')
    coarse <- .parallel_rasterize(coarse, coarse_var_names[2], minres)
  } else if ( add_pycno > 0 ) {
    minres <- min(res(fine))
    pycnolayer <- .pycno( rasterToPolygons(coarse), .as_data_frame_factors(coarse), 0.05, converge=add_pycno, verbose=FALSE )
  }
  if (min(res(fine)) >= min(res(coarse))) { stop('Resolution of fine data should be higher than resolution of coarse data') }
  if (!(data_type == "numeric" || data_type == "count" || data_type == "categorical" )) {
    stop('Data type should be numeric, categorical or count')
  }
  
  library(ff)
  library(ffbase)
  # Store names of coarse data and fine-scale covariates
  nm_coarse <<- names(coarse)
  nm_covariates <<- names(fine)
  
  if(length(nm_coarse)>1){
    stop('Multi corse variables not suported in this version')
  } else{
    nm_coarse <- "layer"
  }
  
  # Get cell numbers of the coarse grid and convert coarse data to data.frame
  if ( !input_polygons ) {
    ids_coarse <- raster(coarse)
    ids_coarse[] <- 1:ncell(coarse)
    names(ids_coarse) <- 'cell'
    coarse_df <- .as_data_frame_factors(coarse, xy = TRUE)
    coarse_df$cell <- 1:nrow(coarse_df)
    coarse_df$cell2 <- coarse_df$cell
  } else {
    coarse_df <- .raster_to_ff_data_frame_xy(coarse)
    coarse_df$cell <- .raster_values_to_ff(ids_coarse)
    
    #TODO: No caso de os dados serem categoricos?
    #coarse_df$cell <- sapply(coarse_df$cell, function(x) if(is.factor(x)) { as.numeric(x) } else { x })
    #result <- ffvecapply(teste[i1:i2]+6, X=teste, RETURN=TRUE, BATCHSIZE=500000, VERBOSE=TRUE)
    #result <- ffvecapply(teste[i1:i2]<- if(teste[i2-i1 + 1]>5) { 7 } else { 8 } , X=teste, BATCHSIZE=500000, VERBOSE=TRUE)
    #result <- ffapply(X=teste,  AFUN=".isfac", RETURN=TRUE, BATCHSIZE=100000 )
    
    coarse_df$cell2 <- .create_ff_indices(nrow(coarse_df))
  }
  
  # Convert fine data to data.frame
  #fine_df <- .as_data_frame_factors(fine, xy = TRUE)
  fine_df <- .raster_to_ff_data_frame_xy(fine)
  
  # Add coarse cell ID to fine data.frame
  #fine_df[['cell']] <- as.integer(.create_lut_fine(ids_coarse, fine))
  fine_df[['cell']] <- .create_lut_fine_ff(ids_coarse, fine)
  
  ids_coarse2 <- raster(coarse)
  ids_coarse2[] <- 1:ncell(coarse)
  
  
  print("here12")
  fine_df[['cell2']] <- .create_lut_fine_ff(ids_coarse2, fine)
  if ( add_pycno > 0 || ( input_polygons && data_type == "count") ) { 
    fine_df[['pycnolayer']] <- .create_lut_fine_ff(pycnolayer, fine)
  }
  
  #debug
  #fine_df[ffwhich(fine_df, cell==10157),]
  
  #fine_df <- na.exclude(fine_df)
  fine_df <- .na_exclude_ffdf(fine_df)
  
  # Resampled model onto fine grid
  #fine_df <- cbind(
  #   fine_df[, c('x', 'y', 'cell', 'cell2', nm_covariates)],
  #   .join_interpol(coarse_df = coarse_df[, c('cell', 'cell2', nm_coarse)], fine_df = fine_df, attr = nm_coarse, by = 'cell2')
  # )
  
  fine_df <- .cbind_ffdf2(
    fine_df[c('x', 'y', 'cell', 'cell2', nm_covariates)],
    .join_interpol_ffdf(coarse_df = coarse_df[c('cell', 'cell2', nm_coarse)], fine_df = fine_df, attr = nm_coarse, by = 'cell2')
  )
  
  #for debug
  #plot(.ffdf_to_raster(fine_df, fine_df$layer, res(fine), crs(fine)))
  
  coarse_df <- .na_exclude_ffdf(coarse_df)
  fine_df <- .na_exclude_ffdf(fine_df)
  
  if (is.null(p)) { p = as.numeric( nrow( coarse_df ) / nrow(fine_df) ) }
  # Sub-sample for modelling
  n_spl <- ceiling(nrow(fine_df) * p)
  if ( !is.null(nmax) && nmax > 0 ) {  n_spl <- min(n_spl, nmax) }
  print("here14")
  if ( is.null(sample_method) || sample_method == 'random' ){
    #id_spl <- sample(1:nrow(fine_df), size = n_spl) # sample random grid cells
    ff_inds <- .create_ff_indices(nrow(fine_df))
    id_spl <- .bigsample_in_chunks(ff_inds, n_spl)
  } else {
    print("here15")
    #TODO: Replace line for secondary memory usage
    id_spl <- SpatialPixelsDataFrame(fine_df[, c('y', 'x')], data.frame(fine_df,cell3=1:nrow(fine_df)), proj4string = CRS(projection(fine)))
    print("here16")
    id_spl <- over( spsample( x = id_spl , type=sample_method , n = n_spl ) , id_spl , fn = median )$cell3 # sample grid cells  
  }
  
  if (verbose) message('Selecting best model parameters')
  #y_aux = fine_df[id_spl, nm_coarse, drop = TRUE]  #Nao deveria estar no else do if ( data_type == "count" )?
  tempd_ff <-  fine_df[nm_coarse]
  y_aux <- ffdfindexget(tempd_ff, id_spl)
  
  ########################
  #Data scaling to range 0 - 1
  #range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
  #fine_df$elev <- range01(fine_df$elev)
  #fine_df$cover <- range01(fine_df$cover)
  
  if ( data_type == "count" ) {
    if ( add_pycno > 0 || input_polygons ) {
      y_aux = fine_df[id_spl, 'pycnolayer', drop = TRUE] 
    } else {
      factor = nrow(fine_df) / nrow( coarse_df )
      y_aux = y_aux / as.numeric( factor )
    }
  }
  
  if( method == 'gwr' ) {
    if (data_type != "count" && data_type != "numeric") {
      stop('Data type should be count or numeric, when performing geographically weighted regression')
    }
  } else {
    fit <- .update_model_ffdf( vars = fine_df[id_spl, nm_covariates], y = y_aux, method = method, control = train_control_init, tune_grid = tune_grid, data_type = data_type , long=fine_df$x[id_spl] , lat=fine_df$y[id_spl] )
    best_params <- fit$bestTune
    if (verbose) {
      best_params_str <- paste( lapply(names(best_params), function(x) paste(x, " = ", best_params[[x]], sep = "")), collapse = ", ")
      message("Parameters retained: ", best_params_str)
    }
  }
  print("here17")
  # Initiate matrix to store performance of disseveration
  perf <- matrix(ncol = 3, nrow = 0, dimnames = list(NULL,c("lower_error", "error", "upper_error")))
  # Initiate dissever result data.frame
  #diss_result <- fine_df[, c('x', 'y', 'cell', 'cell2', nm_coarse)]
  diss_result <- fine_df[c('x', 'y', 'cell', 'cell2', nm_coarse)]
  # Our first approximation is actually the nearest neighbour interpolation
  diss_result$diss <- fine_df[[nm_coarse]]
  if ( data_type == "count" ) {
    if ( add_pycno > 0 || input_polygons ) { diss_result$diss <- fine_df[['pycnolayer']] } else {
      factor = nrow(fine_df) / nrow( coarse_df )
      diss_result$diss = diss_result$diss / as.numeric( factor )
    }
  }
  print("here18")
  # Initiate dissever results data.frame aggregated back to coarse grid
  diss_coarse <- coarse_df
  diss_coarse$diss <- coarse_df[[nm_coarse]]
  # Initialising best model selection
  best_fit <- Inf

  for (k in 1:max_iter){
    if (verbose) message('| - iteration ', k)
    if (verbose) message('| -- computing adjustement factor')
    # Calculate adjustment factor
    diss_coarse$adjust <- diss_coarse[[nm_coarse]] / diss_coarse[['diss']]
    # Resample adjustement factor to fine grid
    #diss_result$adjust <- .join_interpol(diss_coarse, fine_df, attr = 'adjust', by = 'cell2')[, 'adjust']
    diss_result$adjust <- .join_interpol_ffdf (diss_coarse, fine_df, attr = 'adjust', by = 'cell2')[['adjust']]
    # Apply adjustement and replace the current
    if ( !( data_type == "categorical" ) ) { diss_result$diss <- diss_result$adjust * diss_result$diss }
    # Sampling new points
    #id_spl <- sample(1:nrow(fine_df), size = n_spl)
    ff_inds <- .create_ff_indices(nrow(fine_df))
    id_spl <- .bigsample_in_chunks(ff_inds, n_spl)
    
    # Update model and update dissever predictions on fine grid
    if( method != 'gwr' ) {
      if (verbose) message('| -- updating model')
      
      #fit <- .update_model( vars = fine_df[id_spl, nm_covariates], y = diss_result[id_spl, 'diss', drop = TRUE], method = method, control = train_control_iter, tune_grid = best_params, data_type = data_type , latLong=data.frame( long=fine_df$x[id_spl] , lat=fine_df$y[id_spl] ))
      fit <- .update_model_ffdf( vars = fine_df[id_spl, nm_covariates], y = diss_result[id_spl, 'diss', drop = TRUE], method = method, control = train_control_iter, tune_grid = best_params, data_type = data_type , long=fine_df$x[id_spl] , lat=fine_df$y[id_spl] )
      
      if (verbose) message('| -- updating predictions')
      if ( method == 'lme' ) { 
        #diss_result$diss <- .predict_map(fit=fit,data.frame(fine_df,dummy=rep.int(1,nrow(fine_df))), split = split_cores, boot = NULL, data_type=data_type, latLong=data.frame( long=fine_df$x , lat=fine_df$y ))
      } else { 
        #diss_result$diss <- .predict_map(fit=fit, fine_df, split = split_cores, boot = NULL, data_type=data_type)
        fineData <- fine_df[c("pycnolayer",nm_covariates)]
        names(fineData) <- c("x", nm_covariates)
        diss_result$diss <- .predict_map_ffdf(fit=fit, fineData, boot = NULL, data_type=data_type)
      }
    } else {
      #TODO: adaptar para uso de memoria secindaria
      varaux = fine_df[id_spl, nm_covariates]
      varr = diss_result[id_spl, 'diss', drop = TRUE]
      datagwr = SpatialPointsDataFrame(fine_df[id_spl, c('x','y')], data.frame(varr, varaux), proj4string = CRS(projection(fine)))
      coordgwr = SpatialPointsDataFrame(fine_df[, c('x','y')], data.frame(fine_df[nm_covariates]), proj4string = CRS(projection(fine)))
      form = as.formula(paste("varr~",paste(names(fine_df[nm_covariates]), collapse="+")))
      if (verbose) message('| -- tuning GWR bandwidth')
      dMat1 <- .gw.dist(dp.locat=as.matrix(fine_df[id_spl, c('x','y')]), rp.locat=as.matrix(fine_df[, c('x','y')]) )
      dMat2 <- .gw.dist(dp.locat=as.matrix(fine_df[id_spl, c('x','y')]), rp.locat=as.matrix(fine_df[id_spl, c('x','y')]) )
      baux <- bw.gwr(form, data = datagwr, kernel="gaussian", longlat=TRUE, adaptive=TRUE, dMat=dMat2 )
      if (verbose) message('| -- updating model')
      fit <- gwr.predict(form, data = datagwr, predictdata = coordgwr, longlat = TRUE, bw = baux, kernel="gaussian", adaptive=TRUE, dMat1=dMat1 , dMat2=dMat2 )
      if (verbose) message('| -- updating predictions')
      diss_result$diss = fit$SDF$prediction
    }
    
    diss_result <- .mass_preserv_ffdf(diss_result)
    
    if (verbose) message('| -- computing aggregates of predictions on coarse grid')
    summary_metric <- function( data, type ) {
      #print(data)
      if (type == 'count') { return( sum(data) ) }
      if (type == 'categorical') { return( median(data) ) }
      return( mean(data) )
    }
    
    #diss_coarse <- diss_result %>% group_by(cell) %>% summarise(diss = summary_metric(diss,data_type)) %>% inner_join(coarse_df, ., by = "cell")
    #debug
    #diss_coarse_dummy <- diss_result_dummy %>% group_by(cell) %>% summarise(diss = summary_metric(diss,data_type)) %>% inner_join(coarse_df_dummy, ., by = "cell")
    ##############################################
    
    splitby <- as.character.ff(diss_result$cell)
    
    grp_qty <- ffdfdply(x=diss_result[c("cell","diss", nm_coarse)], 
                        split=splitby, 
                        FUN = function(x){
                          ## This happens in RAM - containing **several** split elements so here we can use data.table which works fine for in RAM computing
                          require(data.table)
                          data <- as.data.table(x)
                          result <- data[, list(diss = summary_metric(diss,data_type), count = length(diss), layer=first(layer)), by = list(cell)]
                          as.data.frame(result)
                        })
    
    diss_coarse <- merge(coarse_df,grp_qty[c("cell","diss")], by.x = "cell", by.y = "cell", all.x=FALSE, all.y=FALSE, trace = TRUE)
    
    #Debug
    #grp_qty$corrections[ffwhich(grp_qty, cell==cell)]
    # for(i in 1:nrow(diss_coarse_dummy)){
    #   a <- diss_coarse_dummy[i,6]
    #   b <- diss_coarse$diss[i]
    #   if (a != b){print("a=",a," b=",b)}
    # }
    
    if (verbose) message('| -- computing performance stats')
    # compute error
    n <- nrow(diss_coarse)
    sqe <- (diss_coarse[[nm_coarse]] - diss_coarse[['diss']])^2
    mse <- mean(sqe)
    error <- sqrt(mse) # RMSE
    # Confidence intervals (in this case we use 5% C.I.)
    t_student <- qt(1 - (0.05/2), df = n - 1) # 0.975 quantile from Student t distribution
    var <- ((1)^2/(n * (n - 1))) * sum(sqe) # Variance
    se <- sqrt(var) # Standard error
    ci <- se * t_student
    upper_ci <- sqrt(mse + ci)
    lower_ci <- sqrt(max(0, mse - ci))
    if (data_type == 'categorical' ) {
      aux <- (diss_coarse[[nm_coarse]] == diss_coarse[['diss']])
      error <- 1.0 - ( sum(aux) / as.numeric(n) )
      aux <- replicate( 1000 , 1.0 - ( sum(sample(aux, 0.5 * n )) / ( 0.5 * n ) ) )
      upper_ci = quantile(aux, probs=0.975)
      lower_ci = quantile(aux, probs=0.025)
    }
    # Fill result matrix
    perf <- rbind(perf, c(lower_ci, error, upper_ci))
    if (verbose) message("| -- ERROR = ", round(error, 3))
    # Choose whether we retain the model or not
    if (error < best_fit) {
      best_fit <- error
      best_iteration <- k
      best_model <- fit
    }
    # We only test improvement if more than 3 iterations
    if (k >= min_iter & k >= 3) {
      # Computing stop criterion
      stop_criterion <- mean(
        (perf[k - 2, 2] - perf[k - 1, 2]) + # improvement at last iteration
          (perf[k - 1, 2] - error) + # improvement at this iteration
          (perf[k - 2, 2] - error) # improvement over last 2 iterations
      )
      # If we have reach some kind of pseudo-optimium, we finish the iteration stage
      if (stop_criterion <= thresh) break
    }
  }
  
  if (verbose) message('Retaining model fitted at iteration ', best_iteration)
  if( method == 'gwr' ) { map <- fit$SDF$prediction } else {
    if ( method == 'lme' ) map <- .predict_map(fit=best_model,data.frame(fine_df,dummy=rep.int(1,nrow(fine_df))), boot = boot, level = level, data_type=data_type, latLong=data.frame( long=fine_df$x , lat=fine_df$y ))
    else map <- .predict_map_ffdf(fit=best_model, fine_df, boot = boot, level = level, data_type=data_type)
  }
  
  diss_result$diss <- map
  diss_result <- .mass_preserv_ffdf(diss_result)
  
  #Data scaling to range 0 - 1
  range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
  diss_result$diss <- range01(diss_result$diss)
  diss_result$diss <- diss_result$diss + 1
  
  # if (data_type == 'count') { 
  #   #map[map < 0.0] <- 0 
  #   temp_inds <- ffwhich(map, map< 0.0)
  #   if(!is.null(temp_inds)){map[temp_inds] <- ff(0, length = length(temp_inds))}
  # }
  
  #map <- rasterFromXYZ( data.frame( diss_result[, c('x', 'y')], diss = map ), res = res(fine), crs = projection(fine) )
  map <- .ffdf_to_raster(diss_result=diss_result, map=diss_result$diss, res=res(fine) ,crs=projection(fine))
  if( method == 'gwr' ) {
    res <- list( fit = fit$SDF, map = map, perf = data.frame(perf) )
  } else {
    res <- list( fit = fit, map = map, perf = data.frame(perf) )
  }
  class(res) <- c(class(res), 'dissever')
  

	stopCluster(cl)
  
  # end.time <- Sys.time()
  # timeTotal <- difftime(end.time, start.time, units = "secs")
  # assign("timeTotal",as.numeric(timeTotal),.GlobalEnv)
  
  return(res)
}

#' @name plot.dissever
#' @title Plots a dissever result
#' @description Plots a dissever result. Two modes are possible to visualise either the resulting map or the convergence of the disseveration.
#' @param x object of class \code{dissever}, output from the \code{dissever} function
#' @param type character, type of visualisation to produce. Options are "map", to produce a map of the dissevered coarse map, or "perf", to show the convergence of the error during the disseveration procedure.
#' @param ... Additional arguments passed to plot
#' @author Pierre Roudier
#' @examples
#' # See ?dissever
plot.dissever <- function(x, type = 'map', ...) {
  if (! type %in% c('map', 'perf')) stop('Invalid type of plot.')
  if (type == 'map') {
    plot(x$map, col = viridis(100), ...)
  } else {
    n_iter <- nrow(x$perf)
    plot(1:n_iter, x$perf$error, type = 'l', xlab = 'Iteration', ylab = 'Error', ylim = range(x$perf), ...)
    lines(1:n_iter, x$perf$upper_error, lty = 3)
    lines(1:n_iter, x$perf$lower_error, lty = 3)
    best <- which.min(x$perf$error)
    points(best, x$perf[best, 'error'], pch = 16, col = 2)
  }
}

#' @name print.dissever
#' @title Prints the performance of the dissever procedure
#' @description Prints the performance of the model used to do the dissever procedure.
#' @param x object of class \code{dissever}, output from the \code{dissever} function
#' @param ... Additional arguments passed to print
#' @author Pierre Roudier
print.dissever <- function(x, ...) { print(x$fit, ...) }

#' @name summary.dissever
#' @title Prints summary of the model used in the dissever procedure
#' @description Prints summary of the model used in the dissever procedure.
#' @param object object of class \code{dissever}, output from the \code{dissever} function
#' @param ... Additional arguments passed to summary
#' @author Pierre Roudier
summary.dissever <- function(object, ...) { summary(object$fit, ...) }

if(!isGeneric("generate_ci")) setGeneric("generate_ci", function(object, covariates, ...) { standardGeneric("generate_ci") })

#' @name generate_ci
#' @aliases generate_ci,list,RasterStack-method
#' @title Confidence intervals using bootstraping
#' @description Generates confidence intervals of a dissever output using bootstraping
#' @param object object of class \code{dissever}, output from the \code{dissever} function
#' @param covariates object of class \code{"RasterStack"}, the fine-resolution stack of predictive covariates used to generate the dissever output
#' @param level If this is a numeric value, it is used to derive confidence intervals using quantiles. If it is a function, it is used to derive the uncertainty using this function.
#' @param n the number of bootstrap replicates used to derive the confidence intervals
#' @docType methods
#' @author Pierre Roudier
#' @examples
#' # Load the Edgeroi dataset (see ?edgeroi)
#' data(edgeroi)
#'
#' # Create a dissever output
#' diss <- dissever(
#'   coarse = edgeroi$carbon,
#'   fine = edgeroi$predictors,
#'   method = "lm",
#'   min_iter = 5, max_iter = 10,
#'   p = 0.05
#' )
#'
#' # Generate the confidence intervals
#' \dontrun{
#' ci <- generate_ci(diss, edgeroi$predictors, n = 5)
#'
#' plot(ci)
#' }
setMethod('generate_ci', signature(object = "list", covariates = "RasterStack"), .generate_ci )

if(!isGeneric("dissever")) setGeneric("dissever", function(coarse, fine, ...) { standardGeneric("dissever") })

#' @title Spatial downscaling
#' @name dissever
#' @aliases dissever,RasterStack-method
#' @description Performs spatial downscaling of coarse grid mapping to fine grid mapping using predictive covariates and a model fitted using the caret package.
#' @param coarse, object of class \code{"RasterLayer"}, the coarse-resolution layer that needs to be downscaled, or alternatively an object of class \code{"SpatialPolygonsDataFrame"} with two attributes
#' @param fine, object of class \code{"RasterStack"}, the fine-resolution stack of predictive covariates
#' @param coarse_var_names, names for two attributes from the \code{"SpatialPolygonsDataFrame"} passed as the coarse object, corresponding to IDs for the different regions and to the per-region values that are to be downscaled
#' @param method, a string specifying which classification or regression model to use (via the caret package). Possible values are found using names(caret::getModelInfo()).
#' @param p, numeric proportion of the fine map that is sampled for fitting the dissever model (between 0 and 1, defaults to the ratio between the coarse grid resolution and the fine grid resolution)
#' @param sample_method, string specifying the method used for sampling the fine map through the spsample function  
#' @param nmax numeric maximum number of pixels selected for fitting the dissever model. It will override the number of pixels chosen by the \code{p} option if that number is over the value passed to \code{nmax}.
#' @param thresh numeric, dissever iterations will proceed until the error of the dissever model reaches this value, or until the maximum number of iterations is met (defaults to 0.01)
#' @param min_iter numeric, minimum number of iterations (defaults to 5)
#' @param max_iter numeric, maximum number of iterations (defaults to 100)
#' @param boot numeric, if not NULL (default), the number of bootstrap replicates used to derive the confidence intervals.
#' @param level If this is a numeric value, it is used to derive confidence intervals using quantiles. If it is a function, it is used to derive the uncertainty using this function.
#' @param tune_length numeric, the number of parameters to test to find the optimal parametrisation of the caret model (defaults to 3)
#' @param tune_grid a data frame with possible tuning values
#' @param train_control_init Control parameters for finding the optimal parameters of the caret model (see trainControl)
#' @param train_control_iter Control parameters for fitting the caret model during the iteration phase (see trainControl)
#' @param data_type a string indicating the type of data to be downscaled/disaggregated. Can be 'numeric', 'count' or 'categorical' (defaults to 'numeric')
#' @param add_pycno controls if the results of pycnophylactic interpolation should be used as initialization (> 0)
#' @param verbose controls the verbosity of the output (TRUE or FALSE)
#' @docType methods
#' @author Brendan Malone, Pierre Roudier, Bruno Martins, João Cordeiro
#' @references Malone, B.P, McBratney, A.B., Minasny, B., Wheeler, I., (2011) A general method for downscaling earth resource information. Computers & Geosciences, 41: 119-125. \url{http://dx.doi.org/10.1016/j.cageo.2011.08.021}
#' @examples
#' # Load the Edgeroi dataset (see ?edgeroi)
#' data(edgeroi)
#'
#' # Plot the Edgeroi dataset (using the raster package)
#' library(raster)
#' plot(edgeroi$carbon) # coarse resolution layer
#' plot(edgeroi$predictors) # fine resolution predictors
#'
#' # Run dissever using a simple linear model.
#'
#' # In this instance we are subsampling heavily (p = 0.05) to keep
#' # run time short
#' res_lm <- dissever(
#'   coarse = edgeroi$carbon,
#'   fine = edgeroi$predictors,
#'   method = "lm",
#'   min_iter = 5, max_iter = 10,
#'   p = 0.05
#' )
#'
#' # A lot more models are available through caret:
#' \dontrun{
#' subset(caret::modelLookup(), forReg == TRUE, select = 'model')
#' }
#'
#' # Plot dissever results
#' plot(res_lm, type = 'map', main = "Dissever using GAM")
#' plot(res_lm, type = 'perf', main = "Dissever using GAM")
#'
setMethod( 'dissever', signature(fine = "RasterStack"), .dissever )

if(!isGeneric("pycno")) setGeneric("pycno", function(x, pops, celldim, ...) { standardGeneric("pycno") })

setMethod( 'pycno', signature(x = "SpatialPolygonsDataFrame"), .pycno )
