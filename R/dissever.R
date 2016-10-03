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
  if ( method == 'gwrm' ) { 
    fit <- gw( as.formula(paste("x~",paste(names(vars), collapse="+"))) , data= data.frame( vars , x=y_aux ) )
  } else if ( method == 'lme' ) {
    fit <- lme( fixed=as.formula("x ~ . - dummy - lat - long") , data=data.frame( vars , latLong, x=y_aux , dummy=rep.int( 1 , length(y_aux) ) ) , random = ~ 1 | dummy, method = "ML" , correlation = corGaus(form = ~ lat+long | dummy ) )
    print(fit)
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
    else if ( method == 'lme' ) bootstrap_fit <- lme( fixed=as.formula(".outcome ~ . - dummy - lat - long") , data=data.frame( bootstrap_df , latLong[idx, ], dummy=rep.int( 1 , nrow(bootstrap_df) ) ) , random = ~ 1 | dummy, method = "ML" , correlation = corGaus(form = ~ lat+long | dummy ) )
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
        res <- predict(object=fit , newdata = data[split == i, ])
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

.dissever <- function(
    coarse,
    fine,
    coarse_var_names = NULL,
    method = "rf",
    p = NULL,
    sample_method = 'regular',
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
    data_type = "numeric",
    add_pycno = 0,
    verbose = FALSE
  ) {
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
    if ( add_pycno > 0 ) { pycnolayer <- raster( pycno( coarse, coarse[[coarse_var_names[2]]], min(minres), converge=add_pycno, verbose=FALSE ) ) }
    else if ( data_type == "count" ) { pycnolayer <- raster( pycno( coarse, coarse[[coarse_var_names[2]]], min(minres), converge=0, verbose=FALSE ) ) }    
    ids_coarse <- rasterize(coarse, raster( resolution=minres * 1.01, ext=extent(coarse) ), coarse_var_names[1], fun='first')
    names(ids_coarse) <- 'cell'
    coarse <- rasterize(coarse, raster( resolution=minres * 1.01, ext=extent(coarse) ), coarse_var_names[2], fun='first')    
  } else if ( add_pycno > 0 ) {
    minres <- min(res(fine))
    pycnolayer <- raster( pycno( rasterToPolygons(coarse), .as_data_frame_factors(coarse), 0.05, converge=add_pycno, verbose=FALSE ) )
  }
  if (min(res(fine)) >= min(res(coarse))) { stop('Resolution of fine data should be higher than resolution of coarse data') }
  if (!(data_type == "numeric" || data_type == "count" || data_type == "categorical" )) {
    stop('Data type should be numeric, categorical or count')
  }
  # Store names of coarse data and fine-scale covariates
  nm_coarse <- names(coarse)
  nm_covariates <- names(fine)
  # Get cell numbers of the coarse grid and convert coarse data to data.frame
  if ( !input_polygons ) {
    ids_coarse <- raster(coarse)
    ids_coarse[] <- 1:ncell(coarse)
    names(ids_coarse) <- 'cell'
    coarse_df <- .as_data_frame_factors(coarse, xy = TRUE)
    coarse_df$cell <- 1:nrow(coarse_df)
    coarse_df$cell2 <- coarse_df$cell
  } else {
    coarse_df <- .as_data_frame_factors(coarse, xy = TRUE)
    coarse_df$cell <- .as_data_frame_factors(ids_coarse, xy = TRUE)[['cell']]
    coarse_df$cell <- sapply(coarse_df$cell, function(x) if(is.factor(x)) { as.numeric(x) } else { x })
    coarse_df$cell2 <- 1:nrow(coarse_df)
  } 
  # Convert fine data to data.frame
  fine_df <- .as_data_frame_factors(fine, xy = TRUE)
  # Add coarse cell ID to fine data.frame
  fine_df[['cell']] <- as.integer(.create_lut_fine(ids_coarse, fine))
  ids_coarse2 <- raster(coarse)
  ids_coarse2[] <- 1:ncell(coarse)
  fine_df[['cell2']] <- as.integer(.create_lut_fine(ids_coarse2, fine))
  if ( add_pycno > 0 || ( input_polygons && data_type == "count") ) { fine_df[['pycnolayer']] <- as.integer(.create_lut_fine(pycnolayer, fine)) }
  fine_df <- na.exclude(fine_df)
  # Resampled model onto fine grid
  fine_df <- cbind(
    fine_df[, c('x', 'y', 'cell', 'cell2', nm_covariates)],
    .join_interpol(coarse_df = coarse_df[, c('cell', 'cell2', nm_coarse)], fine_df = fine_df, attr = nm_coarse, by = 'cell2')
  )
  coarse_df <- na.exclude(coarse_df)
  fine_df <- na.exclude(fine_df)
  if (is.null(p)) { p = as.numeric( nrow( coarse_df ) / nrow(fine_df) ) }
  # Sub-sample for modelling
  n_spl <- ceiling(nrow(fine_df) * p)
  if ( !is.null(nmax) && nmax > 0 ) {  n_spl <- min(n_spl, nmax) }                             
  if ( is.null(sample_method) || sample_method == 'random' ) id_spl <- sample(1:nrow(fine_df), size = n_spl) # sample random grid cells
  else {
    id_spl <- SpatialPixelsDataFrame(fine_df[, c('y', 'x')], data.frame(fine_df,cell3=1:nrow(fine_df)), proj4string = CRS(projection(fine)))
    id_spl <- over( spsample( x = id_spl , type=sample_method , n = n_spl ) , id_spl , fn = median )$cell3 # sample grid cells  
  }
  if (verbose) message('Selecting best model parameters')
  y_aux = fine_df[id_spl, nm_coarse, drop = TRUE]  
  if ( data_type == "count" ) { 
     if ( add_pycno > 0 || input_polygons ) { y_aux = fine_df[id_spl, 'pycnolayer', drop = TRUE] } else {
      factor = nrow(fine_df) / nrow( coarse_df )
      y_aux = y_aux / as.numeric( factor )
     }
  }
  if( method == 'gwr' ) {
    if (data_type != "count" && data_type != "numeric") {
      stop('Data type should be count or numeric, when performing geographically weighted regression')
    }
  } else {
    fit <- .update_model( vars = fine_df[id_spl, nm_covariates], y = y_aux, method = method, control = train_control_init, tune_grid = tune_grid, data_type = data_type , latLong=data.frame( long=fine_df$x[id_spl] , lat=fine_df$y[id_spl] ) )
    best_params <- fit$bestTune
    if (verbose) {
      best_params_str <- paste( lapply(names(best_params), function(x) paste(x, " = ", best_params[[x]], sep = "")), collapse = ", ")
      message("Parameters retained: ", best_params_str)
    }
  }
  # Initiate matrix to store performance of disseveration
  perf <- matrix(ncol = 3, nrow = 0, dimnames = list(NULL,c("lower_error", "error", "upper_error")))
  # Initiate dissever result data.frame
  diss_result <- fine_df[, c('x', 'y', 'cell', 'cell2', nm_coarse)]
  # Our first approximation is actually the nearest neighbour interpolation
  diss_result$diss <- fine_df[[nm_coarse]]
  if ( data_type == "count" ) {
    if ( add_pycno > 0 || input_polygons ) { diss_result$diss <- fine_df[['pycnolayer']] } else {
     factor = nrow(fine_df) / nrow( coarse_df )
     diss_result$diss = diss_result$diss / as.numeric( factor )
    }
  }
  # Initiate dissever results data.frame aggregated back to coarse grid
  diss_coarse <- coarse_df
  diss_coarse$diss <- coarse_df[[nm_coarse]]
  # Initialising best model selection
  best_fit <- Inf
  # If parallel computing: Get split vector, assigning each row in fine data to a core (for prediction of models on fine grid)
  if (.has_parallel_backend()) {
    n_cores <- getDoParWorkers()
    split_cores <- .get_split_idx(nrow(fine_df), p = n_cores)
  } else {
    split_cores <- NULL
  }
  for (k in 1:max_iter){
    if (verbose) message('| - iteration ', k)
    if (verbose) message('| -- computing adjustement factor')
    # Calculate adjustment factor
    diss_coarse$adjust <- diss_coarse[[nm_coarse]] / diss_coarse[['diss']]
    # Resample adjustement factor to fine grid
    diss_result$adjust <- .join_interpol(diss_coarse, fine_df, attr = 'adjust', by = 'cell2')[, 'adjust']
    # Apply adjustement and replace the current
    if ( !( data_type == "categorical" ) ) { diss_result$diss <- diss_result$adjust * diss_result$diss }
    # Sampling new points
    id_spl <- sample(1:nrow(fine_df), size = n_spl)
    # Update model and update dissever predictions on fine grid
    if( method != 'gwr' ) {
      if (verbose) message('| -- updating model')
      fit <- .update_model( vars = fine_df[id_spl, nm_covariates], y = diss_result[id_spl, 'diss', drop = TRUE], method = method, control = train_control_iter, tune_grid = best_params, data_type = data_type , latLong=data.frame( long=fine_df$x[id_spl] , lat=fine_df$y[id_spl] ))
      if (verbose) message('| -- updating predictions')
      if ( method == 'lme' ) diss_result$diss <- .predict_map(fit=fit,data.frame(fine_df,dummy=rep.int(1,nrow(fine_df))), split = split_cores, boot = NULL, data_type=data_type, latLong=data.frame( long=fine_df$x , lat=fine_df$y ))
      else diss_result$diss <- .predict_map(fit=fit, fine_df, split = split_cores, boot = NULL, data_type=data_type)
    } else {
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
    if (data_type == 'count') { diss_result$diss[diss_result$diss < 0.0] <- 0 }
    if (verbose) message('| -- computing aggregates of predictions on coarse grid')
    summary_metric <- function( data, type ) {
      if (type == 'count') { return( sum(data) ) }
      if (type == 'categorical') { return( median(data) ) }
      return( mean(data) )
    }
    diss_coarse <- diss_result %>% group_by(cell) %>% summarise(diss = summary_metric(diss,data_type)) %>% inner_join(coarse_df, ., by = "cell")
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
    if ( method == 'lme' ) map <- .predict_map(fit=best_model,data.frame(fine_df,dummy=rep.int(1,nrow(fine_df))), split = split_cores, boot = boot, level = level, data_type=data_type, latLong=data.frame( long=fine_df$x , lat=fine_df$y ))
    else map <- .predict_map(fit=best_model, fine_df, split = split_cores, boot = boot, level = level, data_type=data_type)
  }
  if (data_type == 'count') { map[map < 0.0] <- 0 }
  map <- rasterFromXYZ( data.frame( diss_result[, c('x', 'y')], diss = map ), res = res(fine), crs = projection(fine) )
  if( method == 'gwr' ) {
    res <- list( fit = fit$SDF, map = map, perf = data.frame(perf) )
  } else {
    res <- list( fit = fit, map = map, perf = data.frame(perf) )
  }
  class(res) <- c(class(res), 'dissever')
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

if(!isGeneric("generate_ci")) {
  setGeneric("generate_ci", function(object, covariates, ...) { standardGeneric("generate_ci") })
}

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

if(!isGeneric("dissever")) {
  setGeneric("dissever", function(coarse, fine, ...) { standardGeneric("dissever") })
}

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
