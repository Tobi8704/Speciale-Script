# Functions for clustered standard errors as provided by Ian Gow: 
# http://www.people.hbs.edu/igow/GOT/Code/cluster2.R.html, 2017-12-22.
coeftest.cluster_grøn <- function(data, fm, cluster1 = NULL, cluster2 = NULL, ret = 'test') {
  options(warn = -1)
  
  # Return White (1980) standard errors if no cluster
  # variable is provided
  if (is.null(cluster1)) {
    if (ret == 'cov') {
      return(sandwich::vcovHC(fm, type = 'HC0'))
    } else {
      return(coeftest(fm, vcov = vcovHC(fm, type = 'HC0')))
    }
  }
  
  # Calculation shared by covariance estimates
  est.fun <- sandwich::estfun(fm)
  
  # Need to identify observations used in the regression (i.e.,
  # non-missing) values, as the cluster vectors come from the full
  # data set and may not be in the regression model.
  inc.obs <- !is.na(est.fun[, 1])
  est.fun <- est.fun[inc.obs, ]
  
  # Shared data for degrees-of-freedom corrections
  N <- dim(fm$model)[1]
  NROW <- NROW(est.fun)
  K <- fm$rank
  
  # Calculate the sandwich covariance estimate
  cov <- function(cluster) {
    cluster <- factor(cluster, exclude = NULL)
    
    # Calculate the "meat" of the sandwich estimators
    u <- apply(est.fun, 2, function(x) tapply(x, cluster, sum))
    meat <- crossprod(u) / N
    
    # Calculations for degrees-of-freedom corrections, followed
    # by calculation of the variance-covariance estimate.
    # NOTE: NROW/N is a kluge to address the fact that sandwich
    # uses the wrong number of rows (includes rows omitted from
    # the regression).
    M <- length(levels(cluster))
    dfc <- M / (M-1) * (N-1) / (N-K)
    
    return(dfc * NROW / N * sandwich::sandwich(fm, meat = meat))
  }
  
  # Modified approach to handle cluster variables
  if (is.character(cluster1)) {
    # It's a column name - extract the column from the data frame
    cluster1_vec <- data[[cluster1]]
  } else {
    # It's already a vector
    cluster1_vec <- cluster1
  }
  
  # Make sure we subset to match the observations in the model
  cluster1_vec <- cluster1_vec[inc.obs]
  cov1 <- cov(cluster1_vec)
  
  if (is.null(cluster2)) {
    # If only one cluster supplied, return single cluster results
    if (ret == 'cov') {
      return(cov1)
    } else {
      return(lmtest::coeftest(fm, cov1))
    }
  } else {
    # Handle the second cluster variable the same way
    if (is.character(cluster2)) {
      cluster2_vec <- data[[cluster2]]
    } else {
      cluster2_vec <- cluster2
    }
    
    cluster2_vec <- cluster2_vec[inc.obs]
    cluster12 <- paste(cluster1_vec, cluster2_vec, sep = '')
    
    # Calculate the covariance matrices for cluster2, the "intersection"
    # cluster, then put all the pieces together.
    cov2 <- cov(cluster2_vec)
    cov12 <- cov(cluster12)
    covMCL <- (cov1 + cov2 - cov12)
    
    # Return the output of coeftest using two-way cluster-robust
    # standard errors.
    if (ret == 'cov') {
      return(covMCL)
    } else {
      return(lmtest::coeftest(fm, covMCL))
    }
  }
  options(warn = 0)
}

summary.cluster_grøn <- function( obj , data , cluster1 , cluster2 = NULL , alpha = 0.05 ) {
  # Following based on suggestion from
  # https://stat.ethz.ch/pipermail/r-help/2011-January/264777.html
  # provided by Achim Zeileis.
  options( warn = -1 )
  # Get original summary
  s <- memisc::getSummary( obj , alpha = alpha )
  
  ## replace Wald tests of coefficients
  s$coef[ , 1 : 4 , 1 ] <- coeftest.cluster_grøn( data , obj , cluster1 , cluster2 )
  
  ## replace confidence intervals
  crit <- qt( alpha / 2 , obj$df.residual )
  s$coef[ , 5 , 1 ] <- s$coef[ , 1 , 1 ] + crit * s$coef[ , 2 , 1 ]
  s$coef[ , 6 , 1 ] <- s$coef[ , 1 , 1 ] - crit * s$coef[ , 2 , 1 ]
  
  # Note that some components of s$sumsstat will be inconsistent with
  # the clustered calculations
  
  return( s )
  options( warn = 0 )
}



## RD Functions ####

jump.plot_grøn <- function( data , force.var , yvar , seat.identifier , polynomial ){
  data <- data[ , c( force.var , yvar  , seat.identifier )]
  data <- na.omit( data )
  library( ggplot2 )
  p <- ggplot( ) +
    geom_point( data = data 
                , aes_string( x = force.var, y = yvar , shape = seat.identifier ) 
                , size = 2 ) +
    geom_smooth( data = subset( data , data[ , force.var] < 0 )
                 , aes_string ( x = force.var , y = yvar )
                 , method = 'lm' , formula = y ~ poly( x , polynomial , raw = TRUE )
                 , linetype = 1 , color = 'black' , size = 1 ) +
    geom_smooth( data = subset( data , data[ , force.var] >= 0 )
                 , aes_string ( x = force.var , y = yvar )
                 , method = 'lm' , formula = y ~ poly( x , polynomial , raw = TRUE )
                 , linetype = 1 , color = 'black' , size = 1 ) +
    scale_x_continuous( name = 'Grønne partiers stemmeandel'
                        , limits = c( -5 , 10 )
                        , breaks = seq( -5 , 10 , 2.5 )) +
    scale_y_continuous( name = 'Miljøfokus' 
                        , limits = c( -8 , 8 )
                        , breaks = seq( -8 , 8 , 4 )) +
    scale_shape_manual( values = c( 1 , 19 ) 
                        , labels = c( 'Grønne partier w/o seats     ' , 'Grønne partier w seat(s)     ')) +
    geom_vline( xintercept = 0 , linetype = 2 , size = .6 ) +
    theme( legend.position = 'bottom' , legend.title = element_blank())
  return( p )
  detach( package:ggplot2 )
}


rd.core_grøn <- function( data , force.var , yvar , seat.identifier , fixed.effects 
                     , clust1 , clust2 , polynomial , bws ){
  i <- polynomial
  data <- as.data.frame( data )
  data <- data[ , c( yvar , force.var , fixed.effects , seat.identifier , clust1 , clust2 )]
  
  if( i <= 2 & is.null( bws )){
    h <- rdd::IKbandwidth( X = data[ , force.var ] , Y = data[ , yvar ]
                           , cutpoint = 0 , kernel = 'triangular' )
    data$w <- rdd::kernelwts( X = data[,force.var] , center = 0 
                              , bw = h,  kernel = 'triangular' )
  }
  if( i <= 2 & !is.null( bws )){
    h <- bws
    data$w <- rdd::kernelwts( X = data[ , force.var ] , center = 0
                              , bw = h , kernel = "triangular" )
  }
  
  if( !is.null( clust1 )){ data[ , clust1 ] <- as.factor( data[ , clust1 ])}
  if( !is.null( clust2 )){ data[ , clust2 ] <- as.factor( data[ , clust2 ])}
  data <- na.omit( data )
  
  data$above[ data[ , force.var] >= 0 & !is.na( data[ , force.var ])] <- 1
  data$above[ data[ , force.var] < 0 & !is.na( data[ , force.var ])] <- 0
  data$force_above <- data[ , force.var] * data$above
  
  formula = as.formula( paste( yvar, "~" , seat.identifier , " + poly (" , force.var , "," , i , " , raw = TRUE ) +
                              poly( force_above , " , i , " , raw = TRUE ) + as.factor( ", fixed.effects , " ) |
                              above + poly( " , force.var , "," , i , " , raw = TRUE )+
                              poly( force_above , " , i , " , raw = TRUE ) + as.factor( " , fixed.effects , " )" ))
  if( i <= 2 ){                  
    ivreg <- AER::ivreg( formula = formula
                         , weights = w
                         , data = subset( data , w > 0 ))
    data2 <- subset( data , w > 0 ) 
  }
  if( i > 2 ){ 
    ivreg <- AER::ivreg( formula = formula , data = data )
    data2 <- data
  }
  ivreg.out <- summary( ivreg )
  
  data2[ , clust1 ] <- as.factor( as.character( data2[ , clust1 ]))
  data2[ , clust2 ] <- as.factor( as.character( data2[ , clust2 ]))
  
  coeftest.cluster_grøn( data2 , ivreg , cluster1 = clust1 , cluster2 = clust2 )
  coef <- summary.cluster_grøn( ivreg , data2 , cluster1 = clust1 , cluster2 = clust2 , alpha = 0.05 )
  return( coef )
  
  
}


rd.base_grøn <- function( data , force.var , yvar , seat.identifier , fixed.effects 
                     , clust1 , clust2 , polynomials , bws ){
  data <- as.data.frame( data )
  data <- data[ , c( force.var , yvar , seat.identifier , fixed.effects , clust1 , clust2 )]
  data <- na.omit( data )
  
  for ( i in polynomials ){
    coef <- rd.core_grøn( data = data 
                     , force.var = force.var 
                     , yvar = yvar
                     , seat.identifier = seat.identifier 
                     , fixed.effects = fixed.effects 
                     , clust1 = clust1 
                     , clust2 = clust2 
                     , polynomial = i 
                     , bws = NULL )
    coef <- as.data.frame( t ( coef$coef[ 2 ,  , 1 ] ))
    if( i > 2 ){
      coef$IK_BW <- 'global'
      coef$Estimation <- 'Parametric'
      Nleft <- as.character( nrow( subset( data , data[ , force.var ] < 0 )))
      Nright <- as.character( nrow( subset( data , data[ , force.var ] >= 0 )))
    }
    if( i <= 2 & !is.null( bws )){
      coef$IK_BW <- bws
      coef$IK_BW <- sprintf( '%.3f' , round( coef$IK_BW , 3 ))
      coef$Estimation <- 'Non-Parametric'
      coef$Nleft <- as.character( nrow( subset( data.cut , data.cut[ , force.var ] >= h * -1 &  data.cut[ , force.var ] < 0 )))
      coef$Nright <- as.character( nrow( subset( data.cut , data.cut[ , force.var ] <= h &  data.cut[ , force.var ] >= 0 )))
    }
    if( i <= 2  & is.null( bws )){
      h <- rdd::IKbandwidth( X = data[ , force.var ] , Y = data[ , yvar ]
                             , cutpoint = 0 , kernel = 'triangular' )
      coef$IK_BW <- h
      coef$IK_BW <- sprintf( '%.3f' , round( coef$IK_BW , 3 ))
      coef$Estimation <- 'Non-Parametric'
      Nleft <- as.character( nrow( subset( data , data[ , force.var ] >= h * -1 &  data[ , force.var ] < 0 )))
      Nright <- as.character( nrow( subset( data , data[ , force.var ] <= h &  data[ , force.var ] >= 0 )))
    }
    
    coef$stars[ coef[ , 4 ] > .1 ] <- ''
    coef$stars[ coef[ , 4 ] <= .1 ] <- '*'
    coef$stars[ coef[ , 4 ] <= .05 ] <- '**'
    coef$stars[ coef[ , 4 ] <= .01 ] <- '***'
    coef$est <- sprintf( '%.3f' , round( coef$est , 3 ))
    coef$est <- paste( coef$est , coef$stars , sep = '' )
    coef$Poly <- i
    coef$Poly <- as.character( coef$Poly )
    coef$Nleft <- Nleft
    coef$Nright <- Nright
    coef[ , 'stars' ] <- NULL
    
    if( exists( 'return.ds' )){
      return.ds <- rbind( return.ds , coef )
    }
    if ( !exists( 'return.ds' )){
      return.ds <- coef
    }
  }
  return.ds <- return.ds[ , c( 1:2 , 4 , 7:11 ) ]
  
  colnames( return.ds ) <- c( 'LATE' , 'St. Err.' , 'p-value', 'Bandwith'
                              , 'Approach' , 'Polynomial' , 'N left of c' 
                              , 'N right of c' ) 
  return( return.ds )
}