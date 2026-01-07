#' A k-fold cross-validation procedure for optimal graphical lasso tuning parameter
#'
#' @param data A data matrix containing data of a single component (Tmatrix%*%t(Pmatrix))
#' @param nfolds Number of folds. If missing, 10-fold cross-validation will be performed
#' @param GraphicalLasso The range of graphical lasso tuning parameters. Default is a sequence of 100 numbers from 0.0000001 to 0.05.
#'
#' @returns A graphical lasso tuning parameter that leads to a moddel with PRESS closest to the lowest PRESS + 1SE.
#' @export
#'
cv_graphicalLasso <- function(data, nfolds = 10, GraphicalLasso = exp(seq(from = log(0.0000001), to = log(0.05), length.out = 100))){
  sd_MSE = NULL
  error_x <- NULL
  press = NULL
  id = sample(1:nfolds, nrow(data), replace = TRUE)
  list = 1:nfolds

  for (gl in 1:length(GraphicalLasso)){
    for (k in 1:nfolds){
      train = subset(data, id%in%list[-k])
      test = subset(data, id%in%c(k))

      ridge = matrix(0, nrow = ncol(train), ncol = ncol(train))
      diag(ridge) = 1
      train_ridge = stats::cov(train) + (ridge*1/100000) #add ridge penalty for increased numeric stability
      SigmaINV_hat = glasso::glasso(train_ridge, GraphicalLasso[gl], penalize.diagonal = FALSE)$wi #estimate inverse covariance matrix

      #Make the inverse covariance matrix symmetric
      newINV = matrix(nrow = ncol(SigmaINV_hat), ncol = ncol(SigmaINV_hat))
      for (i in 1:ncol(SigmaINV_hat)){
        for (j in 1:ncol(SigmaINV_hat)){
          newINV[i,j] = newINV[j,i] = mean(c(SigmaINV_hat[i,j], SigmaINV_hat[j,i]))
        }
      }

      if(is.finite(determinant_check(newINV))== TRUE){
      error_x[k] = log(det(newINV))-psych::tr(newINV%*%stats::cov(test))
      }
    }
    press <- rbind(press, c(GraphicalLasso[gl], sum(error_x)/nfolds))
    sd_MSE <- c(sd_MSE, stats::sd(error_x)/sqrt(nfolds))
  }

  lowestPress <- max(press[,2], na.rm = TRUE)
  lowestplus1SE <- lowestPress + sd_MSE[which(press[,2] == lowestPress)]
  indexTuning <- which(abs(press[,2] - lowestplus1SE)==min(abs(press[,2] - lowestplus1SE), na.rm = TRUE))

  bestTuning = GraphicalLasso[indexTuning[1]]

  return(bestTuning)

}




