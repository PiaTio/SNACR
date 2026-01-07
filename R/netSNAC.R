#' Estimate partial correlation network of common component
#'
#' The partial correlation and inverse covariance matrices of the common component(s) are estimated using graphical lasso.
#' Part of Sparse Network And Component (SNAC) analysis, see also [comSNAC].
#'
#' @param Pcommon Component loadings of the common component(s). Should be a matrix of p variables by Rc number of common components
#' @param Tcommon Component scores of the common component(s). Should be a matrix of n participants by Rc number of common components
#' @param lambda_net Range of the lambda tuning parameter for the graphical lasso. Default value is a sequence of 100 numbers from 0.0000001 to 0.05. The numbers are equally spaced on a log scale.
#'
#' @return The partial correlation matrix and inverse covariance matrix of the common component, as well as the
#' \item{graph}{Partial correlation matrix.}
#' \item{inv_cov}{Inverse covariance matrix.}
#' \item{lambda_net}{Lambda tuning parameter used in for graphical lasso that leads to the model with PRESS closest to the lowest PRESS + 1SE.}
#' @references Tio, P., Waldorp, L., & VanDeun, K. (2020). Constructing graphical models for multi-source data: Sparse network and component analysis. In \emph{Advanced Studies in Classification and Data Science} (pp. 275-287). Springer Singapore.
#' @export

netSNAC = function(Pcommon, Tcommon, lambda_net){

  #Create data set based on chosen components
  Xcomp <- Tcommon%*%t(Pcommon)
  #Add ridge penalty for increased numeric stability
  ridge = matrix(0, nrow = ncol(Xcomp), ncol = ncol(Xcomp))
  diag(ridge) = 1
  cov_ridge = stats::cov(Xcomp) + (ridge*1/100000)

  if(missing(lambda_net)){
    bestLambda = cv_graphicalLasso(Xcomp)
    SigmaINVC_hat = glasso::glasso(cov_ridge, rho = bestLambda, penalize.diagonal = FALSE)$wi #estimate inverse covariance matrix
    lambda_net <- bestLambda
  }
  if(length(lambda_net)==1){
    SigmaINVC_hat = glasso::glasso(cov_ridge, rho = lambda_net, penalize.diagonal = FALSE)$wi #estimate inverse covariance matrix
  }
  if(length(lambda_net)>1){
    bestLambda = cv_graphicalLasso(Xcomp, GraphicalLasso = lambda_net)
    SigmaINVC_hat = glasso::glasso(cov_ridge, rho = bestLambda, penalize.diagonal = FALSE)$wi #estimate inverse covariance matrix
    lambda_net <- bestLambda
  }


  #Take average of covariance
  newINV = matrix(nrow = ncol(SigmaINVC_hat), ncol = ncol(SigmaINVC_hat))

  for (i in 1:ncol(SigmaINVC_hat)){
    for (j in 1:ncol(SigmaINVC_hat)){
      newINV[i,j] = newINV[j,i] = mean(c(SigmaINVC_hat[i,j], SigmaINVC_hat[j,i]))
    }
  }

  #calculate partial correlations directly from inverse covariance matrix
  pcor = stats::cov2cor(newINV)*-1
  diag(pcor) = 1

  Results = list(
    graph = pcor,
    inv_cov = newINV,
    lambda_net = lambda_net
  )
return(Results)

}
