#' Estimate partial correlation network of common component
#'
#' The partial correlation and inverse covariance matrices of the common component(s) are estimated using graphical lasso.
#' Part of Sparse Network And Component (SNAC) analysis, see also [comSNAC].
#'
#' @param Pcommon Component loadings of the common component(s). Should be a matrix of p variables by Rc number of common components
#' @param Tcommon Component scores of the common component(s). Should be a matrix of n participants by Rc number of common components
#' @param lasso Lambda parameter for graphical lasso. Default is set at 0.002
#'
#' @return An inverse covariance matrix and partial correlation matrix of common component
#' @references Tio, P., Waldorp, L., & VanDeun, K. (2020). Constructing graphical models for multi-source data: Sparse network and component analysis. In \emph{Advanced Studies in Classification and Data Science} (pp. 275-287). Springer Singapore.
#' @export

netSNAC = function(Pcommon, Tcommon, lasso = 0.002){
  network = list()

  #Create data set based on common components
  Xcom <- Tcommon%*%t(Pcommon)

  #Perform Network analysis step using graphical lasso
  SigmaINVC_hat <- glasso::glasso(s = stats::cov(Xcom), rho = lasso, penalize.diagonal = FALSE)$wi

  #Take average of covariance
  newINV = matrix(nrow = ncol(Xcom), ncol = ncol(Xcom))

  for (i in 1:ncol(Xcom)){
    for (j in 1:ncol(Xcom)){
      newINV[i,j] = newINV[j,i] = mean(c(SigmaINVC_hat[i,j], SigmaINVC_hat[j,i]))
    }
  }

  #calculate partial correlations directly from inverse covariance matrix
  pcor = stats::cov2cor(newINV)*-1
  diag(pcor) = 1

  network$InvCov = SigmaINVC_hat
  network$Pcor = pcor
  return(network)

}
