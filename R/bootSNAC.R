#' Single function to perform two-step Sparse Network And Component (SNAC) analysis to estimate cross-domain partial correlations
#'
#'This function is used to bootstrap SNAC analysis and should not be used on its own.
#'
#' @param data Data file consisting of 2 data blocks with subjects as rows and variables as columns. Variables should be grouped per data block
#' @param b1 The number of variables in data block 1
#' @param b2 The number of variables in data block 2
#' @param R The number of components
#' @param method How is component structure determined? Default is "structured" for which a predetermined component structure in the form of a Target is required.
#' @param Target Required for structured analysis. Should be matrix of 2 rows (for each data block) by R components. Enter 0 if data block should not contribute to that component; 1 is data block should contribute to that component.
#' @param comp The position of the component to bootstrap
#' @param lambda_com Lambda parameter for lasso (component step). Default is determined by 10-fold cross-validation.
#' @param lambda_net Lambda parameter for graphical lasso (network step). Default is determined by 10-fold cross-validation.
#'
#' @returns
#' \item{graph}{A partial correlation matrix}
#' \item{inv_cov}{Inverse covariance matrix.}
#' \item{cv_lambda_com}{Estimated value for lasso lambda.}
#' \item{cv_lambda_net}{Estimated value for graphical lasso lambda.}
#' \item{no_common}{Value indicating whether the component was a common component (1) or not (0).}
#' @export
#'
bootSNAC <- function(data, b1, b2, R, method = "structured", Target, comp, lambda_com, lambda_net){
  no_common = NULL
  Jk <- c(b1, b2)
  component = list()

  if(missing(lambda_com)){
    #Estimate lambda value for Component analysis' lasso using cross-validation
    para <- cv_structuredSCA(DATA = data, Jk = Jk, R = R, Target = Target)
    #Perform Component analysis step
    results = structuredSCA(DATA = data, Jk = Jk, R = R, Target = Target, LASSO = para$RecommendedLasso)
    cv_lambda_com <- para$RecommendedLasso

  } else {
    results = structuredSCA(DATA = data, Jk = Jk, R = R, Target = Target, LASSO = lambda_com)
    cv_lambda_com <- lambda_com
  }
  #Undo shrinkage done on the P (loadings) and T (scores) matrices
  final <- undoShrinkage(DATA = data, R = R, Phat = results$Pmatrix)
  Pmatrix <- final$Pmatrix
  Tmatrix <- final$Tmatrix

  #Determine component structure
  structure <- num_common(Pmatrix, Jk)

  component$Structure <- structure
  component$Pmatrix <- Pmatrix
  component$Tmatrix <- Tmatrix
  #
  # Tar <- NULL
  # for (t in 1:ncol(Target)){
  #   d = Target[,t]
  #   if (sum(d)>= 2) {Tar[t]= "common"} else if (d[1]==1 & d[2]==0) {Tar[t] = "block 1"} else if (d[1]==0 & d[2]==1) {Tar[t] = "block 2"} else {Tar[t] = "NAN"}
  # }

  #Check if the estimated component structure contains common component in indicated spot
  if(sum(component$Structure[comp] == "common") == 1){
    no_common = 1
  }

  #Create data set based on chosen components
  Xcomp <- component$Tmatrix[,comp]%*%t(component$Pmatrix[,comp])
  #Add ridge penalty for increased numeric stability
  ridge = matrix(0, nrow = ncol(Xcomp), ncol = ncol(Xcomp))
  diag(ridge) = 1
  cov_ridge = stats::cov(Xcomp) + (ridge*1/100000)

  if(missing(lambda_net)){
    bestLambda = cv_graphicalLasso(Xcomp)
    SigmaINVC_hat = glasso::glasso(cov_ridge, rho = bestLambda, penalize.diagonal = FALSE)$wi #estimate inverse covariance matrix
    cv_lambda_net <- bestLambda
  } else {
    SigmaINVC_hat = glasso::glasso(cov_ridge, rho = lambda_net, penalize.diagonal = FALSE)$wi #estimate inverse covariance matrix
    cv_lambda_net <- lambda_net
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

  cvResults = list(
    graph = pcor,
    inv_cov = newINV,
    cv_lambda_com = cv_lambda_com,
    cv_lambda_net = cv_lambda_net,
    no_common = no_common
  )
  return(cvResults)
}
