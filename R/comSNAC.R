#' Decompose multi-domain data into common and domain-specific components.
#'
#' A component analysis technique that transforms a data set consisting of two data blocks (domains) into \code{R} components.
#' These components are either domain-specific (component loadings from only one domain) or common (component loadings from both domains); the latter of which can be used in [netSNAC].
#' Part of Sparse Network And Component (SNAC) analysis. For more information on the component analysis see Gu & Van Deun (2016).
#'
#' @param data Data file consisting of 2 data blocks with subjects as rows and variables as columns. Variables should be grouped per data block
#' @param b1 The number of variables in data block 1
#' @param b2 The number of variables in data block 2
#' @param R The number of components
#' @param method How is component structure determined? Default is "structured" for which a predetermined component structure in the form of a Target is required.
#' @param Target Required for structured analysis. Should be matrix of 2 rows (for each data block) by R components. Enter 0 if data block should not contribute to that component; 1 is data block should contribute to that component.
#'
#' @return Component loadings ($Pmartix) and scores ($Tmatrix) of all R components. Additionally returns the tuning parameter lambda (estimated using 10-fold cross validation) and the predetermined component structure.
#'
#' @references Gu, Z., & Van Deun, K. (2016). A variable selection method for simultaneous component based data integration. \emph{Chemometrics and Intelligent Laboratory Systems}, 158, 187-199.
#' @references Tio, P., Waldorp, L., & VanDeun, K. (2020). Constructing graphical models for multi-source data: Sparse network and component analysis. In \emph{Advanced Studies in Classification and Data Science} (pp. 275-287). Springer Singapore.
#' @export
comSNAC = function(data, b1, b2, R, method = "structured", Target){
  if(is.numeric(b1)==FALSE){
    stop("b1 is not a numeric entry.")
  }

  if(is.numeric(b2)==FALSE){
    stop("b2 is not a numeric entry.")
  }
  if(R==1){
    stop("R = 1 is not allowed.")
  }
  # if(missing(method)){
  #   stop("Choose method.")
  # }

  Jk <- c(b1, b2)
  component = list()

  if(method == 'structured'){
    if(missing(Target)){
      stop("Target is missing.")
    }
    if(missing(R)){
      stop("R is missing.")
    }
    if(ncol(Target)!=R){
      stop("Number of columns of Target must be equal to R.")
    }

    #Estimate lambda value for Component analysis' lasso using cross-validation
    para <- cv_structuredSCA(DATA = data, Jk = Jk, R = R, Target = Target)
    #Perform Component analysis step
    results = structuredSCA(DATA = data, Jk = Jk, R = R, Target = Target, LASSO = para$RecommendedLasso)

    #Undo shrinkage done on the P (loadings) and T (scores) matrices
    final <- undoShrinkage(DATA = data, R = R, Phat = results$Pmatrix)
    Pmatrix <- final$Pmatrix
    Tmatrix <- final$Tmatrix

    #Determine component structure
    num_common = function(Pmatrix, Jk){
      common <- NULL
      for (c in 1:ncol(Pmatrix)){ #look per component whether it is a common or distinctive one
        s <- 0
        d <- NULL
        for (j in 1:length(Jk)){
          d[j] <- sum(abs(Pmatrix[(s+1):(s+Jk[j]),c]))
          if(d[j]!=0) {d[j] = 1}
          s <- sum(Jk[1:j])
        }
        if (sum(d)>= 2) {common[c]= "common"} else if (d[1]==1 & d[2]==0) {common[c] = "block 1"} else if (d[1]==0 & d[2]==1) {common[c] = "block 2"} else {common[c] = "NAN"}
        #common component is defined as having variables of at least two blocks
      }

      return(common)
    }

    structure <- num_common(Pmatrix, Jk)

    component$lambda_com <- para$RecommendedLasso
    component$Structure <- structure
  }

  # if(method == 'unstructured'){
  #   #Estimate lambda values for Component analysis' lasso and group lasso using cross-validation
  #   para <- cv_sparseSCA(DATA = data, Jk = Jk, R = R)
  #   #Perform Component analysis step
  #   results <- sparseSCA(DATA = data, Jk = Jk, R = R, LASSO = para$RecommendedLambda[1],
  #                        GROUPLASSO = para$RecommendedLambda[2])
  #
  #   #Undo shrinkage done on the P (loadings) and T (scores) matrices
  #   final <- undoShrinkage(DATA = data, R = R, Phat = results$Pmatrix)
  #   Pmatrix <- final$Pmatrix
  #   Tmatrix <- final$Tmatrix
  #
  #   #Determine component structure
  #   structure <- num_common(Pmatrix, Jk)
  #
  #   component$Parameters <- para$RecommendedLambda
  #   component$Structure <- structure
  # }

  component$Pmatrix <- Pmatrix
  component$Tmatrix <- Tmatrix
  return(component)

}
