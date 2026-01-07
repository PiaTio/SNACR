#' Determine the component structure after the component step of Sparse Network And Component analysis
#'
#' @param Pmatrix P matrix of the component
#' @param Jk A vector with the number of columns of each data block
#'
#' @returns Table with for each component whether they are a common, block 1-specific, or block 2-specific component.
#' @export
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
