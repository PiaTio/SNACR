#' Summarize additional bootstrapped SNAC network results
#'
#' @param Bootstraps A bootnet object containing the bootstrapped network estimation of a SNAC network
#'
#' @returns
#' \item{cocorec}{The common component recovery rate. How often is the component a common component?}
#' \item{lambdaCom}{The values for the lasso's lambda parameter (SNAC's component step)}
#' \item{lambdaNet}{The values for the graphical lasso's lambda parameter (SNAC's network step)}
#' \item{Structure}{A list of the component structures for each bootstrap result}
#' @export
#'
bootSNAC_sum <- function(Bootstraps){
  nBoots = length(Bootstraps$boots)
  lambdaCom = NULL
  lambdaNet = NULL
  compStruc = list()
  numbCommon = 0
  for (i in 1:nBoots){
    lambdaCom[i] = Bootstraps$boots[[i]]$results$cv_lambda_com
    lambdaNet[i] = Bootstraps$boots[[i]]$results$cv_lambda_net
    numbCommon = numbCommon + (Bootstraps$boots[[i]]$results$no_common)
    compStruc[[i]] = Bootstraps$boots[[i]]$result$com_struc
  }
  cocorec = (numbCommon/nBoots*100)

  # Ordereing by node name to make nice paths:
  Result <- list(
    cocorec = cocorec,
    lambdaCom = lambdaCom,
    lambdaNet = lambdaNet,
    Structure = compStruc)
  message(paste("The common component recovery rate is", paste0(cocorec, collapse=", "), "%"))

  return(Result)
}
