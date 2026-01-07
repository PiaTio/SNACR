#' Check whether a determinant of an inverse covariance matrix can be determined without warning or error
#'
#' @param x Inverse covariance matrix
#'
#' @returns TRUE if determinant can be determined without warnings or error
#' @export
#'
determinant_check <- function(x){
  tryCatch(
    expr = {
      suppressWarnings(log(det(x)))
    },
    error = function(e){},
    warning = function(w){},
    finally = {}
  )
}
