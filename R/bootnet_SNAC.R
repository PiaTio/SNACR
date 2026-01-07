#' Bootstrapped SNAC common network estimation
#'
#' Adjusted version of function bootnet() from R package bootnet (Epskamp, Borsboom, & Fried, 2018)
#'
#' @param data Data set
#' @param nBoots Number of bootstraps
#' @param type Bootstrap method to use, default to "nonparametric"
#' @param nCores Number of cores used, default is 1
#' @param statistics Type of statistics calculated, defaults to "edge"
#' @param fun Function used to estimate network
#' @param verbose  Messages on what is being done?
#' @param labels Variable labels
#' @param alpha centrality alpha
#' @param caseMin Minimum proportion to drop
#' @param caseMax Maximum proportion to drop
#' @param caseN Default is 10
#' @param subCases subCases
#' @param computeCentrality Whether to compute centrality TRUE/FALSE
#' @param propBoot M out of N
#' @param replacement Bootstrap with replacement? Default is TRUE
#' @param includeDiagonal Include diagonal? Default is FALSE
#' @param bridgeArgs Argument for bridge statistics
#' @param library Path to save files
#' @param memorysaver Default is TRUE
#' @param ...
#'
#' @returns A bootnet object with the following elements
#' \item{sampleTable}{A data frame containing all estimated values on the real sample.}
#' \item{bootTable}{A data frame containing all estimated values on the bootstrapped samples.}
#' \item{sample}{A bootnetResult object with plot and print method containing the estimated network of the real sample.}
#' \item{boots}{A list of bootnetResult objects containing the raw bootstrap results.}
#' \item{numbCommon}{Number of times the component comp was a common component.}
#' \item{lambdaCom}{String of bootstrapped 10-fold cross-validation values for lasso lambda.}
#' \item{lambdaNet}{String of bootstrapped 10-fold cross-validation values for graphical lasso lambda.}
#' @references Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating psychological networks and their accuracy: A tutorial paper. Behavior Research Methods, 50(1), 195-212.
#' @export
#'
bootnet_SNAC <- function(
    data, # Dataset
    nBoots = 1000, # Number of bootstrap samples.
    type = "nonparametric", # Bootstrap method to use
    nCores = 1,
    statistics = "edge",
    fun,
    # prepFun, # Fun to produce the correlation or covariance matrix
    # prepArgs, # list with arguments for the correlation function
    # estFun, # function that results in a network
    # estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
    # graphFun, # set to identity if missing
    # graphArgs, # Set to null if missing
    # intFun, # Set to null if missing
    # intArgs, # Set to null if missing
    verbose = TRUE, # messages on what is being done?
    # construct = c("default","function","arguments"),
    labels, # if missing taken from colnames
    alpha = 1, # centrality alpha
    caseMin = 0.05, # minimum proportion to DROP
    caseMax = 0.75, # Maximum proportion to DROP
    caseN = 10,
    subCases,
    computeCentrality = TRUE,
    propBoot = 1, # M out of N
    # subsampleSize,
    replacement = TRUE,
    includeDiagonal = FALSE,
    bridgeArgs=list(),
    library = .libPaths(),
    memorysaver = TRUE,
    # datatype = c("normal","graphicalVAR"), # Extracted from object or given
    ... # Other arguments
    # edgeResample = FALSE # If true, only resample edges from original estimate
    # scaleAdjust = FALSE
){

  # Set a NULL sampleResult:
  sampleResult <- NULL
  type <- match.arg(type)
  manual <- FALSE
  datatype <- "normal"
  dots <- list(...)
  N <- ncol(data)
  Np <- nrow(data)
  if (missing(fun)){
    fun <- NULL
  }

  default = "none"

  inputCheck <- checkInput(
    default = default,
    fun = fun,
    # prepFun = prepFun, # Fun to produce the correlation or covariance matrix
    # prepArgs = prepArgs, # list with arguments for the correlation function
    # estFun=estFun, # function that results in a network
    # estArgs=estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
    # graphFun=graphFun, # set to identity if missing
    # graphArgs=graphArgs, # Set to null if missing
    # intFun=intFun, # Set to null if missing
    # intArgs=intArgs, # Set to null if missing
    # sampleSize = Np,
    # construct = construct,
    .dots = dots
  )
  weighted <- TRUE
  signed <- TRUE
  directed <- FALSE


  # First test if data is a data frame:
  if (datatype == "normal" && !manual && !(is.data.frame(data) || is.matrix(data))){
    stop("'data' argument must be a data frame")
  }

  # If matrix coerce to data frame:
  if (!manual && is.matrix(data)){
    data <- as.data.frame(data)
  }

  if (missing(labels)){
    if (manual){
      labels <- colnames(graph)
      if (is.null(labels)){
        labels <- seq_len(ncol(graph))
      }
    } else {
      labels <- colnames(data)
      if (is.null(labels)){
        labels <- seq_len(ncol(data))
      }
    }


  }

  # Estimate sample result:
  # Check the input:


  if (!manual)
  {


    if (is.null(sampleResult)){
      if (verbose){
        message("Estimating sample network...")
      }

      sampleResult <- estimateNetwork_SNAC(data,
                                      default = default,
                                      fun = inputCheck$estimator,
                                      .dots = inputCheck$arguments,
                                      labels = labels,
                                      verbose = verbose,
                                      weighted = weighted,
                                      signed = signed,
                                      .input = inputCheck,
                                      datatype = datatype,
                                      directed = directed)
    }

  } else {

    sampleResult <- list(
      graph = graph,
      intercepts = intercepts,
      labels = labels,
      nNodes = N,
      nPerson = Np,
      estimator = inputCheck$estimator,
      arguments = inputCheck$arguments,
      default = default,
      weighted = weighted,
      signed = signed
    )
    class(sampleResult) <- c("bootnetResult", "list")

  }

  # Extract arguments:
  # default <- sampleResult$input$default
  # prepFun <- sampleResult$input$prepFun
  # prepArgs <- sampleResult$input$prepArgs
  # estFun <- sampleResult$input$estFun
  # estArgs <- sampleResult$input$estArgs
  # graphFun <- sampleResult$input$graphFun
  # graphArgs <- sampleResult$input$graphArgs
  # intFun <- sampleResult$input$intFun
  # intArgs <- sampleResult$input$intArgs


  # if (!isSymmetric(as.matrix(sampleResult[['graph']]))){
  #   stop("bootnet does not support directed graphs")
  # }


  #   ### Observation-wise bootstrapping!
  #   if (type == "observation"){
  # Bootstrap results:
  if (nCores == 1){
    bootResults <- vector("list", nBoots)

    if (verbose){
      message("Bootstrapping...")
      pb <- utils::txtProgressBar(0,nBoots,style = 3)
    }

    for (b in seq_len(nBoots)){

      tryLimit <- 10
      tryCount <- 0
      repeat{

        if (! type %in% c("node","person")){
          nNode <- N
          inSample <- seq_len(N)

         {
            nPerson <- Np

            if (datatype == "normal"){
              bootData <- data[sample(seq_len(Np), round(propBoot*Np), replace=replacement), ]
            } else {
              bootData <- data
              bootSample <- sample(seq_len(Np), round(propBoot*Np), replace=replacement)
              bootData$data_c <- bootData$data_c[bootSample,]
              bootData$data_l <- bootData$data_l[bootSample,]
            }

          }

        } else {
          # Personwise:
          if (length(subCases) == 1){
            nPerson <- subCases
          } else {
            nPerson <- sample(subCases,1)
          }

          inSample <- 1:N
          persSample <- sort(sample(seq_len(Np),nPerson))
          if (datatype == "normal"){
            bootData <- data[persSample,, drop=FALSE]
          } else if (datatype == "graphicalVAR") {
            bootData <- data
            bootData$data_c <- bootData$data_c[persSample,, drop = FALSE]
            bootData$data_l <- bootData$data_l[persSample,, drop = FALSE]
          }
        }

        # Some checks to remove progress bars:
        # if (!missing(prepFun)){
        #   # EBICglasso:
        #   if (!missing(prepArgs) & is.list(prepArgs) & identical(prepFun,qgraph::cor_auto)){
        #     prepArgs$verbose <- FALSE
        #   }
        # }
        #
        res <- suppressWarnings(try({
          # estimateNetwork_SNAC(bootData,
          #                 default = default,
          #                 fun = fun,
          #                 prepFun = prepFun, # Fun to produce the correlation or covariance matrix
          #                 prepArgs = prepArgs, # list with arguments for the correlation function
          #                 estFun = estFun, # function that results in a network
          #                 estArgs = estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
          #                 graphFun = graphFun, # set to identity if missing
          #                 graphArgs = graphArgs, # Set to null if missing
          #                 intFun = intFun, # Set to null if missing
          #                 intArgs = intArgs, # Set to null if missing
          #                 labels = labels[inSample],
          #                 verbose = FALSE,
          #                 construct = construct)
          #

          estimateNetwork_SNAC(bootData,
                          default = default,
                          fun = inputCheck$estimator,
                          .dots = inputCheck$arguments,
                          labels = labels[inSample],
                          verbose = FALSE,
                          weighted = weighted,
                          signed = signed,
                          .input = inputCheck,
                          memorysaver = memorysaver,
                          directed = directed)

        }))
        if (methods::is(res, "try-error")){

          if (tryCount == tryLimit) {
            stop("Maximum number of errors in bootstraps reached")
          }

          # warning("Error in bootstrap; retrying")
          tryCount <- tryCount + 1
        } else {
          break
        }

      }

      bootResults[[b]] <- res

      if (verbose){
        utils::setTxtProgressBar(pb, b)
      }
    }
    if (verbose){
      close(pb)
    }
  } else {
    if (verbose){
      message("Bootstrapping...")
    }

    if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) &&
        Sys.info()["sysname"] == "Darwin" && gsub("\\..*","",getRversion()) == "4") {
      snow::setDefaultClusterOptions(setup_strategy = "sequential")
    }

    nClust <- nCores - 1
    cl <- snow::makeSOCKcluster(nClust)

    # IF graph or data is missing, dummy graph:
    if (missing(graph)){
      graph <- matrix(0,N,N)
    }
    if (missing(data)){
      data <- matrix(0,Np,N)
    }
    if (missing(intercepts)){
      intercepts <- rep(0,N)
    }
    if (missing(sampleSize)){
      sampleSize <- Np
    }

    # Needed arguments to be excluded:
    excl <- c("prepFun", "prepArgs", "estFun", "estArgs", "graphFun",
              "graphArgs", "intFun", "intArgs", "fun")

    snow::clusterExport(cl, ls()[!ls()%in%c(excl,"cl")], envir = environment())
    # clusterExport(cl, export, envir = environment())


    # Run loop:
    bootResults <- pblapply::pblapply(seq_len(nBoots), function(b){
      # Set library:
      .libPaths(library)

      tryLimit <- 10
      tryCount <- 0
      repeat{

        if (! type %in% c("node","person")){
          nNode <- N
          inSample <- seq_len(N)

          if (type == "jackknife"){
            if (datatype == "normal"){
              bootData <- data[-b,,drop=FALSE]
            } else {
              bootData <- data
              bootData$data_c <- bootData$data_c[-b,,drop=FALSE]
              bootData$data_l <- bootData$data_l[-b,,drop=FALSE]
            }

            nPerson <- Np - 1
          } else {
            nPerson <- Np

            if (datatype == "normal"){
              bootData <- data[sample(seq_len(Np), round(propBoot*Np), replace=replacement), ]
            } else {
              bootData <- data
              bootSample <- sample(seq_len(Np), round(propBoot*Np), replace=replacement)
              bootData$data_c <- bootData$data_c[bootSample,]
              bootData$data_l <- bootData$data_l[bootSample,]
            }

          }

        } else if (type == "node") {

          # Nodewise
          nPerson <- Np
          nNode <- sample(subNodes,1)
          inSample <- sort(sample(seq_len(N),nNode))
          if (datatype == "normal"){
            bootData <- data[,inSample, drop=FALSE]
          } else if (datatype == "graphicalVAR") {
            bootData <- data
            bootData$data_c <- bootData$data_c[,data$vars[inSample], drop = FALSE]
            bootData$data_l <- bootData$data_l[,c("1",grep(data$vars[inSample],names(data$data_l),value=TRUE)), drop = FALSE]
          }

        } else {
          # Personwise:
          if (length(subCases) == 1){
            nPerson <- subCases
          } else {
            nPerson <- sample(subCases,1)
          }
          inSample <- 1:N
          persSample <- sort(sample(seq_len(Np),nPerson))
          if (datatype == "normal"){
            bootData <- data[persSample,, drop=FALSE]
          } else if (datatype == "graphicalVAR") {
            bootData <- data
            bootData$data_c <- bootData$data_c[persSample,, drop = FALSE]
            bootData$data_l <- bootData$data_l[persSample,, drop = FALSE]
          }
        }

        # Some checks to remove progress bars:
        # if (!missing(prepFun)){
        #   # EBICglasso:
        #   if (!missing(prepArgs) & is.list(prepArgs) & identical(prepFun,qgraph::cor_auto)){
        #     prepArgs$verbose <- FALSE
        #   }
        # }

        res <- suppressWarnings(try({
          # estimateNetwork_SNAC(bootData,
          #                 default = default,
          #                 fun = fun,
          #                 prepFun = prepFun, # Fun to produce the correlation or covariance matrix
          #                 prepArgs = prepArgs, # list with arguments for the correlation function
          #                 estFun = estFun, # function that results in a network
          #                 estArgs = estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
          #                 graphFun = graphFun, # set to identity if missing
          #                 graphArgs = graphArgs, # Set to null if missing
          #                 intFun = intFun, # Set to null if missing
          #                 intArgs = intArgs, # Set to null if missing
          #                 labels = labels[inSample],
          #                 verbose = FALSE,
          #                 construct = construct)
          #

          estimateNetwork_SNAC(bootData,
                          default = default,
                          fun = inputCheck$estimator,
                          .dots = inputCheck$arguments,
                          labels = labels[inSample],
                          verbose = FALSE,
                          weighted = weighted,
                          signed = signed,
                          .input = inputCheck,
                          memorysaver = memorysaver,
                          directed = directed)

        }))
        if (is(res, "try-error")){

          if (tryCount == tryLimit) {
            stop("Maximum number of errors in bootstraps reached")
          }

          # warning("Error in bootstrap; retrying")
          tryCount <- tryCount + 1
        } else {
          break
        }

      }

      return(res)
    }, cl = cl)
  }

  ### Compute the full parameter table!!
  if (verbose){
    message("Computing statistics...")
  }
  statTableOrig <- statTable(sampleResult,  name = "sample", alpha = alpha, computeCentrality = computeCentrality,statistics=statistics, directed=directed, includeDiagonal=includeDiagonal, bridgeArgs=bridgeArgs)

  if (nCores == 1){
    if (verbose){
      pb <- txtProgressBar(0,nBoots,style = 3)
    }
    statTableBoots <- vector("list", nBoots)
    for (b in seq_len(nBoots)){
      statTableBoots[[b]] <- statTable(bootResults[[b]], name = paste("boot",b), alpha = alpha, computeCentrality = computeCentrality, statistics=statistics, directed=directed,  bridgeArgs=bridgeArgs, includeDiagonal=includeDiagonal)
      if (verbose){
        setTxtProgressBar(pb, b)
      }
    }
    if (verbose){
      close(pb)
    }
  }  else {
    statTableBoots <- pblapply(seq_len(nBoots),function(b){
      # Set library:
      .libPaths(library)

      statTable(bootResults[[b]], name = paste("boot",b), alpha = alpha, computeCentrality = computeCentrality, statistics=statistics, directed=directed, bridgeArgs=bridgeArgs, includeDiagonal=includeDiagonal)
    }, cl = cl)
    # Stop the cluster:
    snow::stopCluster(cl)
  }

  lambdaCom = NULL
  lambdaNet = NULL
  numbCommon = 0
  for (i in 1:nBoots){
    lambdaCom[i] = bootResults[[i]]$cv_lambda_com
    lambdaNet[i] = bootResults[[i]]$cv_lambda_net
    numbCommon = numbCommon + (bootResults[[i]]$no_common)

  }

  # Ordereing by node name to make nice paths:
  Result <- list(
    sampleTable = dplr::ungroup(statTableOrig),
    bootTable =  ungroup(dplyr::bind_rows(statTableBoots)),
    sample = sampleResult,
    boots = bootResults,
    type = type,
    sampleSize = Np,
    numbCommon = numbCommon,
    lambdaCom = lambdaCom,
    lambdaNet = lambdaNet)

  class(Result) <- "bootnet"

  return(Result)

}
