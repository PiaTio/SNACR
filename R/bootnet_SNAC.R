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
#' @param subNodes subNodes
#' @param computeCentrality Whether to compute centrality TRUE/FALSE
#' @param propBoot M out of N
#' @param replacement Bootstrap with replacement? Default is TRUE
#' @param includeDiagonal Include diagonal? Default is FALSE
#' @param bridgeArgs Argument for bridge statistics
#' @param library Path to save files
#' @param memorysaver Default is TRUE
#' @param ... Argument for 'fun' function
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
    subNodes,
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
  # subNodes and subCases:
  if (missing(subNodes)){
    if (datatype == "normal"){
      subNodes <- 2:(N-1)
    }
  }

  checkInput <- function(
    default = c("none", "EBICglasso","ggmModSelect", "pcor","IsingFit","IsingSampler", "huge","adalasso","mgm","relimp",
                "cor","TMFG","ggmModSelect","LoGo","graphicalVAR","piecewiseIsing","SVAR_lavaan",
                "GGMncv"),
    fun, # Estimator function
    # prepFun, # Fun to produce the correlation or covariance matrix
    # prepArgs, # list with arguments for the correlation function
    # estFun, # function that results in a network
    # estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
    # graphFun, # set to identity if missing
    # graphArgs, # Set to null if missing
    # intFun, # Set to null if missing
    # intArgs, # Set to null if missing
    # nSample,
    verbose=TRUE,
    # construct = c("default","function","arguments"),
    .dots = list(),
    ... # Arguments to the estimator function
  ){
    construct <- "function"
    if (default[[1]]=="glasso") default <- "EBICglasso"
    if (default[[1]]=="IsingSampler") default <- "IsingSampler"
    default <- match.arg(default)
    # construct <- match.arg(construct)

    ### DEFAULT OPTIONS ###
    if (missing(fun)){
      fun <- NULL
    }



    # Stop if not compatible:
    dots <- c(.dots,list(...))

    # gather names:
    argNames <- character(0)
    #
    # if (!missing(prepFun)){
    #   argNames <- c(argNames,"prepFun")
    # }
    # if (!missing(prepArgs)){
    #   argNames <- c(argNames,"prepArgs")
    # }
    # if (!missing(estFun)){
    #   argNames <- c(argNames,"estFun")
    # }
    # if (!missing(estArgs)){
    #   argNames <- c(argNames,"estArgs")
    # }
    # if (!missing(graphFun)){
    #   argNames <- c(argNames,"graphFun")
    # }
    # if (!missing(graphArgs)){
    #   argNames <- c(argNames,"graphArgs")
    # }
    # if (!missing(intFun)){
    #   argNames <- c(argNames,"intFun")
    # }
    # if (!missing(intArgs)){
    #   argNames <- c(argNames,"intArgs")
    # }
    #
    # # Not compatible if construct is used:
    # if (length(dots) > 0 && construct == "arguments"){
    #
    #   stop(paste0("Ambiguous argument specification. Old functonality is used (construct = 'arguments') in combination with new functionality arguments (implying construct = 'function'): ",
    #               paste0("'",names(dots),"'",collapse="; "),". These arguments are NOT compatible!"))
    #
    # }
    #
    # # relimp not compatable with old:
    # if (construct == "arguments" & default == "relimp"){
    #   stop("default = 'relimp' not supported with old bootnet style (construct = 'arguments')")
    #
    # }
    #
    # if (length(argNames) > 0 && construct == "function"){
    #
    #   stop(paste0("Ambiguous argument specification. New functonality is used (construct = 'function') in combination with old functionality arguments (implying construct = 'arguments'): ",
    #               paste0("'",argNames,"'",collapse="; "),". These arguments are NOT compatible!"))
    #
    # }
    #
    # # not compatible if both dots are used and arguments are used:
    # if (length(argNames) > 0 & length(dots) > 0){
    #
    #   stop(paste0("Ambiguous argument specification. Both old functionality arguments are used, compatible with construct = 'arguments': ",
    #               paste0("'",argNames,"'",collapse="; "),", as well as new functionality arguments are used, compatible with construct = 'function': ",
    #               paste0("'",names(dots),"'",collapse="; "),". These two types of arguments are NOT compatible!"))
    #
    # }
    #
    #
    # # Check to construct via function or to construct via arguments:
    # # if no default and no fun, use arguments:
    # if (construct == "default"){
    #   construct <- "function"
    #
    #   if (default == "none" && is.null(fun)){
    #     construct <- "arguments"
    #   }
    #
    #   # If fun is missing, default is not none and one argument is not missing, use arguments (backward competability):
    #   if (default != "none" && is.null(fun) && (!missing(prepFun) | !missing(prepArgs) | !missing(estFun) | !missing(estArgs))){
    #     construct <- "arguments"
    #   }
    # }
    #
    # # Check if arguments are not missing:
    # if (default == "none" && construct == "arguments"){
    #   if (missing(prepFun) | missing(prepArgs) | missing(estFun) | missing(estArgs)){
    #     stop("If 'default' is not set and 'fun' is missing, 'prepFun', 'prepArgs', 'estFun' and 'estArgs' may not be missing.")
    #   }
    # }

    ### Construct estimator function via function:
    if (construct == "function"){
      # Arguments:
      Args <- dots
      #
      # # Warn user that arguments are ignored:
      # if (!missing(prepFun)){
      #   warning("'prepFun' argument is ignored as a function is used as arguments. To use 'prepFun', please set construct = 'arguments'")
      # }
      # if (!missing(prepArgs)){
      #   warning("'prepArgs' argument is ignored as a function is used as arguments. To use 'prepArgs', please set construct = 'arguments'")
      # }
      # if (!missing(estFun)){
      #   warning("'estFun' argument is ignored as a function is used as arguments. To use 'estFun', please set construct = 'arguments'")
      # }
      # if (!missing(estArgs)){
      #   warning("'estArgs' argument is ignored as a function is used as arguments. To use 'estArgs', please set construct = 'arguments'")
      # }
      # if (!missing(graphFun)){
      #   warning("'graphFun' argument is ignored as a function is used as arguments. To use 'graphFun', please set construct = 'arguments'")
      # }
      # if (!missing(graphArgs)){
      #   warning("'graphArgs' argument is ignored as a function is used as arguments. To use 'graphArgs', please set construct = 'arguments'")
      # }
      # if (!missing(intFun)){
      #   warning("'intFun' argument is ignored as a function is used as arguments. To use 'intFun', please set construct = 'arguments'")
      # }
      # if (!missing(intArgs)){
      #   warning("'intArgs' argument is ignored as a function is used as arguments. To use 'intArgs', please set construct = 'arguments'")
      # }
      #
      # per default:
      if (default == "none"){
        Function <- fun
      } else if (default == "EBICglasso"){
        Function <- bootnet::bootnet_EBICglasso
      } else if (default == "ggmModSelect"){
        Function <- bootnet::bootnet_ggmModSelect
      } else if (default == "IsingFit"){
        Function <- bootnet::bootnet_IsingFit
      } else if (default == "IsingSampler"){
        Function <- bootnet::bootnet_IsingSampler
      } else if (default == "pcor"){
        Function <- bootnet::bootnet_pcor
      } else if (default == "cor"){
        Function <- bootnet::bootnet_cor
      } else if (default == "adalasso"){
        Function <- bootnet::bootnet_adalasso
      } else if (default == "huge"){
        Function <- bootnet::bootnet_huge
      } else if (default == "mgm"){
        Function <- bootnet::bootnet_mgm
      } else if (default == "relimp"){
        Function <- bootnet::bootnet_relimp
      } else if (default == "TMFG"){
        Function <- bootnet::bootnet_TMFG
      } else if (default == "LoGo"){
        Function <- bootnet::bootnet_LoGo
      } else if (default == "graphicalVAR"){
        Function <- bootnet::bootnet_graphicalVAR
      } else if (default == "piecewiseIsing"){
        Function <- bootnet::bootnet_piecewiseIsing
      } else if (default == "SVAR_lavaan"){
        Function <- bootnet::bootnet_SVAR_lavaan
      } else stop("Currently not supported.")



      # } else {
      #   warning("Arguments (prepFun, estFun, etcetera) used to construct estimator. This functionality is deprecated and will no longer be supported in a future version of bootnet. Please consult the manual or contact the authors.")
      #
      #   # Check dots, and warn user:
      #   if (length(dots) > 0){
      #     dotNames <- names(dots)
      #     warning(paste0("Arguments (prepFun, estFun, etcetera) used to construct estimator. As a result, the following arguments are ignored: ",paste0("'",dotNames,"'", collapse = ", "),". To use these arguments use construct = 'function' and supply a default set or set the 'fun' argument. In addition, do not use the 'prepFun', 'estFun', etcetera arguments."))
      #   }
      #
      #   # Construct via arguments
      #   if (!(default == "none")){
      #     # prepFun:
      #     if (missing(prepFun)){
      #       prepFun <- switch(default,
      #                         EBICglasso = qgraph::cor_auto,
      #                         IsingFit = binarize,
      #                         IsingSampler = binarize,
      #                         pcor = qgraph::cor_auto,
      #                         huge = function(x)huge::huge.npn(na.omit(as.matrix(x)),verbose = FALSE),
      #                         adalasso = identity
      #       )
      #       #       prepFun <- switch(default,
      #       #                         EBICglasso = cor,
      #       #                         IsingFit = binarize,
      #       #                         pcor = cor
      #       #       )
      #     }
      #
      #     # prepArgs:
      #     #     qgraphVersion <- packageDescription("qgraph")$Version
      #     #     qgraphVersion <- as.numeric(strsplit(qgraphVersion,split="\\.|\\-")[[1]])
      #     #     if (length(qgraphVersion)==1) qgraphVersion <- c(qgraphVersion,0)
      #     #     if (length(qgraphVersion)==2) qgraphVersion <- c(qgraphVersion,0)
      #     #     goodVersion <-
      #     #       (qgraphVersion[[1]] >= 1 & qgraphVersion[[2]] >= 3 & qgraphVersion[[3]] >= 1) |
      #     #       (qgraphVersion[[1]] >= 1 & qgraphVersion[[2]] > 3) |
      #     #       qgraphVersion[[1]] > 1
      #
      #     if (missing(prepArgs)){
      #       prepArgs <- switch(default,
      #                          EBICglasso = ifElse(identical(prepFun,qgraph::cor_auto),list(verbose=verbose),
      #                                              ifElse(identical(prepFun,cor),list(use = "pairwise.complete.obs"),list())),
      #                          IsingFit = list(),
      #                          pcor =  ifElse(identical(prepFun,qgraph::cor_auto),list(verbose=verbose),
      #                                         ifElse(identical(prepFun,cor),list(use = "pairwise.complete.obs"),list())),
      #                          IsingSampler = list(),
      #                          huge = list(),
      #                          adalasso = list()
      #       )
      #
      #
      #     }
      #

      # # estFun:
      # if (missing(estFun)){
      #   estFun <- switch(default,
      #                    EBICglasso = qgraph::EBICglasso,
      #                    pcor = corpcor::cor2pcor,
      #                    IsingFit = IsingFit::IsingFit,
      #                    IsingSampler = IsingSampler::EstimateIsing,
      #                    huge = function(x)huge::huge.select(huge::huge(x,method = "glasso",verbose=FALSE), criterion = "ebic",verbose = FALSE),
      #                    adalasso = parcor::adalasso.net
      #   )
      # }

      # # estArgs:
      # if (missing(estArgs)){
      #   estArgs <- switch(default,
      #                     EBICglasso = list(n = nSample, returnAllResults = TRUE),
      #                     IsingFit = list(plot = FALSE, progress = FALSE),
      #                     pcor = list(),
      #                     IsingSampler = list(method = "ll"),
      #                     huge = list(),
      #                     adalasso = list()
      #   )
      # }
      #
      # # graphFun:
      # if (missing(graphFun)){
      #   graphFun <- switch(default,
      #                      EBICglasso = function(x)x[['optnet']],
      #                      IsingFit = function(x)x[['weiadj']],
      #                      pcor = function(x)as.matrix(Matrix::forceSymmetric(x)),
      #                      IsingSampler = function(x)x[['graph']],
      #                      huge = function(x)as.matrix(qgraph::wi2net(as.matrix(x$opt.icov))),
      #                      adalasso = function(x)as.matrix(Matrix::forceSymmetric(x$pcor.adalasso))
      #   )
      # }
      #
      # # graphArgs:
      # if (missing(graphArgs)){
      #   graphArgs <- switch(default,
      #                       EBICglasso = list(),
      #                       IsingFit = list(),
      #                       pcor = list(),
      #                       IsingSampler = list(),
      #                       huge = list(),
      #                       adalasso = list()
      #   )
      # }
      #
      # intFun:
      # if (missing(intFun)){
      #   intFun <- switch(default,
      #                    EBICglasso = null,
      #                    IsingFit = function(x)x[['thresholds']],
      #                    pcor = null,
      #                    IsingSampler = function(x) x[['thresholds']],
      #                    huge = null,
      #                    adalasso = null
      #   )
      # }


      # }
      #
      # if (missing(prepFun)){
      #   prepFun <- identity
      # }
      #
      # if (missing(prepArgs)){
      #   prepArgs <- list()
      # }
      #
      # if (missing(graphFun)){
      #   graphFun <- identity
      # }
      #
      # if (missing(graphArgs)){
      #   graphArgs <- list()
      # }
      #
      # if (missing(intFun)){
      #   intFun <- null
      # }
      #
      # if (missing(intArgs)){
      #   intArgs <- list()
      # }
      #
      # Function:
      # Function <- bootnet_argEstimator
      #
      #   # List of arguents:
      #   Args <- list(
      #     prepFun = prepFun,
      #     prepArgs = prepArgs,
      #     estFun = estFun,
      #     estArgs = estArgs,
      #     graphFun = graphFun,
      #     graphArgs = graphArgs,
      #     intFun = intFun,
      #     intArgs = intArgs
      #   )
      # }
      #

    }

    # Output:
    Output <- list(
      data = data,
      default = default,
      estimator = Function,
      arguments = Args
    )

    return(Output)
  }




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
    bootResults <- pbapply::pbapply(seq_len(nBoots), function(b){
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
      pb <- utils::txtProgressBar(0,nBoots,style = 3)
    }
    statTableBoots <- vector("list", nBoots)
    for (b in seq_len(nBoots)){
      statTableBoots[[b]] <- statTable(bootResults[[b]], name = paste("boot",b), alpha = alpha, computeCentrality = computeCentrality, statistics=statistics, directed=directed,  bridgeArgs=bridgeArgs, includeDiagonal=includeDiagonal)
      if (verbose){
        utils::setTxtProgressBar(pb, b)
      }
    }
    if (verbose){
      close(pb)
    }
  }  else {
    statTableBoots <- pbapply::pbapply(seq_len(nBoots),function(b){
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
    sampleTable = dplyr::ungroup(statTableOrig),
    bootTable =  dplyr::ungroup(dplyr::bind_rows(statTableBoots)),
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
