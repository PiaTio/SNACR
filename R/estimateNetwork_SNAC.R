#' This function takes data as input and produced a network. It is used inside bootnet_SNAC() and should not be used by itself.
#'
#' @param data Data matrix
#' @param default Default network estimation functions
#' @param fun Function that takes data and returns a network or list entitled "graph" and "thresholds" (optional)
#' @param labels Labels of the variables
#' @param verbose Want updates?
#' @param .dots A list
#' @param weighted Is the network weighted? Default is TRUE
#' @param signed Is the network signed? Default is TRUE
#' @param directed Is the network directed?
#' @param datatype What is the data type?
#' @param checkNumeric Internal check
#' @param ... Arguments for the 'fun' function
#' @param .input Skip most steps if supplied
#' @param memorysaver If set to FALSE (default) data, estimator and results are not stored.
#'
#' @export

estimateNetwork_SNAC <- function(
    data,
    default = c("none", "EBICglasso", "pcor","IsingFit","IsingSampler", "huge","adalasso","mgm","relimp", "cor","TMFG",
                "ggmModSelect", "LoGo","graphicalVAR", "piecewiseIsing","SVAR_lavaan",
                "GGMncv"),
    fun, # A function that takes data and returns a network or list entitled "graph" and "thresholds". optional.
    # prepFun, # Fun to produce the correlation or covariance matrix
    # prepArgs, # list with arguments for the correlation function
    # estFun, # function that results in a network
    # estArgs, # arguments sent to the graph estimation function (if missing automatically sample size is included)
    # graphFun, # set to identity if missing
    # graphArgs, # Set to null if missing
    # intFun, # Set to null if missing
    # intArgs, # Set to null if missing
    labels, # if missing taken from colnames
    verbose = TRUE, # Dummy used in cor_auto and in the future in other functions. Set to FALSE in bootnet
    # construct = c("default","function","arguments"),
    .dots = list(),
    weighted = TRUE,
    signed = TRUE,
    directed,
    datatype,
    checkNumeric = FALSE,
    # plot = TRUE, # Plot the network?
    ..., # Arguments to the 'fun' function
    .input, # Skips most of first steps if supplied
    memorysaver = FALSE # If set to FALSE data, estimator and results are not stored.
){
  construct <- "function"
  # Borsboom easter egg:
  if (default[1] == "Borsboom") return(42)

  if (default[[1]]=="glasso") default <- "EBICglasso"
  default <- match.arg(default)

  # datatype test:
  if (missing(datatype)){
    if (methods::is(data,"tsData")){
      datatype <- "graphicalVAR"
    } else {
      datatype <- "normal"
    }

  }

  if (!datatype%in% c("normal","graphicalVAR")){
    stop("Only datatypes 'normal' and 'graphicalVAR' currently supported.")
  }
  #
  # If NAs and default can't handle, stop:
  # if (any(is.na(data)) && default %in% c("huge","adalasso")){
  #   stop(paste0("Missing data not supported for default set '",default,"'. Try using na.omit(data)."))
  # }

  # First test if data is a data frame:
  if (datatype == "normal" && !(is.data.frame(data) || is.matrix(data))){
    stop("'data' argument must be a data frame")
  }

  # If matrix coerce to data frame:
  if (datatype == "normal" && is.matrix(data)){
    data <- as.data.frame(data)
  }

  if (missing(directed)){
    if (default == "graphicalVAR"){
      directed <- list(contemporaneous = FALSE, temporal = TRUE)
    } else  if (default == "SVAR_lavaan"){
      directed <- list(contemporaneous = TRUE, temporal = TRUE)
    } else if (!default %in% c("relimp","DAG")){
      directed <- FALSE
    } else {
      directed <- TRUE
    }
  }

  if (datatype == "normal"){
    N <- ncol(data)
    Np <- nrow(data)
    if (missing(labels)){
      labels <- colnames(data)
    }

    if (checkNumeric){
      # Check and remove any variable that is not ordered, integer or numeric:
      goodColumns <- sapply(data, function(x) is.numeric(x) | is.ordered(x) | is.integer(x))

      if (!all(goodColumns)){
        if (verbose){
          warning(paste0("Removing non-numeric columns: ",paste(which(!goodColumns),collapse="; ")))
        }
        data <- data[,goodColumns,drop=FALSE]
      }
    }


  } else if (datatype == "graphicalVAR"){
    N <- length(data$vars)
    Np <- nrow(data$data_c)
    if (missing(labels)){
      labels <- data$vars
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


  # Compute estimator:
  if (missing(.input)){
    .input <- checkInput(
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
      # construct=construct,
      verbose=verbose,
      .dots=.dots,
      ...
    )
  }


  # Add verbose:
  # Every estimator must have argument verbose:
  if ("verbose" %in% names(formals(.input$estimator))){
    .input$arguments$verbose <- verbose
  }

  # Unlock function:
  # Every estimator must have argument verbose:
  if ("unlock" %in% names(formals(.input$estimator))){
    .input$arguments$unlock <- TRUE
  }

  # Compute network:
  Result <- do.call(.input$estimator, c(list(data),.input$arguments))

  if (!is.list(Result)){
    sampleGraph <- Result
    intercepts <- NULL
    output <- Result
    nNode <- ncol(Result)
  } else if (is.list(Result$graph)){
    sampleGraph <- Result$graph
    intercepts <- Result$intercepts
    output <- Result$results
    nNode <- ncol(Result$graph[[1]])
  } else {
    sampleGraph <- Result$graph
    intercepts <- Result$intercepts
    output <- Result$results
    nNode <- ncol(Result$graph)
  }


  if (!is.matrix(sampleGraph)){
    if (is.list(sampleGraph)){
      if (!is.matrix(sampleGraph[[1]])){
        stop("Estimated result is not a list of matrices encoding networks.")
      }
    } else {
      stop("Estimated result is not a matrix encoding a network.")
    }
  }

  # Special data?
  if (!is.list(Result) || is.null(Result$specialData)){
    outdata <- data
    datatype <- "normal"
  } else {
    outdata <- Result$specialData$data
    datatype <- Result$specialData$type
  }


  sampleResult <- list(
    graph = sampleGraph,
    intercepts = intercepts,
    results = output,
    labels = labels,
    nNode = nNode,
    nPerson = Np,
    estimator = .input$estimator,
    arguments = .input$arguments,
    data = outdata,
    datatype = datatype,
    default = default,
    weighted = weighted,
    signed = signed,
    directed=directed,
    .input = .input,
    thresholded = FALSE,
    cv_lambda_com = Result$cv_lambda_com,
    cv_lambda_net = Result$cv_lambda_net,
    no_common = Result$no_common,
    com_struc = Result$com_struc
  )
  class(sampleResult) <- c("bootnetResult", "list")

  if (default == "graphicalVAR"){
    sampleResult$labels <- output$data$vars
  }
  if (default == "SVAR_lavaan"){
    sampleResult$labels  <- outdata$vars
  }

  # Memory save:
  if(memorysaver)
  {
    sampleResult$results <- NA
    sampleResult$estimator <- NA
    sampleResult$data <- NA
    sampleResult$.input <- NA
  }

  # Plot?
  #   if (plot){
  #     plot(sampleResult,labels=labels,layout = "spring", ...)
  #   }

  # Return network:
  return(sampleResult)
}
