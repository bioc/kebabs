##345678901234567890123456789012345678901234567890123456789012345678901234567890
#' @rdname performGridSearch
#' @aliases
#' gridSearch
#' GridSearch
#' grid.search
#' @title KeBABS Grid Search
#'
#' @description Perform grid search with one or multiple sequence kernels on
#' one or multiple SVMs with one or multiple SVM parameter sets.
#'
#' @param kernel and other parameters see \code{\link{kbsvm}}
#' @usage
#' ## kbsvm(...., kernel=list(kernel1, kernel2), pkg=pkg1, svm=svm1,
#' ##       cost=cost1, ...., cross=0, noCross=1, ....)
#'
#' ## kbsvm(...., kernel=kernel1, pkg=pkg1, svm=svm1,
#' ##       cost=c(cost1, cost2), ...., cross=0, noCross=1, ....)
#'
#' ## kbsvm(...., kernel=kernel1, pkg=c(pkg1, pkg1, pkg1),
#' ##       svm=c(svm1, svm2, svm3), cost=c(cost1, cost2, cost3), ....,
#' ##       cross=0, noCross=1, ....)
#'
#' ## kbsvm(...., kernel=kernel1, pkg=c(pkg1, pkg2, pkg3),
#' ##       svm=c(svm1, svm2, svm3), cost=c(cost1, cost2, cost3), ....,
#' ##       cross=0, noCross=1, ....)
#'
#' ## kbsvm(...., kernel=list(kernel1, kernel2, kernel3), pkg=c(pkg1, pkg2),
#' ##       svm=c(svm1, svm2), cost=c(cost1, cost2), ...., cross=0,
#' ##       noCross=1, ....)
#'
#' ## for details see below
#'
#' @details
#'
#' Overview\cr
#'
#' To simplify the selection of an appropriate sequence kernel (including
#' setting of the kernel parameters), SVM implementation and setting of SVM
#' hyperparameters KeBABS provides grid search functionality. In
#' addition to the possibility of running the same learning tasks for
#' different settings of the SVM hyperparameters the concept of grid search
#' is seen here in the broader context of finding good values for all major
#' variable parts of the learning task which includes:\cr
#' \itemize{
#' \item{selection of the sequence kernel and standard kernel parameters:
#'       spectrum, mismatch, gappy pair or motif kernel}
#' \item{selection of the kernel variant: regular, annotation-specific,
#'       position-specific or distance weighted kernel variants}
#' \item{selection of the SVM implementation via package and SVM}
#' \item{selection of the SVM hyperparameters for the SVM implementation}
#' }
#'
#' KeBABS supports the joint variation of any combination of these learning
#' aspects together with cross validation (CV) to find the best selection based
#' on cross validation performance. After the grid search the performance
#' values of the different settings and the best setting of the grid search
#' run can be retrieved from the KeBABS model with the accessor
#' \code{\link{modelSelResult}}.\cr
#'
#' Grid search is started with the method \code{\link{kbsvm}} by passing
#' multiple values to parameters for which in regular training only a single
#' parameter value is used. Multiple values can be passed for the parameter
#' \code{kernel} as list of kernel objects and for the parameters \code{pkg},
#' \code{svm} and the hyperparameters of the used SVMs as vectors (numeric or
#' integer vector dependent on the hyperparameter). The parameter cost in the
#' usage section above is just one representative of SVM hyperparameters that
#' can be varied in grid search. Following types of grid search are supported
#' (for examples see below):\cr
#'
#' \itemize{
#' \item{variation of one or multiple hyperparameter(s) for a given SVM
#'       implementation and one specific kernel by passing hyperparameter
#'       values as vectors}
#' \item{variation of the kernel parameters of a single kernel:\cr
#'       for the sequence kernels in addition to the standard kernel parameters
#'       like k for spectrum or m for gappy pair analysis can be performed in a
#'       position-independent or position-dependent manner with multiple
#'       distance weighting functions and different parameter settings for the
#'       distance weighting functions (see \code{\link{positionMetadata}}) or
#'       with or without annotation specific functionality (see
#'       \code{\link{annotationMetadata}} using one specific or multiple
#'       annotations resulting in considerable variation possibilities on the
#'       kernel side. The kernel objects for the different parameter settings
#'       of the kernel must be precreated and are passed as list to
#'       \code{\link{kbsvm}}. Usually each kernel has the best performance
#'       at differernt hyperparameter values. Therefore in general just varying
#'       the kernel parameters without varying the hyperparameter values does
#'       not make sense but both must be varied together as described below.}
#' \item{variation of multiple SVMs from the same or different R packages with
#'       identical or different SVM hyperparameters (dependent on the
#'       formulation of the SVM objective) for one specific kernel}
#' \item{combination of the previous three variants as far as runtime allows
#'       (see also runtime hints below)}
#' }
#'
#' For collecting performance values grid search is organized in a matrix
#' like manner with different kernel objects representing the rows and
#' different hyperparameter settings or SVM and hyperparameter settings as
#' columns of the matrix. If multiple hyperparameters are used on a single
#' SVM the same entry in all hyperparameter vectors is used as one parameter
#' set corresponding to a single column in the grid matrix. The same
#' applies to multiple SVMs, i.e. when multiple SVMs are used from the same
#' package the \code{pkg} parameter still must have one entry for each entry
#' in the \code{svm} parameter (see examples below). The best performing
#' setting is reported dependent on the performance objective.\cr
#'
#' Instead of a single training and test cycle for each grid point cross
#' validation should be used to get more representative results. In this case
#' CV is executed for each parameter setting. For larger datasets or kernels
#' with higher complexity the runtime for the full grid search should be
#' limited through adequate selection of the parameter \code{cross}.\cr
#'
#' Performance measures and performance objective\cr
#'
#' The usual performance measure for grid search is the cross validation error
#' which is stored by default for each grid point. For e.g. non-symmetrical
#' class distribution of the dataset other performance measures can be more
#' expressive. For such sitations also the accuracy, the balanced accuracy and
#' the Matthews correlation coefficient can be stored for a grid point (see
#' parameter \code{perfParameters} in \code{\link{kbsvm}}. (The accuracy
#' corresponds fully to the CV error because it is just the inverted measure.
#' It is included for easier comparability with the balanced accuracy). The
#' performance values can be retrieved from the model selection result object
#' with the accessor \code{\link{performance}}. The objective for selecting the
#' best performing paramters settings is by default the CV error. With the
#' parameter \code{perfObjective} in \code{\link{kbsvm}} one of the other
#' above mentioned performance parameters can be chosen as objective for the
#' best settings instead of the cross validation error. \cr
#'
#' Runtime Hints\cr
#'
#' When parameter \code{showCVTimes} in \code{\link{kbsvm}} is set to TRUE the
#' runtime for the individual cross validation runs is shown for each grid
#' point. In  this way quick runtime estimates can be gathered through running
#' the grid search for a reduced grid and extrapolating the runtimes to the
#' full grid. Display of a progress indication in grid search is available
#' with the parameter \code{showProgress} in \code{\link{kbsvm}}.\cr
#'
#' Dependent on the number of sequences, the complexity of the kernel
#' processing, the type of chosen cross validation and the degree of variation
#' of parameters in grid search the runtime can grow drastically.
#' One possible strategy for reducing the runtime could be a stepwise
#' approach searching for areas with good performance in a first coarse grid
#' search run and then refining the areas of good performance with additional
#' more fine grained grid searches.\cr
#'
#' The implementation of the sequence kernels was done with a strong focus on
#' runtime performance which brings a considerable improvement compared to
#' other implementations. In KeBABS also an interface to the very fast SVM
#' implementations in package LiblineaR is available. Beyond these performance
#' improvements KeBABS also supports the generation of sparse explicit
#' representations for every sequence kernel which can be used instead of the
#' kernel matrix for learning. In many cases especially with a large number of
#' samples where the kernel matrix would become too large this alternative
#' provides additional dynamical benefits. The current implementation of grid
#' search does not make use of multi-core infrastructures, the entire
#' processing is done on a single core.\cr
#'
#' @return grid search stores the results in the KeBABS model. They can be 
#' retrieved with the accessor \code{\link{modelSelResult}{KBModel}}.
#'
#' @seealso \code{\link{kbsvm}},
#' \code{\link{spectrumKernel}}, \code{\link{mismatchKernel}},
#' \code{\link{gappyPairKernel}}, \code{\link{motifKernel}},
#' \code{\link{positionMetadata}}, \code{\link{annotationMetadata}},
#' \code{\link{performModelSelection}}
#'
#' @examples
#' ## load transcription factor binding site data
#'
#' data(TFBS)
#' enhancerFB
#' ## The C-svc implementation from LiblineaR is chosen for most of the
#' ## examples because it is the fastest SVM implementation. With SVMs from
#' ## other packages slightly better results could be achievable.
#' ## To get a realistic image of possible performance values, kernel behavior
#' ## and speed of grid search together with 10-fold cross validation a
#' ## resonable number of sequences is needed which would exceed the runtime
#' ## restrictions for automatically executed examples. Therefore the grid
#' ## search examples must be run manually. In these examples we use the full
#' ## dataset for grid search.
#' train <- sample(1:length(enhancerFB), length(enhancerFB))
#'
#' ## grid search with single kernel object and multiple hyperparameter values
#' ## create gappy pair kernel with normalization
#' gappyK1M3 <- gappyPairKernel(k=1, m=3)
#' ## show details of single gappy pair kernel object
#' gappyK1M3
#'
#' ## grid search for a single kernel object and multiple values for cost
#' pkg <- "LiblineaR"
#' svm <- "C-svc"
#' cost <- c(0.01,0.1,1,10,100,1000)
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=gappyK1M3,
#'                pkg=pkg, svm=svm, cost=cost, explicit="yes", cross=3)
#'
#' ## show grid search results
#' modelSelResult(model)
#'
#' \dontrun{
#' ## create the list of spectrum kernel objects with normalization and
#' ## kernel parameters values for k from 1 to 5
#' specK15 <- spectrumKernel(k=1:5)
#' ## show details of the four spectrum kernel objects
#' specK15
#'
#' ## run grid search with several kernel parameter settings for the
#' ## spectrum kernel with a single SVM parameter setting
#' ## ATTENTION: DO NOT USE THIS VARIANT!
#' ## This variant does not bring comparable performance for the different
#' ## kernel parameter settings because usually the best performing
#' ## hyperparameter values could be quite different for different kernel
#' ## parameter settings or between different kernels, grid search for
#' ## multiple kernel objects should be done as shown in the next example
#' pkg <- "LiblineaR"
#' svm <- "C-svc"
#' cost <- 2
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=specK15,
#'                pkg=pkg, svm=svm, cost=cost, explicit="yes", cross=10)
#'
#' ## show grid search results
#' modelSelResult(model)
#'
#' ## grid search with multiple kernel objects and multiple values for
#' ## hyperparameter cost
#' pkg <- "LiblineaR"
#' svm <- "C-svc"
#' cost <- c(0.01,0.1,1,10,50,100,150,200,500,1000)
#' model <- kbsvm(x=enhancerFB, sel=train, y=yFB[train], kernel=specK15,
#'                pkg=pkg, svm=svm, cost=cost, explicit="yes", cross=10,
#'                showProgress=TRUE)
#'
#' ## show grid search results
#' modelSelResult(model)
#'
#' ## grid search for a single kernel object with multiple SVMs
#' ## from different packages
#' ## here with display of cross validation runtimes for each grid point
#' ## pkg, svm and cost vectors must have same length and the corresponding
#' ## entry in each of these vectors are one SVM + SVM hyperparameter setting
#' pkg <- rep(c("kernlab", "e1071", "LiblineaR"),3)
#' svm <- rep("C-svc", 9)
#' cost <- rep(c(0.01,0.1,1),each=3)
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=gappyK1M3,
#'                pkg=pkg, svm=svm, cost=cost, explicit="yes", cross=3,
#'                showCVTimes=TRUE)
#'
#' ## show grid search results
#' modelSelResult(model)
#'
#' ## run grid search for a single kernel with multiple SVMs from same package
#' ## here all from LiblineaR: C-SVM, L2 regularized SVM with L2 loss and
#' ## SVM with L1 regularization and L2 loss
#' ## attention: for different formulation of the SMV objective use different
#' ## values for the hyperparameters even if they have the same name
#' pkg <- rep("LiblineaR", 9)
#' svm <- rep(c("C-svc","l2rl2l-svc","l1rl2l-svc"), each=3)
#' cost <- c(1,150,1000,1,40,100,1,40,100)
#' model <- kbsvm(x=enhancerFB, sel=train, y=yFB[train], kernel=gappyK1M3,
#'                pkg=pkg, svm=svm, cost=cost, explicit="yes", cross=3)
#'
#' ## show grid search results
#' modelSelResult(model)
#'
#' ## create the list of kernel objects for gappy pair kernel
#' gappyK1M15 <- gappyPairKernel(k=1, m=1:5)
#' ## show details of kernel objects
#' gappyK1M15
#'
#' ## run grid search with progress indication with ten kernels and ten
#' ## hyperparameter values for cost and 10 fold cross validation on full
#' ## dataset (500 samples)
#' pkg <- rep("LiblineaR", 10)
#' svm <- rep("C-svc", 10)
#' cost <- c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000,100000)
#' model <- kbsvm(x=enhancerFB, y=yFB, kernel=c(specK15, gappyK1M15),
#'                pkg=pkg, svm=svm, cost=cost, cross=10, explicit="yes",
#'                showCVTimes=TRUE, showProgress=TRUE)
#'
#' ## show grid search results
#' modelSelResult(model)
#' }
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.
#' @keywords kbsvm
#' @keywords grid search
#' @keywords model selection
#' @keywords methods
#'

performGridSearch <- function(object, model, y, explicit, featureWeights,
                              weightLimit, sel, features=NULL, groupBy,
                              perfParameters, perfObjective, addArgs,
                              includeFullModel=FALSE, showProgress,
                              showCVTimes, verbose)
{
    if (is(object, "kernelMatrix"))
        object <- as.KernelMatrix(object)

    ## allow overwritting of explicit for multiple kernels / multiple SVMS
    overwriteExplicit <- FALSE
    singleKernel <- NULL
    explicitType <- model@svmInfo@reqExplicitType

    ## preload SparseM to avoid loading message to interfere with
    ## progress messages
    if (explicitType %in% c("auto","sparse"))
        library(SparseM, warn.conflicts=FALSE)

    if (explicit == "auto")
    {
        if (is(object, "KernelMatrix"))
            explicit <- "no"
        else if (is(object, "ExplicitRepresentation"))
        {
            explicit <- "yes"
            
            if (is(object, "ExplicitRepresentationDense"))
                explicitType <- "dense"
            else
                explicitType <- "sparse"
        }
        else ## XStringSet or BioVector
        {
            if (all(unlist(lapply(model@modelSelResult@gridRows,
                                  supportsExplicitRep))))
            {
                explicit <- "yes"

                if (explicitType == "auto")
                {
                    if (model@ctlInfo@onlyDense == TRUE)
                        explicitType <- "dense"
                    else
                        explicitType <- "sparse"
                }
            }
            else if (length(model@modelSelResult@gridRows) < 2)
            {
                ## one kernel which does not support ER
                explicit <- "no"
            }
            else
            {
                ## leave value auto for multiple kernels
                overwriteExplicit <- TRUE
            }
        }
    }

    if (explicit == "auto")
        overwriteExplicit <- TRUE

    if (explicit == "no" &&
        any(mapply(svmSupportsOnlyExplicitRep,
                   model@svmInfo@reqPackage, model@svmInfo@reqSVM)))
        stop("SVM package does not support a kernel matrix\n")

    if (explicitType == "auto")
    {
        if (is(object, "ExplicitRepresentationDense") ||
            model@ctlInfo@onlyDense)
            explicitType <- "dense"
        else
            explicitType <- "sparse"
    }

    if (featureWeights == "auto")
        featureWeights <- "no"

    if (is(object, "KernelMatrix"))
    {
        if (explicit != "no")
        {
            stop("grid search with explicit representation is\n",
                 "        not possible for call with kernel matrix\n")
        }

        if (length(model@svmInfo@reqKernel) > 0)
        {
            stop("grid search with kernels not possible for call\n",
                 "         with kernel matrix\n")
        }

        noRows <- 1
        rowName <- "No Kernel"
        noKernel <- TRUE

        ## subset kernel matrix immediately
        if (length(sel) > 0)
        {
            object <- subsetSeqRep(x=object, sel=sel)
            sel <- integer(0)
        }

        callWithOrigData <- TRUE
    }
    else
    {
        if (inherits(object, "ExplicitRepresentation"))
            callWithOrigData <- TRUE
        else
            callWithOrigData <- FALSE

        if (is.null(model@svmInfo@reqKernel))
        {
            noRows <- 1
            rowName <- "No Kernel"
            noKernel <- TRUE
            
            if (!is(object, "ExplicitRepresentation"))
                stop("kernel is missing for grid search on sequences\n")
        }
        else
        {
            noRows <- length(model@modelSelResult@gridRows)
            rowName <- "Kernel"
            noKernel <- FALSE

            if (noRows == 1)
                singleKernel <- model@modelSelResult@gridRows[[1]]
            else
            {
                if (length(features) > 0)
                {
                    stop("feature subsetting not supported for multiple",
                         " kernels\n")
                }
            }
        }

    }


    if (length(sel) > 0)
        numSamples <- length(sel)
    else
        numSamples <- getNoOfElementsOfSeqRep(object)

    if (model@modelSelResult@cross >= numSamples)
        stop("'cross' must be smaller than the number of samples\n")

    if (!is.null(model@modelSelResult@gridCols))
        noCols <- nrow(model@modelSelResult@gridCols)
    else
        noCols <- 1

    cvErrors <- matrix(0, nrow=noRows, ncol=noCols)
    colnames(cvErrors) <- paste("GridCol", 1: ncol(cvErrors), sep="_")
    sumAlphas <- matrix(0, nrow=noRows, ncol=noCols)
    colnames(sumAlphas) <- paste("GridCol", 1: ncol(sumAlphas), sep="_")
    noSV <- matrix(0, nrow=noRows, ncol=noCols)
    colnames(noSV) <- paste("GridCol", 1: ncol(noSV), sep="_")
    cvTimes <- matrix(0, nrow=noRows, ncol=noCols)
    colnames(cvTimes) <- paste("GridCol", 1: ncol(cvTimes), sep="_")

    collectACC <- "ACC" %in% model@modelSelResult@perfParameters
    collectBACC <- "BACC" %in% model@modelSelResult@perfParameters
    collectMCC <- "MCC" %in% model@modelSelResult@perfParameters
    collectAUC <- "AUC" %in% model@modelSelResult@perfParameters

    if (collectACC)
    {
        gridACC <- matrix(NA, nrow=noRows, ncol=noCols)
        colnames(gridACC) <- paste("GridCol", 1: ncol(cvErrors), sep="_")
    }

    if (collectBACC)
    {
        gridBACC <- matrix(NA, nrow=noRows, ncol=noCols)
        colnames(gridBACC) <- paste("GridCol", 1: ncol(cvErrors), sep="_")
    }

    if (collectMCC)
    {
        gridMCC <- matrix(NA, nrow=noRows, ncol=noCols)
        colnames(gridMCC) <- paste("GridCol", 1: ncol(cvErrors), sep="_")
    }

    if (collectAUC)
    {
        gridAUC <- matrix(NA, nrow=noRows, ncol=noCols)
        colnames(gridAUC) <- paste("GridCol", 1: ncol(cvErrors), sep="_")
    }

    if (noKernel)
    {
        rownames(cvErrors) <- rowName
        rownames(sumAlphas) <- rowName
        rownames(noSV) <- rowName
        rownames(cvTimes) <- rowName

        if (collectACC)
            rownames(gridACC) <- rowName

        if (collectBACC)
            rownames(gridBACC) <- rowName

        if (collectMCC)
            rownames(gridMCC) <- rowName
    }
    else
    {
        rownames(cvErrors) <- paste(rowName, 1: nrow(cvErrors), sep="_")
        rownames(sumAlphas) <- paste(rowName, 1: nrow(sumAlphas), sep="_")
        rownames(noSV) <- paste(rowName, 1: nrow(noSV), sep="_")
        rownames(cvTimes) <- paste(rowName, 1: nrow(cvTimes), sep="_")

        if (collectACC)
            rownames(gridACC) <- paste(rowName, 1: nrow(gridACC), sep="_")

        if (collectBACC)
            rownames(gridBACC) <- paste(rowName, 1: nrow(gridBACC), sep="_")

        if (collectMCC)
            rownames(gridMCC) <- paste(rowName, 1: nrow(gridMCC), sep="_")

        if (collectAUC)
            rownames(gridAUC) <- paste(rowName, 1: nrow(gridAUC), sep="_")
    }

    if (perfObjective =="ACC")
        bestResult <- Inf  ## because error rate is used instead
    else
        bestResult <- - Inf

    bestRow <- NA
    bestCol <- NA
    precompObj <- NULL

    model@svmInfo@selPackage <- model@svmInfo@reqPackage
    model@svmInfo@selSVM <- model@svmInfo@reqSVM
    model@ctlInfo@classification <-
        isClassification(model@svmInfo@selPackage,model@svmInfo@selSVM)

    if (model@ctlInfo@classification == TRUE && model@numClasses == 0)
    {
        ## check for multiclass to get the no of classes reliably
        model@numClasses <- length(unique(y))

        if (model@numClasses > 2)
            model@ctlInfo@multiclassType <- getMulticlassType(model)
    }

    if (explicit == "yes")
    {
        model@svmInfo@selExplicit <- TRUE
        model@ctlInfo@selMethod <- "explicitRep"

        ## precompute explicit representation
        if (!is.null(singleKernel) &&
            (is(object, "XStringSet") || is(object, "BioVector")))
        {
            model@ctlInfo@sparse <- (explicitType == "sparse")
            
            obj <- getExRep(x=object, kernel=singleKernel,
                            sparse=model@ctlInfo@sparse, selx=sel,
                            features=features)
            
            model@trainingFeatures <- colnames(obj)
            sel <- integer(0)
        }
    }
    else if (explicit == "no")
    {
        model@svmInfo@selExplicit <- FALSE
        model@ctlInfo@selMethod <- "KernelMatrix"

        ## precompute kernel matrix
        if (!is.null(singleKernel) && !is(object, "KernelMatrix"))
        {
            if (inherits(object, "ExplicitRepresentation"))
                obj <- linearKernel(x=object, selx=sel)
            else
                obj <- singleKernel(x=object, selx=sel)

            sel <- integer(0)
        }
    }

    selOrig <- sel
    yOrig <- y
    groupByOrig <- groupBy

    if (explicit == "yes")
        model@svmInfo@explicitKernel <- model@svmInfo@reqFeatureType

    if (includeFullModel)
        fullModel <- model

    kernel <- NULL

    if (is.null(model@modelSelResult@gridCols))
    {
        ## convert parameters to target SVM
        model@svmInfo <- convertSVMParameters(model)
        convertSVMParameters <- FALSE
        yRange <- 1
        noSVMPar <- 0
    }
    else
    {
        convertSVMParameters <- TRUE
        yRange <- nrow(model@modelSelResult@gridCols)
        noSVMPar <- ncol(model@modelSelResult@gridCols)
        startIndex <- 0

        if (colnames(model@modelSelResult@gridCols)[1] == "pkg")
        {
            noSVMPar <- noSVMPar - 2
            startIndex <- 2
        }
    }

    for (i in 1:noRows)
    {
        sel <- selOrig
        y <- yOrig
        groupBy <- groupByOrig

        if (noKernel)
        {
            currRowName <- rowName
        }
        else
        {
            currRowName <- names(model@modelSelResult@gridRows)[i]
            kernel <- model@modelSelResult@gridRows[[i]]
        }

        if (showProgress)
        {
            if (i == 1)
                message("\nGrid Search Progress:\n", appendLF=FALSE)

            if (noKernel)
                message("\nNo Kernel  ", appendLF=FALSE)
            else
                message("\n", currRowName, "  ", appendLF=FALSE)
        }

        if (!noKernel && is.null(singleKernel))
        {
            if (overwriteExplicit)
            {
                ## only XStringSet or BioVector with mult. kernels
                if (supportsExplicitRep(kernel))
                {
                    explicit <- "yes"
                    model@svmInfo@selExplicit <- TRUE
                    model@ctlInfo@selMethod <- "explicitRep"
                    model@svmInfo@explicitKernel <- model@svmInfo@reqFeatureType
                }
                else
                {
                    explicit <- "no"
                    model@svmInfo@selExplicit <- FALSE
                    model@ctlInfo@selMethod <- "KernelMatrix"
                }
            }
            
            if (is(object, "BioVector") || is(object, "XStringSet"))
            {
                if (model@ctlInfo@selMethod == "KernelMatrix")
                    obj <- kernel(x=object, selx=sel)
                else
                {
                    obj <- getExRep(x=object, kernel=kernel,
                                    sparse=(explicitType=="sparse"), selx=sel)

                    model@trainingFeatures <- colnames(obj)
                }

                sel <- NULL
            }
            else if (inherits(object, "ExplicitRepresentation"))
            {
                if (model@ctlInfo@selMethod == "KernelMatrix")
                    obj <- linearKernel(x=object, selx=sel)
                else
                {
                    if (length(sel) > 0)
                        obj <- object[sel,]
                    else
                        obj <- object
                }

                sel <- integer(0)
            }
        }

        colCount <- 0

        for (j in 1:yRange)
        {
            tempModel <- model
            tempModel@svmInfo@selKernel <- kernel
            tempModel@ctlInfo@sparse <- (explicitType == "sparse")

            if (!is.null(model@modelSelResult@gridCols) &&
                startIndex == 2)
            {
                tempModel@svmInfo@selPackage <-
                    as.character(model@modelSelResult@gridCols[j,1])
                tempModel@svmInfo@selSVM <-
                    as.character(model@modelSelResult@gridCols[j,2])

                ## for quadratic kernel CV is always called with linear ex rep
                ## and internally retrieves the quadratic ex rep if necessary
            }

            tempModel@ctlInfo@classification <- isClassification(
                                                tempModel@svmInfo@selPackage,
                                                tempModel@svmInfo@selSVM)

            if (noSVMPar > 0)
            {
                for (j1 in 1:noSVMPar)
                {
                    tempModel@svmInfo <- addOrReplaceSVMParameters(
                        colnames(model@modelSelResult@gridCols)[startIndex+j1],
                        model@modelSelResult@gridCols[j,startIndex+j1],
                        tempModel)
                }

                if (convertSVMParameters)
                    tempModel@svmInfo <- convertSVMParameters(tempModel)
            }

            if (callWithOrigData)
            {
                if (inherits(object, "ExplicitRepresentationSparse") &&
                    tempModel@svmInfo@selPackage == "kernlab")
                {
                    cvTimes[i,j] <- system.time(tempModel <-
                            performCrossValidation(
                                object=as(object,"ExplicitRepresentationDense"),
                                x=NULL, y=y, sel=integer(0), model=tempModel,
                                cross=tempModel@modelSelResult@cross,
                                noCross=tempModel@modelSelResult@noCross,
                                groupBy=groupBy, perfParameters=
                                model@modelSelResult@perfParameters,
                                verbose=verbose))[3]   ## elapsed time
                }
                else
                {
                    cvTimes[i,j] <- system.time(tempModel <-
                            performCrossValidation(object=object,
                                x=NULL, y=y, sel=integer(0), model=tempModel,
                                cross=tempModel@modelSelResult@cross,
                                noCross=tempModel@modelSelResult@noCross,
                                groupBy=groupBy, perfParameters=
                                model@modelSelResult@perfParameters,
                                verbose=verbose))[3]   ## elapsed time
                }
            }
            else
            {
                if (inherits(obj, "ExplicitRepresentationSparse") &&
                    tempModel@svmInfo@selPackage == "kernlab")
                {
                    cvTimes[i,j] <- system.time(tempModel <-
                            performCrossValidation(
                                object=as(obj, "ExplicitRepresentationDense"),
                                x=object, y=y, sel=integer(0), model=tempModel,
                                cross=tempModel@modelSelResult@cross,
                                noCross=tempModel@modelSelResult@noCross,
                                groupBy=groupBy, perfParameters=
                                model@modelSelResult@perfParameters,
                                verbose=verbose))[3]   ## elapsed time
                }
                else
                {
                    cvTimes[i,j] <- system.time(tempModel <-
                            performCrossValidation(
                                object=obj, x=object, y=y,
                                sel=integer(0), model=tempModel,
                                cross=tempModel@modelSelResult@cross,
                                noCross=tempModel@modelSelResult@noCross,
                                groupBy=groupBy, perfParameters=
                                model@modelSelResult@perfParameters,
                                verbose=verbose))[3]   ## elapsed time
                }
            }

            ## get CV result
            cvErrors[i, j] <- mean(tempModel@cvResult@cvError)
            noSV[i,j] <- mean(tempModel@cvResult@noSV)
            sumAlphas[i,j] <- mean(tempModel@cvResult@sumAlphas)

            if (collectACC)
                gridACC[i,j] <- mean(tempModel@cvResult@ACC)

            if (collectBACC)
                gridBACC[i,j] <- mean(tempModel@cvResult@BACC)

            if (collectMCC)
                gridMCC[i,j] <- mean(tempModel@cvResult@MCC)

            if (collectAUC)
                gridAUC[i,j] <- mean(tempModel@cvResult@AUC)

            if (perfObjective == "ACC")
            {
                if (!is.na(cvErrors[i,j]) && (cvErrors[i, j] < bestResult))
                {
                    bestResult <- cvErrors[i, j]
                    bestRow <- i
                    bestCol <- j
                }
            }
            else if (perfObjective == "BACC")
            {
                if (!is.na(gridBACC[i,j]) && (gridBACC[i, j] > bestResult))
                {
                    bestResult <- gridBACC[i, j]
                    bestRow <- i
                    bestCol <- j
                }
            }
            else if (perfObjective == "MCC")
            {
                if (!is.na(gridMCC[i,j]) && (gridMCC[i, j] > bestResult))
                {
                    bestResult <- gridMCC[i, j]
                    bestRow <- i
                    bestCol <- j
                }
            }
            else if (perfObjective == "AUC")
            {
                if (!is.na(gridAUC[i,j]) && (gridAUC[i, j] > bestResult))
                {
                    bestResult <- gridAUC[i, j]
                    bestRow <- i
                    bestCol <- j
                }
            }

            if (showProgress)
            {
                colCount <- colCount + 1

                if (colCount > 60)
                {
                    colCount <- 0
                    message("\n            ", appendLF=FALSE)
                }
                message(".", appendLF=FALSE)
            }
        }
    }

    if (showProgress)
        message("\n", appendLF=FALSE)

    if (showCVTimes)
    {
        cat("\nGrid Search Times:\n\n")
        print(cvTimes)
        cat("\n")
    }

    model@modelSelResult@gridErrors <- cvErrors
    model@modelSelResult@gridNoSV <- noSV
    model@modelSelResult@gridSumAlphas <- sumAlphas
    model@modelSelResult@smallestCVError <- bestResult

    if (collectACC)
        model@modelSelResult@gridACC <- gridACC

    if (collectBACC)
        model@modelSelResult@gridBACC <- gridBACC

    if (collectMCC)
        model@modelSelResult@gridMCC <- gridMCC

    if (collectAUC)
        model@modelSelResult@gridAUC <- gridAUC

    if (!is.null(model@svmInfo@reqKernel))
    {
        model@modelSelResult@selGridRow <-
                                model@modelSelResult@gridRows[[bestRow]]
    }
    else
        model@modelSelResult@selGridRow <- NULL

    if (!is.null(model@modelSelResult@gridCols))
    {
        model@modelSelResult@selGridCol <-
                    model@modelSelResult@gridCols[bestCol, ,
                                                  drop=FALSE]

        colnames(model@modelSelResult@selGridCol) <-
                    colnames(model@modelSelResult@gridCols)
    }

    if (includeFullModel)
    {
        ## set up info for selected row and col
        if (!noKernel)
        {
            selKernel <- model@modelSelResult@selGridRow
            fullModel@svmInfo@selKernel <- selKernel
        }

        ## set up SVM parameters
        if (!is.null(model@modelSelResult@gridCols))
        {
            if (colnames(model@modelSelResult@selGridCol)[1] == "pkg")
            {
                fullModel@svmInfo@selPackage <-
                        as.character(model@modelSelResult@selGridCol[1,1])
                fullModel@svmInfo@selSVM <-
                        as.character(model@modelSelResult@selGridCol[1,2])
            }
        }

        if (noSVMPar > 0)
        {
            for (j1 in 1:noSVMPar)
            {
                fullModel@svmInfo <- addOrReplaceSVMParameters(
                    colnames(model@modelSelResult@selGridCol)[startIndex+j1],
                    model@modelSelResult@selGridCol[startIndex+j1],
                    fullModel)
            }

            if (convertSVMParameters)
                fullModel@svmInfo <- convertSVMParameters(fullModel)
        }

        ## set up explicit and explicit type parameter
        if (overwriteExplicit == TRUE)
        {
            if (supportsExplicitRep(selKernel))
            {
                explicit <- "yes"
                fullModel@svmInfo@selExplicit <- TRUE
                fullModel@ctlInfo@selMethod <- "explicitRep"
                fullModel@svmInfo@explicitKernel <- model@svmInfo@reqFeatureType

                if (explicitType == "auto")
                {
                    if (fullModel@svmInfo@selPackage == "kernlab")
                        explicitType <- "dense"
                    else
                        explicitType <- "sparse"
                }
            }
            else
            {
                explicit <- "no"
                fullModel@svmInfo@selExplicit <- FALSE
                fullModel@ctlInfo@selMethod <- "KernelMatrix"
            }
        }

        if (fullModel@svmInfo@selPackage == "kernlab")
        {
            fullModel@ctlInfo@sparse <- FALSE
            
            if (callWithOrigData)
            {
                if (is(object, "ExplicitRepresentationSparse"))
                    object <- as(object, "ExplicitRepresentationDense")
            }
            else
            {
                if (is(obj, "ExplicitRepresentationSparse"))
                    obj <- as(obj, "ExplicitRepresentationDense")
            }
        }
        else
        {
            if (callWithOrigData && is(object, "ExplicitRepresentationDense"))
                fullModel@ctlInfo@sparse <- FALSE
            else
                fullModel@ctlInfo@sparse <- TRUE
        }
        
        if (noRows > 1 && (is(object, "BioVector") ||
                           is(object, "XStringSet")))
        {
            if (explicit == "yes")
            {
                obj <- getExRep(x=object, kernel=selKernel,
                                sparse=fullModel@ctlInfo@sparse,
                                features=features)

                fullModel@trainingFeatures <- colnames(obj)
            }
            else
                obj <- selKernel(object)
        }

        ## invoke selected svm
        allArgs <- fullModel@svmInfo@selSVMPar

        if (callWithOrigData)
            allArgs[["x"]] <- object
        else
            allArgs[["x"]] <- obj

        allArgs[["y"]] <- y
        allArgs[["svmInfo"]] <- fullModel@svmInfo
        allArgs[["verbose"]] <- verbose
        fullModel@svmModel <- do.call(trainSVM, allArgs)

        if (is.null(fullModel@svmModel))
        {
            stop(paste(fullModel@svmInfo@selPackage, "-",
                       fullModel@svmInfo@selSVM, "did not return a model\n"))
        }

        if (fullModel@ctlInfo@classification == TRUE)
        {
            fullModel@numClasses <- getSVMSlotValue("numClasses", fullModel)
            fullModel@classNames <- getSVMSlotValue("classNames", fullModel)
        }

        ind <- getSVMSlotValue("svIndex", fullModel)

        if (length(ind) > 0)
        {
            if (callWithOrigData)
            {
                fullModel@SV <- object[ind]
                fullModel@svIndex <- ind
                fullModel@alphaIndex <- getSVMSlotValue("alphaIndex", fullModel)

            }
            else
            {
                if (is(obj, "KernelMatrix"))
                    fullModel@SV <- object[ind]
                else
                    fullModel@SV <- obj[ind]

                fullModel@svIndex <- ind
                fullModel@alphaIndex <- getSVMSlotValue("alphaIndex", fullModel)
            }
        }

        if (fullModel@svmInfo@featureWeights != "no")
        {
            if (length(ind) > 0 && inherits(obj, "ExplicitRepresentation"))
                exRepSV <- obj[ind]
            else
                exRepSV <- NULL

            fullModel@featureWeights <-
                getFeatureWeights(model=fullModel,
                                  exrep=exRepSV,
                                  weightLimit=fullModel@svmInfo@weightLimit)

            fullModel@b <- getSVMSlotValue("b", fullModel)
            fullModel@SV <- NULL

            if (length(fullModel@featureWeights) < 1)
                stop("no feature weights available\n")
        }

        model@modelSelResult@fullModel <- fullModel
    }

    return(model)
}
