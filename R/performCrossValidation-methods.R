##345678901234567890123456789012345678901234567890123456789012345678901234567890
performCrossValidation.KernelMatrix <- function(object, x, y, sel, model,
                            cross, noCross, groupBy, perfParameters, verbose)
{
    ## k-fold cross validation
    ## for cross = -1 run LooCV, for cross > 0 k-fold cv

    if (model@svmInfo@featureWeights != "no" && length(x) > 0 &&
        supportsExplicitRep(model@svmInfo@selKernel))
    {
        exRep <- getExRep(x=x, selx=sel, sparse=model@ctlInfo@sparse,
                          kernel=model@svmInfo@selKernel)
    }
    else
        exRep <- NULL

    ## subset immediately because of cv
    if (length(sel) > 0)
    {
        object <- object[sel,sel]
        selOrig <- sel
        sel <- NULL
    }
    else
        selOrig <- NULL

    if (cross == -1)
        numFolds <- nrow(object)
    else
    {
        if (cross >= nrow(object))
            stop("'cross' must be smaller than the number of samples\n")

        numFolds <- cross
    }

    if (!is.null(groupBy))
    {
        ## groupBy is a numeric vector or a factor
        numGroups <- length(table(as.numeric(groupBy)))
    }
    else
        numGroups <- 0

    model@cvResult <- new("CrossValidationResult")
    model@cvResult@cross <- cross
    model@cvResult@noCross <- noCross
    model@cvResult@groupBy <- groupBy
    model@cvResult@perfParameters <- perfParameters
    model@cvResult@folds <- as.list(rep(NA_real_, noCross * numFolds))
    model@cvResult@foldErrors <- matrix(Inf, nrow=noCross, ncol=numFolds)

    collectACC <- "ACC" %in% perfParameters
    collectBACC <- "BACC" %in% perfParameters
    collectMCC <- "MCC" %in% perfParameters
    collectAUC <- "AUC" %in% perfParameters

    if (collectACC)
        model@cvResult@foldACC <- matrix(Inf, nrow=noCross, ncol=numFolds)

    if (collectBACC)
        model@cvResult@foldBACC <- matrix(Inf, nrow=noCross, ncol=numFolds)

    if (collectMCC)
        model@cvResult@foldMCC <- matrix(Inf, nrow=noCross, ncol=numFolds)

    if (collectAUC)
        model@cvResult@foldAUC <- matrix(Inf, nrow=noCross, ncol=numFolds)

    if (model@ctlInfo@classification == TRUE && model@numClasses == 0)
    {
        ## check for multiclass to get the no of classes reliably
        model@numClasses <- length(unique(y))

        if (model@numClasses > 2)
            model@ctlInfo@multiclassType <- getMulticlassType(model)
    }


    for (i in 1:noCross)
    {
        ## split training data into folds

        randomSamples <- sample(1:nrow(object), nrow(object))

        if (numGroups > 0)
        {
            ## assign full groups to folds
            ## for cross equal to numGroups => Leave Group Out CV
            groupsInFolds <-
                    suppressWarnings(split(sample(1:numGroups), 1:cross))
            grouping <- as.numeric(groupBy)[randomSamples]
            folds <- lapply(1:cross,
                        function(i){randomSamples[which(grouping %in%
                                                    groupsInFolds[[i]])]})
        }
        else
        {
            ## k-fold or Leave One Out CV
            folds <- suppressWarnings(split(randomSamples, 1:numFolds))
        }

        numFolds <- length(folds)

        if (numFolds < 2)
            stop("less than two folds generated\n")

        model@cvResult@folds[((i-1)*numFolds+1):((i)*numFolds)] <- folds
        noOfSVs <- rep(NA_real_, numFolds)
        sumAlphas <- rep(NA_real_, numFolds)
        foldError <- rep(NA_real_, numFolds)
        predDecValues <- NA

        if (collectACC)
            foldACC <- rep(NA_real_, numFolds)

        if (collectBACC)
            foldBACC <- rep(NA_real_, numFolds)

        if (collectMCC)
            foldMCC <- rep(NA_real_, numFolds)

        if (collectAUC)
            foldAUC <- rep(NA_real_, numFolds)

        for (j in 1:numFolds)
        {
            trainSuccess <- TRUE
            tempModel <- model
            trainIndices <- unlist(folds[-j], use.names=FALSE)

            allArgs <- tempModel@svmInfo@selSVMPar
            allArgs[["x"]] <- object[trainIndices, trainIndices]
            allArgs[["y"]] <- y[trainIndices]
            allArgs[["svmInfo"]] <- tempModel@svmInfo
            allArgs[["verbose"]] <- verbose

            tempModel@svmModel <- tryCatch({do.call(trainSVM, allArgs)},
                warning =
                  function(w)
                  {
                    message(paste("\nwarning during training for CV fold ",
                                  j,":"))
                    message(w)
                  },
                error =
                  function(e)
                  {
                    message(paste("\nerror during training for CV fold ",
                                  j,":"))
                    message(e)
                    trainSuccess <- FALSE
                  })

            if (!is.logical(tempModel@svmModel))
            {
                if (is.null(tempModel@svmModel))
                {
                    stop(paste(tempModel@svmInfo@selPackage, "-",
                        tempModel@svmInfo@selSVM, "did not return a model\n"))
                }

                if (tempModel@ctlInfo@classification == TRUE)
                {
                    tempModel@numClasses <-
                        getSVMSlotValue("numClasses", tempModel)
                    tempModel@classNames <-
                        getSVMSlotValue("classNames", tempModel)
                }

                svIndices <- getSVMSlotValue("svIndex", tempModel)

                if (length(svIndices) > 0)
                {
                    alphaIndices <- getSVMSlotValue("alphaIndex", tempModel)

                    if (length(selOrig) > 0)
                    {
                        tempModel@SV <- x[selOrig[trainIndices[svIndices]]]
                        tempModel@svIndex <- selOrig[trainIndices[svIndices]]

                        if (is.list(alphaIndices))
                        {
                            tempModel@alphaIndex <- lapply(alphaIndices,
                                function(x) selOrig[trainIndices[x]])
                        }
                        else
                        {
                            tempModel@alphaIndex <-
                                selOrig[trainIndices[alphaIndices]]
                        }
                    }
                    else
                    {
                        tempModel@SV <- x[trainIndices[svIndices]]
                        tempModel@svIndex <- trainIndices[svIndices]

                        if (is.list(alphaIndices))
                        {
                            tempModel@alphaIndex <- lapply(alphaIndices,
                                function(x) trainIndices[x])
                        }
                        else
                            tempModel@alphaIndex <- trainIndices[alphaIndices]
                    }
                }
                else
                {
                    print("no support vectors found")
                    ## no kernel matrix in LiblineaR
                    next
                }

                if (tempModel@svmInfo@featureWeights != "no" && !is.null(exRep)
                    && tempModel@numClasses == 2)
                {
                    tempModel@featureWeights <-
                         getFeatureWeights(model=tempModel,
                                    exrep=exRep[trainIndices[svIndices],],
                                    weightLimit=tempModel@svmInfo@weightLimit)

                    tempModel@b <- getSVMSlotValue("b", tempModel)
                    tempModel@SV <- NULL

                    pred <- predict(object=tempModel,
                                    x=exRep[folds[[j]],],
                                    predictionType="response", verbose=verbose)
                    
                    if (collectAUC)
                    {
                        predDecValues <- predict(object=tempModel,
                                                 x=exRep[folds[[j]],],
                                                 predictionType="decision",
                                                 verbose=verbose)
                    }
                }
                else
                {
                    if (tempModel@svmInfo@selPackage == "kernlab")
                        currInd <- trainIndices[svIndices]
                    else
                        currInd <- trainIndices
                    
                    ## predict with rectangular kernel matrix
                    pred <- predict(object=tempModel,
                                    x=object[folds[[j]], currInd],
                                    predictionType="response", verbose=verbose)
                    
                    if (collectAUC)
                    {
                        predDecValues <- predict(object=tempModel,
                                                 x=object[folds[[j]],
                                                          currInd],
                                                 predictionType="decision",
                                                 verbose=verbose)
                    }
                }

                if (is.null(pred))
                    stop("Prediction returned NULL\n")

                if (tempModel@ctlInfo@classification)
                {
                    foldError[j] <- sum(pred != y[folds[[j]]]) /
                                    length(folds[[j]])

                    if (length(perfParameters) > 0)
                    {
                        foldPerformance <-
                            evaluatePrediction(pred, y[folds[[j]]],
                                               allLabels=tempModel@classNames,
                                               decValues=predDecValues,
                                               print=FALSE)

                        if (collectACC)
                            foldACC[j] <- foldPerformance$ACC / 100

                        if (collectBACC)
                            foldBACC[j] <- foldPerformance$BAL_ACC / 100

                        if (collectMCC)
                            foldMCC[j] <- foldPerformance$MAT_CC

                        if (collectAUC)
                            foldAUC[j] <- foldPerformance$AUC
                    }
                }
                else
                {
                    ## calc MSE
                    foldError[j] <- sum((as.numeric(pred) -
                                    y[folds[[j]]])^2) / length(folds[[j]])
                }

                if (tempModel@svmInfo@selPackage != "LiblineaR")
                {
                    noOfSVs[j] <-
                        length(getSVMSlotValue("svIndex", tempModel))

                    ## not yet implemented
                    if (tempModel@svmInfo@selSVM != "mc-natW")
                    {
                        sumAlphas[j] <-
                            sum(abs(getSVMSlotValue("coef", tempModel)))
                    }
                    else
                        sumAlphas[j] <- NA
                }
                else
                {
                    noOfSVs[j] <- NA
                    sumAlphas[j] <- NA
                }
            }
            else
            {
                foldError[j] <- NA
                noOfSVs[j] <- NA
                sumAlphas[j] <- NA
            }
        }

        nonNA <- sum((!is.na(foldError)))

        if (nonNA > 0)
        {
            model@cvResult@cvError[[i]] <- sum(foldError, na.rm=TRUE) /
                                           nonNA
            model@cvResult@foldErrors[i,] <- foldError
            model@cvResult@noSV[((i-1)*numFolds+1):
                                ((i)*numFolds)] <- noOfSVs
            model@cvResult@sumAlphas[((i-1)*numFolds+1):
                                     ((i)*numFolds)] <- sumAlphas
        }
        else
        {
            model@cvResult@cvError[[i]] <- NA
            model@cvResult@foldErrors[i,] <- NA
            model@cvResult@noSV[((i-1)*numFolds+1):
                                ((i)*numFolds)] <- rep(NA, numFolds)
            model@cvResult@sumAlphas[((i-1)*numFolds+1):
                                     ((i)*numFolds)] <- rep(NA, numFolds)
        }

        if (collectACC)
        {
            nonNA <- sum((!is.na(foldACC)))

            if (nonNA > 0)
            {
                model@cvResult@ACC[[i]] <- sum(foldACC, na.rm=TRUE) / nonNA
                model@cvResult@foldACC[i,] <- foldACC
            }
            else
            {
                model@cvResult@ACC[[i]] <- NA
                model@cvResult@foldACC[i,] <- NA
            }
        }

        if (collectBACC)
        {
            nonNA <- sum((!is.na(foldBACC)))

            if (nonNA > 0)
            {
                model@cvResult@BACC[[i]] <- sum(foldBACC, na.rm=TRUE) / nonNA
                model@cvResult@foldBACC[i,] <- foldBACC
            }
            else
            {
                model@cvResult@BACC[[i]] <- NA
                model@cvResult@foldBACC[i,] <- NA
            }
        }

        if (collectMCC)
        {
            nonNA <- sum((!is.na(foldMCC)))

            if (nonNA > 0)
            {
                model@cvResult@MCC[[i]] <- sum(foldMCC, na.rm=TRUE) / nonNA
                model@cvResult@foldMCC[i,] <- foldMCC
            }
            else
            {
                model@cvResult@MCC[[i]] <- NA
                model@cvResult@foldMCC[i,] <- NA
            }
        }
        
        if (collectAUC)
        {
            nonNA <- sum((!is.na(foldAUC)))
            
            if (nonNA > 0)
            {
                model@cvResult@AUC[[i]] <- sum(foldAUC, na.rm=TRUE) / nonNA
                model@cvResult@foldAUC[i,] <- foldAUC
            }
            else
            {
                model@cvResult@AUC[[i]] <- NA
                model@cvResult@foldAUC[i,] <- NA
            }
        }
    }

    return(model)
}

#' @rdname performCrossValidation-methods
#' @aliases
#' performCrossValidation
#' crossValidation
#' CrossValidation
#' cross.validation
#' @title KeBABS Cross Validation
#'
#' @description Perform cross validation as k-fold cross validation,
#' Leave-One-Out cross validation(LOOCV) or grouped cross validation (GCV).
#'
#' @param object a kernel matrix or an explicit representation
#' @param x an optional set of sequences
#' @param y a response vector
#' @param sel sample subset for which cross validation should be performed
#' @param model KeBABS model
#' @param cross an integer value K > 0 indicates that k-fold cross validation
#' should be performed. A value -1 is used for Leave-One-Out (LOO) cross
#' validation. (see above) Default=0
#' @param noCross an integer value larger than 0 is used to specify the number
#' of repetitions for cross validation. This parameter is only relevant if
#' 'cross' is different from 0. Default=1
#' @param groupBy allows a grouping of samples during cross validation. The
#' parameter is only relevant when 'cross' is larger than 1. It is an integer
#' vector or factor with the same length as the number of samples used for
#' training and specifies for each sample to which group it belongs. Samples
#' from the same group are never spread over more than one fold. Grouped
#' cross validation can also be used in grid search for each grid point.
#' Default=NULL
#' @param perfParameters a character vector with one or several values from
#' the set "ACC" , "BACC", "MCC", "AUC" and "ALL". "ACC" stands for accuracy,
#' "BACC" for balanced accuracy, "MCC" for Matthews Correlation Coefficient,
#' "AUC" for area under the ROC curve and "ALL" for all four. This parameter
#' defines which performance parameters are collected in cross validation
#' for display purpose. The summary values are computed as mean of the
#' fold values. AUC computation from pooled decision values requires a
#' calibrated classifier output and is currently not supported. Default=NULL
#' @param verbose boolean value that indicates whether KeBABS should print
#' additional messages showing the internal processing logic in a verbose
#' manner. The default value depends on the R session verbosity option.
#' Default=getOption("verbose")
#'
#' this parameter is not relevant for cross validation because
#' the method \code{performCrossValidation} should not be called directly.
#' Cross validation is performed with the method \code{\link{kbsvm}} and the
#' parameters \code{cross} and \code{numCross} are described there
#' @usage ## kbsvm(......, cross=0, noCross=1, .....)
#'
#' ## please use kbsvm for cross validation and do not call the
#' ## performCrossValidation method directly
#'
#' @details
#'
#' Overview\cr
#'
#' Cross validation (CV) provides an estimate for the generalization
#' performance of a model based on repeated training on different subsets of
#' the data and evaluating the prediction performance on the remaining data
#' not used for training. Dependent on the strategy of splitting the data
#' different variants of cross validation exist. KeBABS implements k-fold cross
#' validation, Leave-One-Out cross validation and Leave-Group-Out cross
#' validation which is a specific variant of k-fold cross validation. Cross
#' validation is invoked with \code{\link{kbsvm}} through setting the
#' parameters \code{cross} and \code{noCross}. It can either
#' be used for a given kernel and specific values of the SVM hyperparameters to
#' compute the cross validation error of a single model or in conjuction with
#' grid search (see \link{gridSearch}) and model selection (see
#' \link{modelSelection}) to determine the performance of multiple models.\cr\cr
#'
#' k-fold Cross Validation and Leave-One-Out Cross Validation(LOOCV)\cr
#'
#' For k-fold cross validation the data is split into k roughly equal sized
#' subsets called folds. Samples are assigned to the folds randomly. In k
#' successive training runs one of the folds is kept in round-robin manner
#' for predicting the performance while using the other k-1 folds together as
#' training data. Typical values for the number of folds k are 5 or 10
#' dependent on the number of samples used for CV. For LOOCV the fold size
#' decreases to 1 and only a single sample is kept as hold out fold for
#' performance prediction requiring the same number of training runs in one
#' cross validation run as the number of sequences used for CV.\cr\cr
#'
#' Grouped Cross Validation (GCV)\cr
#'
#' For grouped cross validation samples are assigned to groups by the
#' user before running cross validation, e.g. via clustering the sequences.
#' The predefined group assignment is passed to CV with the parameter
#' \code{groupBy} in \code{\link{kbsvm}}. GCV is a special version of k-fold
#' cross validation which respects group boundaries by avoiding to distribute
#' samples of one group over multiple folds. In this way the group(s) in the
#' test fold do not occur during training and learning is forced to concentrate
#' on more complex features instead of the simple features splitting the
#' groups. For GCV the parameter cross must be smaller than or equal to the
#' number of groups.\cr\cr
#'
#' Cross Validation Result\cr
#'
#' The cross validation error, which is the average of the predicition errors
#' in all held out folds, is used as an estimate for the generalization error
#' of the model assiciated with the cross validation run. For classification
#' the fraction of incorrectly classified samples and for regression the mean
#' squared error (MSE) is used as prediction error. Multiple cross validation
#' runs can be performed through setting the parameter \code{noCross}. The
#' cross validation result can be extracted from the model object returned by
#' cross validation with the \code{\link{cvResult}} accessor. It contains the
#' mean CV error over all runs, the CV errors of the single runs and the
#' CV error for each fold. The CV result object can be plotted with the method
#' \code{\link{plot}} showing the variation of the CV error for the different
#' runs as barplot. With the parameter \code{perfParameters} in
#' \code{\link{kbsvm}} the accuracy, the balanced accuracy and the Matthews
#' correlation coefficient can be requested as additional performance
#' parameters to be recorded in the CV result object which might be of interest
#' especially for unbalanced datasets.\cr\cr
#'
#' @return cross validation stores the cross validation results in the 
#' KeBABS model object returned by . They can be retrieved with the accessor 
#' \code{\link{cvResult}} returned by \code{\link{kbsvm}}.
#'
#' @seealso \code{\link{kbsvm}}, \code{\link{cvResult}},
#' \code{\link{plot}}
#'
#' @examples
#' ## load transcription factor binding site data
#' data(TFBS)
#' enhancerFB
#' ## select a few samples for training - here for demonstration purpose
#' ## normally you would use 70 or 80% of the samples for training and
#' ## the rest for test
#' ## train <- sample(1:length(enhancerFB), length(enhancerFB) * 0.7)
#' ## test <- c(1:length(enhancerFB))[-train]
#' train <- sample(1:length(enhancerFB), 50)
#' ## create a kernel object for the gappy pair kernel with normalization
#' gappy <- gappyPairKernel(k=1, m=4)
#' ## show details of kernel object
#' gappy
#'
#' ## run cross validation with the kernel on C-svc in LiblineaR for cost=10
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=gappy,
#'                pkg="LiblineaR", svm="C-svc", cost=10, cross=3)
#'
#' ## show cross validation result
#' cvResult(model)
#'
#' \dontrun{
#' ## perform tive cross validation runs
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=gappy,
#'                pkg="LiblineaR", svm="C-svc", cost=10, cross=10, noCross=5)
#'
#' ## show cross validation result
#' cvResult(model)
#'
#' ## plot cross validation result
#' plot(cvResult(model))
#'
#'
#' ## run Leave-One-Out cross validation
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=gappy,
#'                pkg="LiblineaR", svm="C-svc", cost=10, cross=-1)
#'
#' ## show cross validation result
#' cvResult(model)
#'
#' ## run gouped cross validation with full data
#' ## on coiled coil dataset
#' ##
#' ## In this example the groups were determined through single linkage 
#' ## clustering of sequence similarities derived from ungapped heptad-specific 
#' ## pairwise alignment of the sequences. The variable {\tt ccgroup} contains 
#' ## the pre-calculated group assignments for the individual sequences. 
#' data(CCoil)
#' ccseq
#' head(yCC)
#' head(ccgroups)
#' gappyK1M6 <- gappyPairKernel(k=1, m=4)
#' 
#' ## run k-fold CV without groups
#' model <- kbsvm(x=ccseq, y=as.numeric(yCC), kernel=gappyK1M6, 
#' pkg="LiblineaR", svm="C-svc", cost=10, cross=3, noCross=2, 
#' perfObjective="BACC",perfParameters=c("ACC", "BACC"))
#'
#' ## show result without groups
#' cvResult(model)
#'
#' ## run grouped CV
#' model <- kbsvm(x=ccseq, y=as.numeric(yCC), kernel=gappyK1M6, 
#' pkg="LiblineaR", svm="C-svc", cost=10, cross=3, 
#' noCross=2, groupBy=ccgroups, perfObjective="BACC",
#' perfParameters=c("ACC", "BACC"))
#'
#' ## show result with groups
#' cvResult(model)
#'
#' ## For grouped CV the samples in the held out fold are from a group which
#' ## is not present in training on the other folds. The simimar CV error
#' ## with and without groups shows that learning is not just assigning
#' ## labels based on similarity within the groups but is focusing on features
#' ## that are indicative for the class also in the CV without groups. For the
#' ## GCV no information about group membership for the samples in the held
#' ## out fold is present in the model. This example should show how GCV
#' ## is performed. Because of package size limitations no specific dataset is
#' ## available in this package where GCV is necessary.
#' }
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs/}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics}, 31(15):2574-2576, 2015.
#' DOI: \doi{10.1093/bioinformatics/btv176}.
#' @keywords kbsvm
#' @keywords cross validation
#' @keywords grid search
#' @keywords model selection
#' @keywords methods
#'

#' @rdname performCrossValidation-methods
#' @aliases
#' performCrossValidation
#' performCrossValidation,KernelMatrix-method


setMethod("performCrossValidation",signature(object="KernelMatrix"),
          performCrossValidation.KernelMatrix)

performCrossValidation.kernelMatrix <- function(object, x, y, sel, model,
                                                cross, noCross, groupBy,
                                                perfParameters, verbose)
{
#    library(kernlab)

    performCrossValidation(object=as.KernelMatrix(object), x=x, y=y, sel=sel,
                           model=model, cross=cross, noCross=noCross,
                           groupBy=groupBy, perfParameters=perfParameters,
                           verbose=verbose)
}

#setMethod("performCrossValidation",signature(object="kernelMatrix"),
#          performCrossValidation.kernelMatrix)

performCrossValidation.ExplicitRep <- function(object, x, y, sel, model,
                                               cross, noCross, groupBy,
                                               perfParameters, verbose)
{
    ## subset immediately because of cv
    if (length(sel) > 0)
    {
        object <- object[sel,]
        selOrig <- sel
        sel <- NULL
    }
    else
        selOrig <- NULL

    if (model@svmInfo@selPackage == "LiblineaR" &&
        model@svmInfo@reqFeatureType == "quadratic")
    {
        exRepLin <- object
        object <- getExRepQuadratic(object)
    }

    if (cross == -1)
        numFolds <- nrow(object)
    else
    {
        if (cross >= nrow(object))
            stop("'cross' must be smaller than the number of samples\n")

        numFolds <- cross
    }

    if (!is.null(groupBy))
    {
        ## groupBy is a numeric vector or a factor
        numGroups <- length(table(as.numeric(groupBy)))
    }
    else
        numGroups <- 0

    model@cvResult <- new("CrossValidationResult")
    model@cvResult@cross <- cross
    model@cvResult@noCross <- noCross
    model@cvResult@groupBy <- groupBy
    model@cvResult@perfParameters <- perfParameters
    model@cvResult@folds <- as.list(rep(NA_real_, noCross * numFolds))
    model@cvResult@foldErrors <- matrix(NA_real_, nrow=noCross, ncol=numFolds)

    collectACC <- "ACC" %in% perfParameters
    collectBACC <- "BACC" %in% perfParameters
    collectMCC <- "MCC" %in% perfParameters
    collectAUC <- "AUC" %in% perfParameters

    if (collectACC)
        model@cvResult@foldACC <- matrix(Inf, nrow=noCross, ncol=numFolds)

    if (collectBACC)
        model@cvResult@foldBACC <- matrix(Inf, nrow=noCross, ncol=numFolds)

    if (collectMCC)
        model@cvResult@foldMCC <- matrix(Inf, nrow=noCross, ncol=numFolds)

    if (collectAUC)
        model@cvResult@foldAUC <- matrix(Inf, nrow=noCross, ncol=numFolds)

    if (model@ctlInfo@classification == TRUE && model@numClasses == 0)
    {
        ## check for multiclass to get the no of classes reliably
        model@numClasses <- length(unique(y))

        if (model@numClasses > 2)
            model@ctlInfo@multiclassType <- getMulticlassType(model)
    }

    for (i in 1:noCross)
    {
        ## split training data into folds
        randomSamples <- sample(1:nrow(object), nrow(object))

        if (numGroups > 0)
        {
            ## assign full groups to folds
            ## for cross equal to numGroups => Leave Group Out CV
            groupsInFolds <-
                    suppressWarnings(split(sample(1:numGroups), 1:cross))
            grouping <- as.numeric(groupBy)[randomSamples]
            folds <- lapply(1:cross,
                            function(i){randomSamples[which(grouping %in%
                                                        groupsInFolds[[i]])]})
        }
        else
        {
            ## k-fold or Leave One Out CV
            folds <- suppressWarnings(split(randomSamples, 1:numFolds))
        }

        numFolds <- length(folds)
        model@cvResult@folds[((i-1)*numFolds+1):((i)*numFolds)] <- folds
        noOfSVs <- rep(NA_real_, numFolds)
        sumAlphas <- rep(NA_real_, numFolds)
        foldError <- rep(NA_real_, numFolds)
        predDecValues <- NA

        if (collectACC)
            foldACC <- rep(NA_real_, numFolds)

        if (collectBACC)
            foldBACC <- rep(NA_real_, numFolds)

        if (collectMCC)
            foldMCC <- rep(NA_real_, numFolds)

        if (collectAUC)
            foldAUC <- rep(NA_real_, numFolds)

        for (j in 1:numFolds)
        {
            trainSuccess <- TRUE
            tempModel <- model
            trainIndices <- unlist(folds[-j], use.names=FALSE)

            allArgs <- tempModel@svmInfo@selSVMPar
            allArgs[["x"]] <- object[trainIndices,]
            allArgs[["y"]] <- y[trainIndices]
            allArgs[["svmInfo"]] <- tempModel@svmInfo
            allArgs[["verbose"]] <- verbose

            tempModel@svmModel <- tryCatch({do.call(trainSVM, allArgs)},
                warning =
                  function(w)
                  {
                    message(paste("\nwarning during cross validation for fold ",
                                  j,":"))
                    message(w)
                  },
                error =
                  function(e)
                  {
                    message(paste("\nerror during cross validation for fold ",
                                  j,":"))
                    message(e)
                    trainSuccess <- FALSE
                  })

            if (!is.logical(tempModel@svmModel))
            {
                if (is.null(tempModel@svmModel))
                {
                    stop(paste(tempModel@svmInfo@selPackage, "-",
                        tempModel@svmInfo@selSVM, "did not return a model\n"))
                }

                if (tempModel@ctlInfo@classification == TRUE)
                {
                    tempModel@numClasses <-
                        getSVMSlotValue("numClasses", tempModel)
                    tempModel@classNames <-
                        getSVMSlotValue("classNames", tempModel)
                }

                svIndices <- getSVMSlotValue("svIndex", tempModel)

                if (length(svIndices) > 0)
                {
                    alphaIndices <- getSVMSlotValue("alphaIndex", tempModel)

                    if (length(selOrig) > 0)
                    {
                        tempModel@SV <- x[selOrig[trainIndices[svIndices]]]
                        tempModel@svIndex <- selOrig[trainIndices[svIndices]]

                        if (is.list(alphaIndices))
                        {
                            tempModel@alphaIndex <- lapply(alphaIndices,
                                function(x) selOrig[trainIndices[x]])
                        }
                        else
                        {
                            tempModel@alphaIndex <-
                                selOrig[trainIndices[alphaIndices]]
                        }
                    }
                    else
                    {
                        tempModel@SV <- x[trainIndices[svIndices]]
                        tempModel@svIndex <- trainIndices[svIndices]

                        if (is.list(alphaIndices))
                        {
                            tempModel@alphaIndex <- lapply(alphaIndices,
                                        function(x) trainIndices[x])
                        }
                        else
                            tempModel@alphaIndex <- trainIndices[alphaIndices]
                    }

                    exRepSV <- object[trainIndices[svIndices],]
                }
                else
                {
                    ## $$$ TODO remove
                    if (tempModel@svmInfo@selPackage != "LiblineaR")
                    {
                        print("no support vectors found\n")
                        next
                    }

                    exRepSV <- NULL
                }

                if (model@svmInfo@featureWeights != "no")
                {
                    tempModel@featureWeights <-
                        getFeatureWeights(model=tempModel,
                                    exrep=exRepSV,
                                    weightLimit=tempModel@svmInfo@weightLimit)

                    tempModel@b <- getSVMSlotValue("b", tempModel)
                    tempModel@SV <- NULL

                }

                if (model@svmInfo@selPackage == "LiblineaR" &&
                    model@svmInfo@reqFeatureType == "quadratic")
                {
                    pred <- predict(object=tempModel, x=exRepLin[folds[[j]],],
                                    predictionType="response", verbose=verbose)

                    if (collectAUC)
                    {
                        predDecValues <- predict(object=tempModel,
                                                 x=exRepLin[folds[[j]],],
                                                 predictionType="decision",
                                                 verbose=verbose)
                    }
                }
                else
                {
                    pred <- predict(object=tempModel, x=object[folds[[j]],],
                                    predictionType="response", verbose=verbose)
                                    
                    if (collectAUC)
                    {
                        predDecValues <- predict(object=tempModel,
                                                 x=object[folds[[j]],],
                                                 predictionType="decision",
                                                 verbose=verbose)
                    }
                }

                if (is.null(pred))
                    stop("Prediction returned NULL\n")

                if (tempModel@ctlInfo@classification)
                {
                    foldError[j] <- sum(pred != y[folds[[j]]]) /
                                    length(folds[[j]])

                    if (length(perfParameters) > 0)
                    {
                        foldPerformance <-
                            evaluatePrediction(pred, y[folds[[j]]],
                                               allLabels=tempModel@classNames,
                                               decValues=predDecValues,
                                               print=FALSE)

                        if (collectACC)
                            foldACC[j] <- foldPerformance$ACC / 100

                        if (collectBACC)
                            foldBACC[j] <- foldPerformance$BAL_ACC / 100

                        if (collectMCC)
                            foldMCC[j] <- foldPerformance$MAT_CC
                            
                        if (collectAUC)
                            foldAUC[j] <- foldPerformance$AUC
                    }
                }
                else
                {
                    ## calc MSE
                    foldError[j] <- sum((as.numeric(pred) -
                                      y[folds[[j]]])^2) / length(folds[[j]])
                }

                if (tempModel@svmInfo@selPackage != "LiblineaR")
                {
                    noOfSVs[j] <-
                        length(getSVMSlotValue("svIndex", tempModel))

                    ## not yet implemented
                    if (tempModel@svmInfo@selSVM != "mc-natW")
                    {
                        sumAlphas[j] <-
                            sum(abs(getSVMSlotValue("coef", tempModel)))
                    }
                    else
                        sumAlphas[j] <- NA
                }
                else
                {
                    noOfSVs[j] <- NA
                    sumAlphas[j] <- NA
                }
            }
            else
            {
                foldError[j] <- NA
                noOfSVs[j] <- NA
                sumAlphas[j] <- NA
            }
        }

        nonNA <- sum((!is.na(foldError)))

        if (nonNA > 0)
        {
            model@cvResult@cvError[[i]] <- sum(foldError, na.rm=TRUE) / nonNA
            model@cvResult@foldErrors[i,] <- foldError

            model@cvResult@noSV[((i-1)*numFolds+1):
                                ((i)*numFolds)] <- noOfSVs
            model@cvResult@sumAlphas[((i-1)*numFolds+1):
                                     ((i)*numFolds)] <- sumAlphas
        }
        else
        {
            model@cvResult@cvError[[i]] <- NA
            model@cvResult@foldErrors[i,] <- NA

            model@cvResult@noSV[((i-1)*numFolds+1):
                                ((i)*numFolds)] <- rep(NA, numFolds)
            model@cvResult@sumAlphas[((i-1)*numFolds+1):
                                     ((i)*numFolds)] <- rep(NA, numFolds)
        }

        if (collectACC)
        {
            nonNA <- sum((!is.na(foldACC)))

            if (nonNA > 0)
            {
                model@cvResult@ACC[[i]] <- sum(foldACC, na.rm=TRUE) / nonNA
                model@cvResult@foldACC[i,] <- foldACC
            }
            else
            {
                model@cvResult@ACC[[i]] <- NA
                model@cvResult@foldACC[i,] <- NA
            }
        }

        if (collectBACC)
        {
            nonNA <- sum((!is.na(foldBACC)))

            if (nonNA > 0)
            {
                model@cvResult@BACC[[i]] <- sum(foldBACC, na.rm=TRUE) / nonNA
                model@cvResult@foldBACC[i,] <- foldBACC
            }
            else
            {
                model@cvResult@BACC[[i]] <- NA
                model@cvResult@foldBACC[i,] <- NA
            }
        }

        if (collectMCC)
        {
            nonNA <- sum((!is.na(foldMCC)))

            if (nonNA > 0)
            {
                model@cvResult@MCC[[i]] <- sum(foldMCC, na.rm=TRUE) / nonNA
                model@cvResult@foldMCC[i,] <- foldMCC
            }
            else
            {
                model@cvResult@MCC[[i]] <- NA
                model@cvResult@foldMCC[i,] <- NA
            }
        }
        
        if (collectAUC)
        {
            nonNA <- sum((!is.na(foldAUC)))
            
            if (nonNA > 0)
            {
                model@cvResult@AUC[[i]] <- sum(foldAUC, na.rm=TRUE) / nonNA
                model@cvResult@foldAUC[i,] <- foldAUC
            }
            else
            {
                model@cvResult@AUC[[i]] <- NA
                model@cvResult@foldAUC[i,] <- NA
            }
        }
    }

    return(model)
}

#' @rdname performCrossValidation-methods
#' @aliases
#' performCrossValidation,ExplicitRepresentation-method
#'

setMethod("performCrossValidation",signature(object="ExplicitRepresentation"),
          performCrossValidation.ExplicitRep)

