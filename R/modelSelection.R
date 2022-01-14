##345678901234567890123456789012345678901234567890123456789012345678901234567890
#' @rdname performModelSelection
#' @aliases
#' modelSelection
#' ModelSelection
#' model.selection
#' @title KeBABS Model Selection
#'
#' @description Perform model selection with one or multiple sequence kernels
#' on one or multiple SVMs with one or multiple SVM parameter sets.
#'
#' @param nestedCross for this and other parameters see \code{\link{kbsvm}}
#' @usage
#' ## kbsvm(...., kernel=..., pkg=..., svm=..., cost=..., ....,
#' ##       cross=0, noCross=1, ...., nestedCross=0, noNestedCross=1, ....)
#'
#' ## For details see below. With parameter nestedCross > 1 model selection is
#' ## performed, the other parameters are handled identical to grid search.
#'
#' @details
#'
#' Overview\cr
#'
#' Model selection in KeBABS is based on nested k-fold cross validation (CV)
#' (for details see \link{performCrossValidation}). The inner cross
#' validation is used to determine the best parameters settings (kernel
#' parameters and SVM parameters) and the outer cross validation to verify
#' the performance on data that was not included in the selection of the
#' best model. The training folds of the outer CV are used to run a grid
#' search with the inner cross validation running for each point of the
#' grid (see \code{\link{performGridSearch}} to find the best performing model.
#' Once this model is selected the performance of this model on the held out
#' fold of the outer CV is determined. Different model parameters settings
#' could occur for different held out folds of the outer CV. This means that
#' model selection does not deliver a performance estimate for a single
#' best model but for the complete model selection process.\cr
#'
#' For each run of the outer CV KeBABS stores the selected parameter setting
#' for the best performing model. The default performance objective for
#' selecting the best parameters setting is based on minimizing the CV error
#' on the inner CV. With the parameter \code{perfObjective} in
#' \code{\link{kbsvm}} the balanced accuracy or the Matthews correlation
#' coefficient can be used instead for which the parameter setting with the
#' maximal value is selected. The parameter setting of the best performing
#' model for each fold in the outer CV can be retrieved from the KeBABS model
#' with the accessor \code{\link{modelSelResult}}. The performance values on
#' the outer CV are retrieved from the model with the accessor
#' \code{\link{cvResult}}.\cr
#'
#' Model selection is invoked through the method \code{\link{kbsvm}} through
#' setting parameter \code{nestedCross} > 1. For the parameters \code{kernel,
#' pkg, svm} and SVM hyperparameters the handling is identical to grid search
#' (see \code{\link{performGridSearch}}). The parameter cost in the usage
#' section above is just one representative of SVM hyperparameters to indicate
#' their relevance for model selection. The complete model selection process
#' can be repeated multiple times through setting \code{noNestedCross} to the
#' number of desired repetitions. Nested cross validation used in model
#' selection is dynamically more demanding than grid search. Concerning runtime
#' please see the runtime hints for \code{\link{performGridSearch}}.\cr
#'
#' @return model selection stores the results in the KeBABS model. They can be
#' retrieved with the accessor \code{\link{modelSelResult}{KBModel}}. Results
#' from the outer cross validation are extracted from the model with the
#' accessor\code{\link{cvResult}}.
#'
#' @seealso \code{\link{kbsvm}}, \code{\link{performGridSearch}},
#' \code{\link{modelSelResult}},
#' \code{\link{cvResult}}
#'
#' @examples
#'
#' ## load transcription factor binding site data
#' data(TFBS)
#' enhancerFB
#' ## The C-svc implementation from LiblineaR is chosen for most of the
#' ## examples because it is the fastest SVM. With SVMs from other packages
#' ## slightly better results could be achievable. Because of the higher
#' ## runtime needed for nested cross validation please run the examples
#' ## below manually. All samples of the data set are used in the examples.
#' train <- sample(1:length(enhancerFB), length(enhancerFB))
#'
#' ## model selection with single kernel object and multiple
#' ## hyperparameter values, 5 fold inner CV and 3 fold outer CV
#' ## create gappy pair kernel with normalization
#' gappyK1M3 <- gappyPairKernel(k=1, m=3)
#' ## show details of single gappy pair kernel object
#' gappyK1M3
#'
#' pkg <- "LiblineaR"
#' svm <- "C-svc"
#' cost <- c(50,100,150,200,250,300)
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=gappyK1M3,
#'                pkg=pkg, svm=svm, cost=cost, explicit="yes", cross=3,
#'                nestedCross=2, showProgress=TRUE)
#'
#' ## show best parameter settings
#' modelSelResult(model)
#'
#' ## show model selection result which is the result of the outer CV
#' cvResult(model)
#
#' \dontrun{
#' ## repeated model selection
#' pkg <- "LiblineaR"
#' svm <- "C-svc"
#' cost <- c(50,100,150,200,250,300)
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=gappyK1M3,
#'                pkg=pkg, svm=svm, cost=cost, explicit="yes", cross=10,
#'                nestedCross=3, noNestedCross=3, showProgress=TRUE)
#'
#' ## show best parameter settings
#' modelSelResult(model)
#'
#' ## show model selection result which is the result of the outer CV
#' cvResult(model)
#'
#' ## plot CV result
#' plot(cvResult(model))
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
#' @keywords grid search
#' @keywords model selection
#' @keywords methods

performModelSelection <- function(object, model, y, explicit, featureWeights,
                                  weightLimit, sel, features=NULL, groupBy,
                                  perfParameters, perfObjective, addArgs,
                                  showProgress, showCVTimes, verbose)
{
    if (is(object, "kernelMatrix"))
        object <- as.KernelMatrix(object)

    if (model@modelSelResult@nestedCross == 0)
    {
        return(performGridSearch(object=object, model=model, y=y,
                        explicit=explicit, featureWeights=featureWeights,
                        weightLimit=weightLimit, sel=sel, features=features,
                        groupBy=groupBy, perfParameters=perfParameters,
                        perfObjective=perfObjective, includeFullModel=FALSE,
                        showProgress=showProgress, showCVTimes=showCVTimes,
                        addArgs=addArgs, verbose=verbose))
    }

    if (is(object, "KernelMatrix"))
    {
        stop("model selection with kernel matrix only is not\n",
             "        supported\n")
    }

    ## subset immediately to avoid repeated subsetting
    ## y is subsetted on caller side
    if (length(sel) > 0)
    {
        object <- subsetSeqRep(x=object, sel=sel)
        sel <- integer(0)
    }

    numSamples <- getNoOfElementsOfSeqRep(object)

    if (length(y) != numSamples)
        stop("length of y not consistent with number of samples\n")

    if (model@modelSelResult@cross >= numSamples)
        stop("'cross' must be smaller than the number of samples\n")

    model@svmInfo@selPackage <- model@svmInfo@reqPackage
    model@svmInfo@selSVM <- model@svmInfo@reqSVM
    model@ctlInfo@classification <-
        isClassification(model@svmInfo@selPackage, model@svmInfo@selSVM)

    singleKernel <- NULL
    precompTestdata <- NULL
    explicitType <- model@svmInfo@reqExplicitType

    if (length(model@modelSelResult@gridRows) < 2)
    {
        precompute <- TRUE

        if (length(model@modelSelResult@gridRows) == 1)
            model@svmInfo@selKernel <- model@modelSelResult@gridRows[[1]]
        else
            model@svmInfo@selKernel <- model@svmInfo@reqKernel

        singleKernel <- model@svmInfo@selKernel
        
        if (isUserDefined(singleKernel))
            explicit <- "no"
    }

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

    ## preload SparseM to avoid loading message to interfere with
    ## progress messages
    if (explicitType %in% c("auto","sparse"))
    {
        if (!requireNamespace("SparseM", quietly=TRUE))
            stop("package SparseM could not be loaded\n")
    }

    if ((explicit != "no") && is(object, "KernelMatrix"))
        stop("processing via explicit representation is not possible\n")

    if (explicit != "no")
    {
        model@svmInfo@selExplicit <- TRUE
        model@ctlInfo@selMethod <- "explicitRep"

        if (!is.null(singleKernel))
        {
            if ((is(object, "BioVector") || is(object, "XStringSet")))
            {
                precompTestdata <- getExRep(x=object, kernel=singleKernel,
                                            sparse=(explicitType=="sparse"),
                                            features=features)
            }
            else if(is(object, "ExplicitRepresentation"))
                precompTestdata <- object
        }
        else
        {
            if (length(features) > 0)
            {
                stop("feature subsetting is not supported for multiple",
                     " kernels\n")
            }
        }
    }
    else
    {
        model@svmInfo@selExplicit <- FALSE
        model@ctlInfo@selMethod <- "KernelMatrix"

        if (!is.null(singleKernel))
        {
            if ((is(object, "BioVector") || is(object, "XStringSet")))
                precompTestdata <- singleKernel(x=object)
            else if(is(object, "KernelMatrix"))
                precompTestdata <- object
        }
    }

    if (model@svmInfo@selExplicit)
        model@svmInfo@explicitKernel <- model@svmInfo@reqFeatureType

    ## $$$ TODO stratified CV and balanced stratified CV

    if (model@modelSelResult@nestedCross == -1)
        numFolds <- numSamples
    else
        numFolds <- model@modelSelResult@nestedCross

    if (!is.null(groupBy))
    {
        ## groupBy is a numeric vector or a factor
        numGroups <- length(table(as.numeric(groupBy)))
    }
    else
        numGroups <- 0

    model@cvResult <- new("CrossValidationResult")
    model@cvResult@outerCV <- TRUE
    model@cvResult@cross <- model@modelSelResult@nestedCross
    model@cvResult@noCross <- model@modelSelResult@noNestedCross
    model@cvResult@perfParameters <- perfParameters
    model@cvResult@groupBy <- groupBy
    model@cvResult@folds <- as.list(rep(NA_real_,
                                        model@cvResult@noCross * numFolds))
    model@cvResult@foldErrors <-
                matrix(Inf, nrow=model@modelSelResult@noNestedCross,
                       ncol=numFolds)
    model@cvResult@noSV <- numeric(0)
    model@cvResult@sumAlphas <- numeric(0)

    collectACC <- "ACC" %in% perfParameters
    collectBACC <- "BACC" %in% perfParameters
    collectMCC <- "MCC" %in% perfParameters
    collectAUC <- "AUC" %in% perfParameters

    if (collectACC)
        model@cvResult@foldACC <- matrix(Inf,
                   nrow=model@modelSelResult@noNestedCross, ncol=numFolds)

    if (collectBACC)
        model@cvResult@foldBACC <- matrix(Inf,
                   nrow=model@modelSelResult@noNestedCross, ncol=numFolds)

    if (collectMCC)
        model@cvResult@foldMCC <- matrix(Inf,
                   nrow=model@modelSelResult@noNestedCross, ncol=numFolds)

    if (collectAUC)
        model@cvResult@foldAUC <- matrix(Inf,
                   nrow=model@modelSelResult@noNestedCross, ncol=numFolds)

    if (model@ctlInfo@classification == TRUE && model@numClasses == 0)
    {
        ## check for multiclass to get the no of classes reliably
        model@numClasses <- length(unique(y))

        if (model@numClasses > 2)
        model@ctlInfo@multiclassType <- getMulticlassType(model)
    }


    for (i in 1:model@modelSelResult@noNestedCross)
    {
        if (showProgress && model@modelSelResult@noNestedCross > 1)
            message("\nNested CV Run:", i)

        ## split training data into folds
        randomSamples <- sample(1:numSamples, numSamples)

        if (numGroups > 0)
        {
            ## assign full groups to folds
            ## for nestedCross equal to numGroups => Leave Group Out CV
            groupsInFolds <-
                suppressWarnings(split(sample(1:numGroups),
                                        1:model@modelSelResult@nestedCross))
            grouping <- as.numeric(groupBy)[randomSamples]
            folds <- lapply(1:model@modelSelResult@nestedCross,
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

        model@cvResult@folds[((i-1)*numFolds+1):(i*numFolds)] <- folds
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

        for (j in 1:length(folds))
        {
            tempModel <- model
            modelSelIndices <- unlist(folds[-j], use.names=FALSE)

            if (!is.null(groupBy))
            {
                ## remove unused factor levels
                tempGroupBy <- factor(groupBy[modelSelIndices], exclude=NULL)
                tempModel@modelSelResult@groupBy <- tempGroupBy

                tempNumGroups <- length(table(as.numeric(tempGroupBy)))

                if (tempNumGroups < tempModel@modelSelResult@cross)
                {
                 stop("the number of groups must be larger than or equal to\n",
                         "        the number of folds\n")
                }
            }
            else
                tempGroupBy <- NULL

            selModel <- performGridSearch(
                            object=subsetSeqRep(x=object, sel=modelSelIndices),
                            model=tempModel, y=y[modelSelIndices],
                            explicit=explicit, featureWeights=featureWeights,
                            weightLimit=weightLimit, sel=integer(0),
                            features=features, groupBy=tempGroupBy,
                            perfParameters=perfParameters,
                            perfObjective=perfObjective,
                            showProgress=showProgress, includeFullModel=TRUE,
                            showCVTimes=showCVTimes, addArgs=addArgs,
                            verbose=verbose)

            ## store selected grid row and grid col
            if (length(model@modelSelResult@gridRows) > 0)
            {
                model@modelSelResult@selGridRow[
                    length(model@modelSelResult@selGridRow)+1] <-
                                list(selModel@modelSelResult@selGridRow)
            }

            if (!is.null(model@modelSelResult@gridCols) &&
                (nrow(model@modelSelResult@gridCols) > 0))
            {
                if (length(model@modelSelResult@selGridCol) > 0)
                {
                    model@modelSelResult@selGridCol <-
                        rbind(model@modelSelResult@selGridCol,
                            data.frame(selGridCol=
                                rownames(selModel@modelSelResult@selGridCol),
                                selModel@modelSelResult@selGridCol))
                }
                else
                {
                    model@modelSelResult@selGridCol <-
                        data.frame(selGridCol=
                            rownames(selModel@modelSelResult@selGridCol),
                            selModel@modelSelResult@selGridCol)
                }

                rownames(model@modelSelResult@selGridCol) <- NULL
            }

            if (!is.null(singleKernel))
            {
                if (is(precompTestdata, "KernelMatrix"))
                {
                    svInd <- getSVMSlotValue("svIndex",
                                             selModel@modelSelResult@fullModel)

                    if (length(svInd) < 1)
                        stop("model selection via kernel matrix is not\n",
                             "        supported\n")

                    if (model@svmInfo@selPackage == "e1071")
                    {
                        testData <- new("KernelMatrix", .Data=
                                        matrix(0, nrow=length(folds[[j]]),
                                               ncol=length(modelSelIndices)))
                        testData@.Data[,svInd] <-
                            precompTestdata[folds[[j]],
                                            modelSelIndices[svInd],
                                            drop=FALSE]
                    }
                    else
                    {
                        testData <- precompTestdata[folds[[j]],
                                                    modelSelIndices[svInd],
                                                    drop=FALSE]
                    }
                }
                else
                    testData <- precompTestdata[folds[[j]],]
            }
            else
            {
                selKernel <- selModel@modelSelResult@selGridRow
                testData <- subsetSeqRep(x=object, sel=folds[[j]])
                sparse <- selModel@modelSelResult@fullModel@ctlInfo@sparse

                if (selModel@modelSelResult@fullModel@svmInfo@selExplicit)
                {
                    if (is(object, "BioVector") || is(object, "XStringSet"))
                    {
                        ## call predict with linear ex rep also for quadratic
                        testData <- getExRep(x=testData,
                                             kernel=selKernel,
                                             sparse=sparse)
                    }
                }
                else
                {
                    svInd <- getSVMSlotValue("svIndex",
                                          selModel@modelSelResult@fullModel)

                    if (length(svInd) < 1)
                        stop("model selection via kernel matrix is not\n",
                             "        supported\n")

                    if (is(object, "BioVector") || is(object, "XStringSet"))
                    {
                        testData <- selKernel(x=testData,
                                        y=object[modelSelIndices[svInd]])
                    }
                    else if (inherits(object, "ExplicitRepresentation"))
                    {
                        testData <- linearKernel(x=testData,
                                        y=object[modelSelIndices[svInd]])
                    }
                }
            }

            pred <- predict(object=selModel@modelSelResult@fullModel,
                            x=testData, predictionType="response",
                            verbose=verbose)

            if (collectAUC)
            {
                predDecValues <-
                    predict(object=selModel@modelSelResult@fullModel,
                            x=testData, predictionType="decision",
                            verbose=verbose)
            }

            if (is.null(pred))
                stop("Prediction returned NULL\n")

            if (selModel@modelSelResult@fullModel@ctlInfo@classification)
            {
                foldError[j] <- sum(pred != y[folds[[j]]]) /
                                    length(folds[[j]])

                if (length(perfParameters) > 0)
                {
                    foldPerformance <-
                        evaluatePrediction(pred, y[folds[[j]]],
                            allLabels=
                            selModel@modelSelResult@fullModel@classNames,
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
        }

        nonNA <- sum((!is.na(foldError)))

        if (nonNA > 0)
        {
            model@cvResult@cvError[[i]] <-
                     sum(foldError, na.rm=TRUE) / nonNA
            model@cvResult@foldErrors[i,] <- foldError
        }
        else
        {
            model@cvResult@cvError[[i]] <- NA
            model@cvResult@foldErrors[i,] <- NA
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

