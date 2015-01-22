##345678901234567890123456789012345678901234567890123456789012345678901234567890

unifyPredictionOutput <- function(model, pred, predictionType)
{
    ## pred is already converted in low level predict method
    if (is.matrix(pred) && ncol(pred) == 1)
        pred <- as.numeric(pred)

    if (predictionType == "response")
    {
        if (is.integer(model@levels) && !is.integer(pred))
        {
            if (is.factor(pred))
                pred <- as.integer(levels(pred)[pred])
            else
            {
                if (!(is.numeric(pred) && !model@ctlInfo@classification))
                    warning("prediction values are of unexpected type\n")
            }
        }
        else if (is.numeric(model@levels) && !is.numeric(pred))
        {
            if (is.factor(pred))
                pred <- as.numeric(levels(pred)[pred])
            else
                warning("prediction values are of unexpected type\n")
        }
        else if (is.character(model@levels) && !is.character(pred))
        {
            if (is.factor(pred))
                pred <- levels(pred)[pred]
            else
                warning("prediction values are of unexpected type\n")
        }
        else if (is.factor(model@levels) && !is.factor(pred))
        {
            if (is.character(pred))
                pred <- factor(pred, levels=model@levels)
            else
                warning("prediction values are of unexpected type\n")
        }
    }
    else if (predictionType == "decision")
    {
        if (!(is.matrix(pred) || is.numeric(pred)))
            warning("decision values are of unexpected type\n")

        if (model@svmInfo@selPackage == "e1071" &&
                 model@ctlInfo@classification &&
                 model@svmModel$nclasses == 2 &&
                 model@svmModel$labels[1] == 1)
            pred <- - pred
        else if (model@svmInfo@selPackage == "LiblineaR" &&
                 model@svmModel$NbClass == 2 &&
                 model@ctlInfo@classification)
        {
            pred <- pred[,1]
            classNamesRaw <- getSVMSlotValue("classNames", model, raw=TRUE)

            if (any(sort(classNamesRaw, decreasing=TRUE) != classNamesRaw))
                pred <- - pred
        }
    }

    return(pred)
}

## Determination of probabilities from decision values via fitted logistic
## function according to implementation in e1071 in vectorized form
getProbability <- function(decisionValue, A, B)
{
    fApB = decisionValue * A + B
    probabilities <- rep(0, length(decisionValue))

    ## avoid catastrophic cancellation
    ## see H.T. Lin, C.J. Lin, R.C Weng, A Note on Platts Probabilistic
    ## Outputs for Support Vector Machines
    useAlt <- fApB >= 0
    probabilities[useAlt] <- exp(-fApB[useAlt]) / (1 + exp(-fApB[useAlt]))
    probabilities[!useAlt] <- 1 / (1 + exp(-fApB[!useAlt]))
    return(probabilities)
}

predict.PositionDependent <- function(model, x, predictionType, sel,
                                      verbose, ...)
{
    if (verbose)
    {
        classifierType <- kebabsInfo@classifierMap[model@svmInfo@selSVM,
                                                   model@svmInfo@selPackage]

        verbM(paste("predict - kbsvm with position specific feature weights:"),
              classifierType, list(...))
    }

    offsetX <- mcols(x)[["offset"]]

    if (is.null(offsetX))
    {
        minPos <- 1
        maxPos <- max(width(x))
        offsetX <- integer(0)
    }
    else
    {
        minPos <- min(offsetX)
        maxPos <- max(offsetX + width(x) - 1)
    }

    if (!is(model@svmInfo@selKernel, "MotifKernel"))
    {
        k <- kernelParameters(model@svmInfo@selKernel)$k
        motifs <- NULL
        motifLengths <- NULL
        maxMotifLength <- 0
        maxPatternLength <- 0
        nodeLimit <- 0
    }
    else
    {
        motifs <- kernelParameters(model@svmInfo@selKernel)$motifs
        motifLengths <- kernelParameters(model@svmInfo@selKernel)$motifLengths
        maxMotifLength <- max(motifLengths)
        maxPatternLength <- max(nchar(motifs))

        ## rough limit for no of nodes in motif tree from no of
        ## chars and no of substitution groups, add one for root
        nodeLimit <- sum(motifLengths) + 1 +
                     sum(sapply(gregexpr("[", motifs, fixed=TRUE),
                                function(x) length(unlist(x))))
        k <- 0
    }


    if (is(model@svmInfo@selKernel, "MismatchKernel") ||
        is(model@svmInfo@selKernel, "GappyPairKernel"))
        m <- kernelParameters(model@svmInfo@selKernel)$m
    else
        m <- 0

    kernelType <- 1:6
    names(kernelType) <- c("SpectrumKernel",
                           "MixedSpectrumKernel",
                           "MismatchKernel",
                           "MotifKernel",
                           "WeightedDegreeKernel",
                           "GappyPairKernel")

    distWeight <- kernelParameters(model@svmInfo@selKernel)$distWeight
    maxDist <- maxPos - minPos

    if (is.function(distWeight))
    {
        ## precompute distance weight vector
        ## terminate on stop and warning
        ## assuming that all distances are partially overlapping
        distWeight <- tryCatch(distWeight(0:(maxDist - k + 1)),
                               warning=function(w) {stop(w)},
                               error=function(e) {stop(e)})

        if (!(is.numeric(distWeight) && length(distWeight) ==
              maxDist - k + 2))
        {
            stop("distWeight function did not return a numeric vector\n",
                 "       of correct length\n")
        }

        ## limit to values larger than .Machine$double.eps
        ## for non-monotonic decreasing functions search from end
        for (i in (maxDist - k + 1):0)
        {
            if (distWeight[i] > .Machine$double.eps)
                break
        }

        distWeight <- distWeight[0:i]
    }

    if (isTRUE(all.equal(distWeight, c(1, rep(0, length(distWeight)-1)))))
    {
        posSpec <- TRUE
        distWeight <- numeric(0)
    }
    else
        posSpec <- FALSE

    normalized <- kernelParameters(model@svmInfo@selKernel)$normalized
    ignoreLower <- kernelParameters(model@svmInfo@selKernel)$ignoreLower
    reverseComplement <- kernelParameters(model@svmInfo@selKernel)$revComplement
    maxSeqLength <- max(width(x))
    bioCharset <- getBioCharset(x, TRUE)
    unmapped <- is(x, "DNAStringSet") || is(x, "RNAStringSet")
    isXStringSet <- inherits(x, "XStringSet")

    if (length(model@ctlInfo@multiclassType) > 0 &&
        model@ctlInfo@multiclassType == "pairwise")
    {
        pred <- matrix(NA, length(sel), choose(model@numClasses, 2))

        for (i in 1:choose(model@numClasses, 2))
        {
            pos1 <- as.numeric(colnames(model@featureWeights[[i]])[1])

            pred[,i] <-
                .Call("getPosDepPredOrProfC", model@featureWeights[[i]],
                      as.double(model@svmInfo@weightLimit), model@b[i], x,
                      is.logical(isXStringSet), length(sel), sel - 1, offsetX,
                      as.integer(maxSeqLength), as.integer(bioCharset[[2]]),
                      as.integer(kernelType[class(model@svmInfo@selKernel)]),
                      as.integer(k), as.integer(m), motifs, motifLengths,
                      as.integer(maxMotifLength), as.integer(maxPatternLength),
                      as.integer(nodeLimit), as.logical(posSpec), distWeight,
                      as.logical(ignoreLower), as.logical(unmapped),
                      as.logical(reverseComplement), as.logical(normalized),
                      as.logical(FALSE), as.integer(pos1), as.integer(minPos),
                      as.integer(maxPos))

            if (model@ctlInfo@classification == TRUE)
            {
                if (predictionType == "response")
                {
                    pred[,i] <- 1.5 - 0.5 * sign(pred[,i])
                    pred <- model@classNames[pred[,i]]

                    if (is.integer(model@levels))
                        pred[,i]  <- as.integer(pred[,i])
                    else if (is.factor(model@levels))
                        pred[,i] <- factor(pred[,i], levels=model@levels)
                }
                else if (predictionType == "probabilities")
                {
                    pred[,i] <- getProbability(pred[,i], model@probA,
                                               model@probB)
                }
            }
        }

        ## $$$ TODO combine result from pairwise classifiers
        ##          check how this is handled in kernlab, e1071 and LiblineaR
    }
    else
    {
        pos1 <- as.numeric(colnames(model@featureWeights)[1])

        pred <- .Call("getPosDepPredOrProfC", model@featureWeights,
                      as.double(model@svmInfo@weightLimit), model@b,
                      x, is.logical(isXStringSet),
                      sel - 1, as.integer(length(sel)), offsetX,
                      as.integer(maxSeqLength), as.integer(bioCharset[[2]]),
                      as.integer(kernelType[class(model@svmInfo@selKernel)]),
                      as.integer(k), as.integer(m), motifs, motifLengths,
                      as.integer(maxMotifLength), as.integer(maxPatternLength),
                      as.integer(nodeLimit), as.logical(posSpec), distWeight,
                      as.logical(ignoreLower), as.logical(unmapped),
                      as.logical(reverseComplement), as.logical(normalized),
                      as.logical(FALSE), as.integer(pos1), as.integer(minPos),
                      as.integer(maxPos))

        if (model@ctlInfo@classification == TRUE)
        {
            if (predictionType == "response")
            {
                pred <- 1.5 - 0.5 * sign(pred)
                pred <- model@classNames[pred]

                if (is.integer(model@levels))
                    pred  <- as.integer(pred)
                else if (is.factor(model@levels))
                    pred <- factor(pred, levels=model@levels)
            }
            else if (predictionType == "probabilities")
                pred <- getProbability(pred, model@probA, model@probB)
        }
    }

    return(pred);
}

predict.FeatureWeights <- function(model, x, predictionType, sel, exrep=NULL,
                                   raw, verbose, ...)
{
    ## in the current implementation multiclass prediction is not
    ## relevant for this function but is always performed natively in the SVM

    ## $$$ TODO  For large number of samples predict individually or in blocks
    if (is.null(x) && is.null(exrep))
        stop("missing data for prediction\n")

    if (predictionType == "probabilities" && !model@svmInfo@probModel)
        stop("missing probability model in model\n")

    if (length(kernelParameters(model@svmInfo@selKernel)$distWeight) > 0)
    {
        pred <- predict.PositionDependent(model, x, predictionType, sel,
                                          verbose, ...)
        return(pred)
    }

    ## calculate predictions with precalculated feature weights
    if (verbose)
    {
        classifierType <- kebabsInfo@classifierMap[model@svmInfo@selSVM,
                                                   model@svmInfo@selPackage]

        verbM(paste("predict - kbsvm with feature weights:"),
              classifierType, list(...))
    }

    ## get relevant features for prediction
    if (model@svmInfo@reqFeatureType == "quadratic")
    {
        relevantFeatures <- colnames(model@featureWeights)

        ## determine single features from pairs
        relevantFeatures <- strsplit(relevantFeatures, split="_")
        relevantFeatures <- sort(unique(unlist(relevantFeatures)))
    }
    else
    {
        if (is.list(model@featureWeights))
        {
            relevantFeatures <-
                    sort(unique(unlist(lapply(model@featureWeights,
                                           function(x) colnames(x)))))
        }
        else
            relevantFeatures <- colnames(model@featureWeights)
    }

    ## subset ER to / generate ER for relevant features
    if (!is.null(exrep))
    {
        if (exrep@quadratic == TRUE)
        {
            stop("please call predict via feature weights with\n",
                 "        linear explicit representation\n")
        }

        if (!(model@svmInfo@reqFeatureType == "quadratic" &&
              model@svmInfo@selPackage == "LiblineaR"))
        {
            relevantFeatures <- intersect(relevantFeatures, colnames(exrep))

            if(!isTRUE(all.equal(colnames(exrep), relevantFeatures)))
                exrep <- exrep[, relevantFeatures]
        }

        if (nrow(exrep) > 1000)
            blockSize <- 1000
        else
            blockSize <- nrow(exrep)

        numSamples <- nrow(exrep)
    }
    else
    {
        if (length(sel) > 1000)
            blockSize <- 1000
        else
            blockSize <- length(sel)

        numSamples <- length(sel)
    }

#    for (i in 1:ceiling(numSamples / blockSize))
#    {
#        range <- as.integer(((i - 1) * blockSize + 1):(i * blockSize))
#
#        if (range[blockSize] > numSamples)
#            range <- as.integer(((i - 1) * blockSize + 1):numSamples)
#
#    }

    ## for prediction only linear exrep is used
    ## quadratic will be generated later if quadratic kernel
    if (is.null(exrep))
    {
        if (class(model@svmInfo@selKernel) == "MotifKernel")
        {
            ## motif kernel does not support feature subsets
            exrep <- getExRep(x=x, selx=sel, kernel=model@svmInfo@selKernel,
                              sparse=model@ctlInfo@sparse)

            relevantFeatures <- intersect(relevantFeatures, colnames(exrep))

            if (!isTRUE(all.equal(colnames(exrep), relevantFeatures)))
                exrep <- exrep[, relevantFeatures]
        }
        else
        {
            exrep <- getExRep(x=x, selx=sel, kernel=model@svmInfo@selKernel,
                              sparse=model@ctlInfo@sparse,
                              features=relevantFeatures, zeroFeatures=TRUE)
        }
    }

    if (model@svmInfo@reqFeatureType == "quadratic")
    {
            erq <- getExRepQuadratic(exrep)

            relevantFeatures <- intersect(colnames(erq),
                                          colnames(model@featureWeights))

            if (!isTRUE(all.equal(colnames(erq), relevantFeatures)))
                erq <- erq[, relevantFeatures]

            if (model@ctlInfo@classification && (predictionType == "response"))
            {
                pred <- 1.5 - 0.5 * sign(model@b + erq %*%
                                    model@featureWeights[1,relevantFeatures])
            }
            else if (!model@ctlInfo@classification ||
                     (model@ctlInfo@classification &&
                      (predictionType == "decision")))
            {
                pred <- model@b + erq %*%
                        model@featureWeights[1,relevantFeatures]
            }
            else
                pred <- NULL
    }
    else
    {
        ##    the subsetting of feature weights is only local to this
        ##  function and not returned to the calling level
        if (!isTRUE(all.equal(colnames(model@featureWeights),
                              relevantFeatures)))
            model@featureWeights <- model@featureWeights[,relevantFeatures,
                                                         drop=FALSE]

        ## could be explicit with linear kernel or featureWeights
        ## for training with kernel matrix
        if (ncol(exrep) == 0)
            pred <- rep(model@b, nrow(exrep))
        else
        {
            # pred <- model@b + exrep %*% model@featureWeights
            ## work around - missing matrix multipl of dgRMatrix
            pred <- model@b + (as(exrep, "matrix") %*%
                               t(model@featureWeights))
        }

        if (model@ctlInfo@classification && (predictionType == "response"))
            pred <- 1.5 - 0.5 * sign(pred)
    }

    if (is.matrix(pred) && ncol(pred) == 1)
        pred <- as.numeric(pred)

    if (model@ctlInfo@classification == TRUE)
    {
        if (predictionType == "response")
        {
            pred <- model@classNames[pred]

            if (is.integer(model@levels))
                pred  <- as.integer(pred)
            else if (is.factor(model@levels))
                pred <- factor(pred, levels=model@levels)
        }
        else if (predictionType == "probabilities")
            pred <- getProbability(pred, model@probA, model@probB)
    }

    return(pred)
}

predict.NonFeatureWeights <- function(model, x, predictionType, sel, raw,
                                      verbose, ...)
{
    if (is(x, "KernelMatrix"))
    {
        if (length(sel) > 0)
            x <- x[sel,]

        pred <- predictSVM(x=x, model=model, predictionType=predictionType,
                           verbose=verbose, ...)

        if (raw != TRUE)
            pred <- unifyPredictionOutput(model=model, pred=pred,
                                          predictionType=predictionType)

        return(pred)
    }

    if (inherits(x, "ExplicitRepresentation"))
    {
        if (model@svmInfo@selPackage == "LiblineaR" &&
            model@svmInfo@reqFeatureType == "quadratic" &&
            x@quadratic != TRUE)
            x <- getExRepQuadratic(x)

        pred <- predictSVM(x=x, model=model, predictionType=predictionType,
                           verbose=verbose, ...)
    }
    else
    {
        if (model@ctlInfo@selMethod == "explicitRep")
        {
            if (class(model@svmInfo@selKernel) == "MotifKernel")
            {
                x <- getExRep(x=x, kernel=model@svmInfo@selKernel,
                              sparse=model@ctlInfo@sparse, selx=sel)

                if (model@svmInfo@reqFeatureType == "quadratic" &&
                    model@svmInfo@selPackage == "LiblineaR")
                    x <- getExRepQuadratic(x, zeroFeatures=TRUE)

                ## subset to features used in training
                if (!is(x, "ExplicitRepresentationSparse") &&
                    !isTRUE(all.equal(colnames(x), model@trainingFeatures)))
                    x <- x[,model@trainingFeatures]
            }
            else
            {
                if (model@svmInfo@reqFeatureType == "quadratic" &&
                    model@svmInfo@selPackage == "LiblineaR")
                {
                    linearFeatures <-
                        getSingleFeaturesFromQuadratic(model@trainingFeatures)

                    x <- getExRep(x=x, kernel=model@svmInfo@selKernel,
                                  sparse=model@ctlInfo@sparse, selx=sel,
                                  features=linearFeatures,
                                  zeroFeatures=TRUE)

                    x <- getExRepQuadratic(x)

                    if(!isTRUE(all.equal(colnames(x), model@trainingFeatures)))
                        x <- x[,model@trainingFeatures]
                }
                else
                {
                    x <- getExRep(x=x, kernel=model@svmInfo@selKernel,
                                  sparse=model@ctlInfo@sparse, selx=sel,
                                  features=model@trainingFeatures,
                                  zeroFeatures=TRUE)
                }
            }
        }
        else
        {
            if (length(model@SV) < 2)
                stop("support vectors missing in model")

            x <- getKernelMatrix(kernel=model@svmInfo@selKernel, x=x,
                                 y=model@SV, selx=sel)
        }

        pred <- predictSVM(x=x, model=model, predictionType=predictionType,
                           verbose=verbose, ...)
    }

    if (raw != TRUE)
    {
        pred <- unifyPredictionOutput(model=model, pred=pred,
                                      predictionType=predictionType)
    }

    return(pred)
}

predict.KBModel <- function(object, x, predictionType="response", sel=NULL,
                            raw=FALSE, native=FALSE, predProfiles=FALSE,
                            verbose = getOption("verbose"), ...)
{

   ## $$$ TODO
#    check compatibility of predict parameters with chosen SVM
#    addArgs <- match.call(expand.dots=FALSE)$...
#    svmInfo <- getSVMParams(object, ...)

    if (!is.null(predProfiles) && !is.logical(predProfiles))
        stop("'predProfiles' must be TRUE or FALSE\n")

    if(!(predictionType %in% c("response", "decision", "probabilities")))
        stop("wrong value for 'predictionType'\n")

    if (is.null(object@svmModel))
        stop("missing svm model for prediction\n")

    if (missing(x))
    {
        pred <- predictSVM(model=object, predictionType=predictionType,
                           verbose=verbose, ...)
        pred <- unifyPredictionOutput(model=object, pred=pred,
                                      predictionType=predictionType)
        return(pred)
    }

    x <- switch(class(x),
                "DNAString"    = DNAStringSet(x),
                "RNAString"    = RNAStringSet(x),
                "AAString"     = AAStringSet(x),
                x)

    if (!inherits(x, "XStringSet") && !inherits(x, "BioVector") &&
        !inherits(x, "ExplicitRepresentation") &&
        !(class(x) %in% c("kernelMatrix", "KernelMatrix")))
        stop("wrong class of x\n")

    if (is(x, "kernelMatrix"))
        x <- as.KernelMatrix(x)

    if (inherits(x, "XStringSet") || inherits(x, "BioVector"))
    {
        if (length(sel) > 0)
        {
            if (!is.numeric(sel) || min(abs(sel) < 1) ||
                max(abs(sel) > length(x)))
                stop("'sel' must contain indices into 'x'\n")

                if (all(sel < 0))
                    sel <- c(1:length(x))[-sel]
        }
        else
            sel <- 1L:length(x)
    }
    else
    {
        if (length(sel) > 0)
            stop("parameter 'sel' only allowed for sequences\n")
    }

    ## last part to allow prediction for regression via feature weigths
    if (length(object@featureWeights) > 0 && !(is(x, "KernelMatrix")) &&
        (!(object@numClasses > 2) || !object@ctlInfo@classification) &&
        !native)
    {
        ## predict allows dispatching only on single parameter
        if (inherits(x, "ExplicitRepresentation"))
        {
            pred <- predict.FeatureWeights(model=object, x=NULL,
                                           predictionType=predictionType,
                                           sel=sel, exrep=x, raw=raw,
                                           verbose=verbose, ...)
        }
        else
        {
            pred <- predict.FeatureWeights(model=object, x=x,
                                           predictionType=predictionType,
                                           sel=sel, raw=raw,
                                           verbose=verbose, ...)
        }
    }
    else
    {
        pred <- predict.NonFeatureWeights(model=object, x=x,
                                          predictionType=predictionType,
                                          sel=sel, raw=raw,
                                          verbose=verbose, ...)
    }

    if (!predProfiles)
        return(pred)
    else
    {
        if (!(inherits(x, "XStringSet") ||
              inherits(x, "BioVector")))
        {
            warning("Prediction profiles can only be computed for\n",
                    "        sequence data\n")
            return(pred)
        }

        if (is.null(object@featureWeights))
        {
            stop("Prediction profiles cannot be computed for models without\n",
                 "        feature weights\n")
        }

        predictionProfiles <- getPredictionProfile(x, object@svmInfo@selKernel,
                                                  object@featureWeights,
                                                  object@b, sel=sel)

        return(list(predictions=pred, predictionProfiles=predictionProfiles))
    }
}

#' @rdname predict-methods
#' @aliases predict predict.kbsvm predict.KBModel
#'
#' @title KeBABS Prediction Methods
#'
#' @description predict response values for new biological sequences from a
#' model trained with \code{kbsvm}
#'
#' @param object model object of class \code{\linkS4class{KBModel}}
#' created by \code{\link{kbsvm}}.
#'
#' @param x multiple biological sequences in the form of a
#' \code{\linkS4class{DNAStringSet}}, \code{\linkS4class{RNAStringSet}},
#' \code{\linkS4class{AAStringSet}} (or as \code{\linkS4class{BioVector}}). Also
#' a precomputed kernel matrix (see \code{\link{getKernelMatrix}} or a
#' precomputed explicit representation (see \code{\link{getExRep}} can be used
#' instead. The same type of input that was used for training the model should
#' also be used for prediction. If the parameter \code{x} is missing the
#' response is computed for the sequences used for SVM training.
#'
#' @param predictionType one character string of either "response",
#' "probabilities" or "decision" which indicates the type of data returned by
#' prediction: predicted response, class probabilities or decision values. Class
#' probabilities can only be computed if a probability model was generated
#' during the training (for details see parameter \code{probModel} in
#' \link{kbsvm}). Default=\code{"response"}
#'
#' @param sel subset of indices into \code{x}. When this parameter is present
#' the training is performed for the specified subset of samples only.
#' Default=\code{integer(0)}
#'
#' @param raw when setting this boolean parameter to TRUE the prediction result
#' is returned in raw form, i.e. in the SVM specific format. Default=FALSE
#'
#' @param native when setting this boolean parameter to TRUE the prediction is
#' not preformed via feature weights in the KeBABS model but native in the SVM.
#' Default=FALSE
#'
#' @param predProfiles when this boolean parameter is set to TRUE the
#' prediction profiles are computed for the samples passed to \code{predict}.
#' Default=FALSE
#'
#' @param verbose boolean value that indicates whether KeBABS should print
#' additional messages showing the internal processing logic in a verbose
#' manner. The default value depends on the R session verbosity option.
#' Default=getOption("verbose")
#'
#' @param ... additional parameters which are passed to SVM prediction
#' transparently.
#'
#' @details
#'
#' Prediction for KeBABS models\cr
#'
#' For the samples passed to the \code{predict} method the response (which
#' corresponds to the predicted label in case of classification or the predicted
#' target value in case of regression), the decision value (which is the value
#' of decision function separating the classes in classification) or the
#' class probability (probability for class membership in classification) is
#' computed for the given model of class \code{\linkS4class{KBModel}}. (see
#' also parameter \code{predictionType}). For sequence data this includes the
#' generation of an explicit representation or kernel matrix dependent on the
#' processing variant that was chosen for the training of the model. When
#' feature weights  were computed during training (see parameter
#' \code{featureWeights} in \code{\link{kbsvm}}) the response is computed
#' entirely in KeBABS via the feature weights in the model object. The
#' prediction performance can be evaluated with the function
#' \code{\link{evaluatePrediction}}.\cr\cr
#' If feature weights are not available in the model then native prediction
#' is performed via the SVM which was used for training. The parameter
#' \code{native} enforces native prediction even when feature weights are
#' available. Instead of sequence data also a precomputed kernel matrix or a
#' precomputed explicit representation can be passed to \code{predict}.
#' Prediction via feature weights is not supported for kernel variants which
#' do not support the generation of an explicit representation, e.g. the
#' position dependent kernel variants.\cr\cr
#'
#' Prediction with precomputed kernel matrix
#'
#' When training was performed with a precomputed kernel matrix also in
#' prediction a precomputed kernel matrix must be passed to the \code{predict}
#' method. In contrast to the quadratic and symmetric kernel matrix used
#' in training the kernel matrix for prediction is rectangular and contains
#' the similarities of test samples (rows) against support vectors (columns).
#' support vector indices can be read from the model with the accessor SVindex.
#' Please not that these indices refer to the sample subset used in training.
#' An example for training and prediction via precomputed kernel matrix is
#' shown below.
#'
#' Generation of prediction profiles
#'
#' The parameter \code{predProfiles} controls whether prediction profiles
#' (for details see \code{\link{getPredictionProfile}}) are generated during
#' the prediction process for all predicted samples. They show the contribution
#' of the individual sequence positions to the response value. For a subset of
#' sequences prediction profiles can also be computed independent from
#' predicition via the function \code{\link{getPredictionProfile}}.
#'
#'
#' @return
#' predict.kbsvm: upon successful completion, dependent on the parameter
#' \code{predictionType} the function returns either response values,
#' decision values or probability values for class membership. When prediction
#' profiles are also generated a list containing predictions and prediction
#' profiles is passed back to the user.
#'
#' @seealso \code{\linkS4class{KBModel}}, \code{\link{evaluatePrediction}},
#' \code{\link{kbsvm}}, \code{\link{getPredictionProfile}},
#' \code{\linkS4class{PredictionProfile}}
#'
#'
#' @examples
#'
#' ## load transcription factor binding site data
#' data(TFBS)
#' enhancerFB
#' ## select 70% of the samples for training and the rest for test
#' train <- sample(1:length(enhancerFB), length(enhancerFB) * 0.7)
#' test <- c(1:length(enhancerFB))[-train]
#' ## create the kernel object for gappy pair kernel with normalization
#' gappy <- gappyPairKernel(k=1, m=1)
#' ## show details of kernel object
#' gappy
#'
#' ## run training with explicit representation
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=gappy,
#'                pkg="LiblineaR", svm="C-svc", cost=10)
#'
#' ## show feature weights in KeBABS model
#' featureWeights(model)[1:8]
#'
#' ## predict the test sequences
#' pred <- predict(model, enhancerFB[test])
#' evaluatePrediction(pred, yFB[test], allLabels=unique(yFB))
#' pred[1:10]
#'
#' ## output decision values instead
#' pred <- predict(model, enhancerFB[test], predictionType="decision")
#' pred[1:10]
#'
#' \dontrun{
#' ## example for training and prediction via precomputed kernel matrix
#'
#' ## compute quadratic kernel matrix of training samples
#' kmtrain <- getKernelMatrix(gappy, x=enhancerFB, selx=train)
#'
#' ## train model with kernel matrix
#' model <- kbsvm(x=kmtrain, y=yFB[train], kernel=gappy,
#'                pkg="kernlab", svm="C-svc", cost=10)
#'
#' ## compute rectangular kernel matrix of test samples versus
#' ## support vectors
#' kmtest <- getKernelMatrix(gappy, x=enhancerFB, y=enhancerFB,
#'                           selx=test, sely=train[SVindex(model)])
#'
#' ## predict with kernel matrix
#' pred <- predict(model, kmtest)
#' evaluatePrediction(pred, yFB[test], allLabels=unique(yFB))
#'
#' ## example for probability model generation during training
#'
#' ## compute probability model via Platt scaling during training
#' ## and predict class membership probabilities
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=gappy,
#'                pkg="e1071", svm="C-svc", cost=10, probModel=TRUE)
#'
#' ## show parameters of the fitted probability model which are the parameters
#' ## probA and probB for the fitted sigmoid function in case of classification
#' ## and the value sigma of the fitted Laplacian in case of a regression
#' probabilityModel(model)
#'
#' ## predict class probabilities
#' prob <- predict(model, enhancerFB[test], predictionType="probabilities")
#' prob[1:10]
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}
#' @keywords predict
#' @keywords prediction
#' @keywords feature weights
#' @keywords prediction profile
#' @keywords methods
#'

#' @rdname predict-methods
#' @export


setMethod("predict", signature(object = "KBModel"), predict.KBModel)
