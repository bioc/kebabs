##345678901234567890123456789012345678901234567890123456789012345678901234567890
#' @rdname SVMAccess
#' @name predictSVM
#' @aliases
#' predictSVM,KernelMatrix-method

predictSVM.KernelMatrix <- function(x, model, predictionType, verbose, ...)
{
    addArgs <- list(...)

    classifierType <- kebabsInfo@classifierMap[model@svmInfo@selSVM,
                                               model@svmInfo@selPackage]

    if (model@svmInfo@selPackage == "kernlab")
    {
        if (is(model@svmModel, "ksvm"))
        {
            if (verbose)
            {
                verbM(paste("predict - ksvm with kernelMatrix:"),
                      classifierType, addArgs)
            }

            if (!is(x, "kernelMatrix"))
                x <- as.kernelMatrix(x)

            pred <- predict(object=model@svmModel, newdata=x,
                            type=predictionType, ...)

            if (is.matrix(pred) && ncol(pred) == 1)
                pred <- as.numeric(pred)

            return(pred)
        }
        else
            stop("SVM model is not ksvm\n")
    }
    else if (model@svmInfo@selPackage == "e1071")
    {
        if (is(model@svmModel, "svm"))
        {
            if (verbose)
            {
                verbM(paste("predict - svmd with kernelMatrix:"),
                classifierType, addArgs)
            }

            pred <- predict.svmd(object=model@svmModel, newdata=x,
                        decision.values=(predictionType=="decision"), ...)

            if (predictionType=="decision")
                pred <- attr(pred, "decision.values")
            
            return(pred)
        }
        else
            stop("SVM model is not svm\n")
    }
    else if (model@svmInfo@selPackage == "LiblineaR")
        stop("kernel matrix via LiblineaR is not supported\n")
    else
        stop(paste("unsupported package:", model@svmInfo@selPackage, "\n"))
}

setMethod("predictSVM", signature(x="KernelMatrix"),
          predictSVM.KernelMatrix)

predictSVM.missing <- function(x, model, predictionType, verbose, ...)
{
    addArgs <- list(...)

    classifierType <- kebabsInfo@classifierMap[model@svmInfo@selSVM,
                                               model@svmInfo@selPackage]

    if (model@svmInfo@selPackage == "kernlab")
    {
        if (is(model@svmModel, "ksvm"))
        {
            if (verbose)
            {
                verbM(paste("predict of training samples - ksvm :"),
                      classifierType, addArgs)
            }

            pred <- predict(object=model@svmModel, type=predictionType, ...)

            return(pred)
        }
        else
            stop("SVM model is not ksvm\n")
    }
    else if (model@svmInfo@selPackage == "e1071")
    {
        if (is(model@svmModel, "svm"))
        {
            if (verbose)
            {
                verbM(paste("predict of training samples - svm :"),
                      classifierType, addArgs)
            }

            ## if precomputed kernel matrix
            if (model@svmModel$kernel ==4)
            {
                pred <- predict.svmd(object=model@svmModel,
                            decision.values=(predictionType=="decision"), ...)
            }
            else
            {
                pred <- predict(object=model@svmModel,
                            decision.values=(predictionType=="decision"), ...)
            }

            if (predictionType=="decision")
               pred <- attr(pred, "decision.values")

            return(pred)
        }
        else
            stop("SVM model is not svm\n")
    }
    else
        stop(paste("unsupported package:", model@svmInfo@selPackage, "\n"))
}

#' @rdname SVMAccess
#' @aliases
#' predictSVM,missing-method
#'

setMethod("predictSVM", signature(x="missing"),
          predictSVM.missing)

predictSVM.ExRep <- function(x, model, predictionType, verbose, ...)
{
    addArgs <- list(...)

    classifierType <- kebabsInfo@classifierMap[model@svmInfo@selSVM,
                                               model@svmInfo@selPackage]

    erType <- " linear "
    sparse <- " (dense) "

    if (is(x, "ExplicitRepresentationSparse"))
        sparse=" (sparse) "

    if (model@svmInfo@selPackage == "kernlab")
    {
        ## $$$ Remove when kernlab is supporting dgRMatrix
        if (is(x, "ExplicitRepresentationSparse"))
            stop("ksvm currently does not support sparse data\n")

        if (is(model@svmModel, "ksvm"))
        {
            if (verbose)
            {
                verbM(paste("predict ksvm with", erType,
                            "explicit representation", sparse, sep=""),
                      classifierType, addArgs)
            }

            pred <- predict(object=model@svmModel, newdata=x,
                            type=predictionType, ...)

            return(pred)
        }
        else
            stop("SVM model is not ksvm\n")
    }
    else if (model@svmInfo@selPackage == "e1071")
    {
        if (is(model@svmModel, "svm"))
        {
            if (verbose)
            {
                verbM(paste("predict svm with", erType,
                            "explicit representation", sparse, sep=""),
                      classifierType, addArgs)
            }

            ## $$$ Remove conversion when e1071 is supporting dgRMatrix
            if (is(x, "ExplicitRepresentationSparse"))
                x <- as(x, "matrix.csr")

            pred <- predict(object=model@svmModel, newdata=x,
                            decision.values=(predictionType=="decision"), ...)

            if (predictionType=="decision")
                pred <- attr(pred, "decision.values")

            return(pred)
        }
        else
            stop("SVM model is not svm\n")
    }
    else if (model@svmInfo@selPackage == "LiblineaR")
    {
        if (class(model@svmModel) == "LiblineaR")
        {
            if (verbose)
            {
                if (x@quadratic)
                    erType <- " quadratic "

                verbM(paste("predict LiblineaR with", erType,
                            "explicit representation", sparse, sep=""),
                      classifierType, addArgs)
            }

            ## $$$ Remove conversion when LiblineaR is supporting dgRMatrix
            if (is(x, "ExplicitRepresentationSparse"))
                x <- as(x, "matrix.csr")

            pred <- predict(object=model@svmModel, newx=x,
                            decisionValues=(predictionType=="decision"), ...)

            if (predictionType=="decision")
                pred <- pred$decisionValues
            else
                pred <- pred$predictions

            return(pred)
        }
        else
            stop("SVM model is not LiblineaR\n")
    }
    else
        stop(paste("unsupported package:", model@svmInfo@selPackage, "\n"))
}

#' @rdname SVMAccess
#' @aliases
#' predictSVM,ExpicitRepresentation-method
#'

setMethod("predictSVM", signature(x="ExplicitRepresentation"),
          predictSVM.ExRep)
