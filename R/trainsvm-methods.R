##345678901234567890123456789012345678901234567890123456789012345678901234567890
trainSVM.KernelMatrix <- function(x, y=NULL, svmInfo, verbose, ...)
{
    addArgs <- list(...)

    classifierType <- kebabsInfo@classifierMap[svmInfo@selSVM,
                                               svmInfo@selPackage]

    if (nchar(classifierType) < 1)
    {
        stop(paste("Classifier", svmInfo@selSVM, "is not available\n"),
             paste("       in package", svmInfo@selPackage, "\n"))
    }

    if (svmInfo@selPackage == "kernlab")
    {
#        library("kernlab")

        if (verbose)
        {
            if (svmInfo@reqFeatureType == "quadratic")
            {
                verbM(paste("train ksvm with squared kernelMatrix:"),
                            classifierType, addArgs)
            }
            else
            {
                verbM(paste("train ksvm with kernelMatrix:"),
                            classifierType, addArgs)
            }
        }

        if (svmInfo@reqFeatureType == "quadratic")
            x <- x^2

        if (!is(x, "kernelMatrix"))
            x <- as.kernelMatrix(x)

        ## set scaling always to FALSE - needed also for kernel matrix
        scaling <- FALSE

        return(ksvm(x=x, y=y, type=classifierType,
                    prob.model=svmInfo@probModel, scaled=scaling, ...))
    }
    else if (svmInfo@selPackage == "e1071")
        stop("kernel matrix via e1071 is currently not supported\n")
    else if (svmInfo@selPackage == "LiblineaR")
        stop("kernel matrix via LiblineaR is not supported\n")
    else
        stop(paste("unsupported package:", svmInfo@selPackage, "\n"))
}

#' @rdname SVMAccess
#' @name trainSVM
#' @aliases
#' trainSVM,KernelMatrix-method
#' @title SVM Access for Training and Prediction
#'
#' @description
#' Functions for SVM access (used only for testing purpose)
#
#' @param x kernel matrix or explicit representation
#' @param y label vector
#' @param svmInfo SVM related info
#' @param model KeBABS model
#' @param predictionType type of prediction
#' @param verbose controlling verbosity
#' @param ... additional arguments to be passed to the selected SVM
#' @details
#' These methods are exported only for test purpose and are not meant
#' to be generally used.
#' @return \code{trainSVM:} returns the SVM specific model\cr\cr
#' \code{predictSVM:} returns the prediction in native format
#' @examples
#' ## this function is exported only for testing purpose
#' ## use function kbsvm instead for examples see help page of kbsvm
#' data(TFBS)

setMethod("trainSVM", signature(x="KernelMatrix"),
          trainSVM.KernelMatrix)

##setMethod("trainSVM", signature(x="kernelMatrix"),
##          trainSVM.KernelMatrix)


trainSVM.explicitRep <- function(x, y=NULL, svmInfo, verbose, ...)
{
    addArgs <- list(...)

    classifierType <- kebabsInfo@classifierMap[svmInfo@selSVM,
                                               svmInfo@selPackage]

    if (nchar(classifierType) < 1)
    {
        stop(paste("Classifier", svmInfo@selSVM, "is not available\n"),
             paste("       in package", svmInfo@selPackage, "\n"))
    }

    ## set scaling always to FALSE to not up-/downweight single features
    scaling <- FALSE

    erType <- " linear "
    sparse <- " (dense) "

    if (is(x, "ExplicitRepresentationSparse"))
        sparse=" (sparse) "

    ## only linear and quadratic kernel relevant
    ## kernels are set in the SVM parameter routine
    if (svmInfo@selPackage == "kernlab")
    {
#        library("kernlab")

        if (is(x, "ExplicitRepresentationSparse"))
            stop("ksvm currently does not support sparse data\n")

        if (svmInfo@explicitKernel == "linear")
        {
            if (verbose)
            {
                verbM(paste("train ksvm with", erType,
                            "explicit representation", sparse,
                            "\nand vanilladot:", sep=""),
                      classifierType, addArgs)
            }

            return(ksvm(x=x, y=y, type=classifierType,
                        kernel=vanilladot(), kpar=list(), scaled=scaling,
                        prob.model=svmInfo@probModel, ...))
        }
        else if (svmInfo@explicitKernel == "quadratic")
        {
            if (verbose)
            {
                verbM(paste("train ksvm with", erType,
                            "explicit representation", sparse,
                            "\nand polydot(degree=2):", sep=""),
                      classifierType, addArgs)
            }

            return(ksvm(x=x, y=y, type=classifierType,
                        kernel=polydot(degree=2, offset=0), scaled=scaling,
                        prob.model=svmInfo@probModel, ...))
        }
        else
            stop("Wrong type explicit kernel\n")
    }
    else if (svmInfo@selPackage == "e1071")
    {
        if (svmInfo@explicitKernel == "linear")
        {
            if (verbose)
            {
                verbM(paste("train svm with", erType,
                            "explicit representation", sparse,
                            "\nand linear kernel:", sep=""),
                      classifierType, addArgs)
            }

            ## $$$ Remove conversion when e1071 is supporting dgRMatrix
            if (is(x, "ExplicitRepresentationSparse"))
            {
                library(SparseM, warn.conflicts=FALSE)
                x <- as(x, "matrix.csr")
            }

            return(svm(x=x, y=y, type=classifierType,
                       kernel="linear", scale=scaling,
                       probability=svmInfo@probModel, ...))
        }
        else if (svmInfo@explicitKernel == "quadratic")
        {
            if (verbose)
            {
                verbM(paste("train svm with", erType,
                            "explicit representation", sparse,
                            "\nand polynomial kernel(degree=2):", sep=""),
                      classifierType, addArgs)
            }

            ## $$$ Remove conversion when e1071 is supporting dgRMatrix
            if (is(x, "ExplicitRepresentationSparse"))
            {
                library(SparseM, warn.conflicts=FALSE)
                x <- as(x, "matrix.csr")
            }

            return(svm(x=x, y=y, type=classifierType,
                       kernel="polynomial", gamma=1, coef0=0, degree=2,
                       scale=scaling, probability=svmInfo@probModel,  ...))
        }
        else
            stop("Wrong type explicit kernel\n")
    }
    else if (svmInfo@selPackage == "LiblineaR")
    {
        if (verbose)
        {
            if (x@quadratic)
                erType <- " quadratic "

            verbM(paste("train LiblineaR with", erType,
                        "explicit representation", sparse, sep=""),
                  classifierType, addArgs)
        }

        ## $$$ Remove conversion when LiblineaR is supporting dgRMatrix
        if (is(x, "ExplicitRepresentationSparse"))
        {
            featureNames <- colnames(x)
            library(SparseM, warn.conflicts=FALSE)
            x <- as(x, "matrix.csr")
        }
        
        ## check version because of interface change with 1.94-1
        liblinearVersion <- packageVersion("LiblineaR")
        if (liblinearVersion$major == 1 &&
            liblinearVersion$minor < 94)
        {
            model <- LiblineaR(data=x, labels=y,
                               type=as.integer(classifierType), ...)
        }
        else
        {
            model <- LiblineaR(data=x, target=y,
                               type=as.integer(classifierType), ...)
        }

        ## $$$ Remove name assignment when LiblineaR is supporting dgRMatrix
        ## matrix.csr does not support names =>
        ## set up names for feature weights in model for sparse ER
        if (is(x, "matrix.csr") && length(featureNames) > 0 &&
            is(model, "LiblineaR"))
        {
            if (model$Bias)
                colnames(model$W) <- c(featureNames, "Bias")
            else
                colnames(model$W) <- featureNames
        }

        return(model)
     }
    else
        stop(paste("unsupported package:", svmInfo@selPackage, "\n"))
}

#' @rdname SVMAccess
#' @aliases
#' trainSVM,ExplicitRepresentation-method
#'

setMethod("trainSVM", signature(x="ExplicitRepresentation"),
          trainSVM.explicitRep)

verbM <- function(text, classifierType, addArgs)
{
    message(text)
    message("  svm type:", classifierType)

    if (length(addArgs) > 0)
    {
        message("  svm parameters:")
        message("  ", listToString(addArgs))
    }
}
