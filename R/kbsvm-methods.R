##345678901234567890123456789012345678901234567890123456789012345678901234567890
kbsvm.seqs <- function(x, y, kernel=NULL, pkg="auto", svm="C-svc",
            explicit="auto", explicitType="auto", featureType="linear",
            featureWeights="auto", weightLimit=.Machine$double.eps,
            classWeights = numeric(0), cross=0, noCross=1, groupBy=NULL,
            nestedCross=0, noNestedCross=1, perfParameters=character(0),
            perfObjective="ACC", probModel=FALSE, sel=integer(0),
            features=NULL, showProgress=FALSE, showCVTimes=FALSE,
            runtimeWarning=TRUE, verbose = getOption("verbose"), ...)
{
    addArgs <- list(...)

    if (length(perfParameters) >= 1 && perfParameters[1] == "ALL")
        perfParameters <- c("ACC", "BACC", "MCC")

    ## check parameters and allocate model
    model <- checkKBSVMParams(x=x, y=y, kernel=kernel, svm=svm, pkg=pkg,
                    explicit=explicit, explicitType=explicitType,
                    featureType=featureType, featureWeights=featureWeights,
                    weightLimit=weightLimit, classWeights=classWeights,
                    probModel=probModel, sel=sel, features=features,
                    cross=cross, noCross=noCross, groupBy=groupBy,
                    nestedCross=nestedCross, noNestedCross=noNestedCross,
                    perfParameters=perfParameters, perfObjective=perfObjective,
                    runtimeWarning=runtimeWarning, addArgs=addArgs)

    model@call <- deparse(sys.call(-1))

    if (!is.null(model@modelSelResult))
    {
        return(performModelSelection(object=x, model=model, y=y, sel=sel,
                    explicit=explicit, features=features, groupBy=groupBy,
                    featureWeights=featureWeights, weightLimit=weightLimit,
                    perfParameters=perfParameters, perfObjective=perfObjective,
                    showProgress=showProgress, showCVTimes=showCVTimes,
                    addArgs=addArgs, verbose=verbose))
    }

    ## select SVM and method
    result <- selectSVMMethod(model, kebabsInfo@supportedPkgs, x)

    if (length(result) == 2)
    {
        model@svmInfo <- result[[1]]
        model@ctlInfo <- result[[2]]
    }
    else
        stop("wrong result from selectSVMMethod\n")

    ## generate runtime warnings
    model@ctlInfo <- generateRuntimeMessage(model, x)

    ## convert parameters to target SVM
    ## and remove unneeded parameters
    model@svmInfo <- convertSVMParameters(model)

    # change parameters if necessary
    ## $$$ TODO modify dotArgs
    ## finally do.call(<function name>, named pars, dotArgs)


    # part of the control info will go into the model
    # a part should be kept in an internal structure - maybe linked to the
    # model but not displayed with show

    if (model@ctlInfo@selMethod == "KernelMatrix")
    {
        if (length(features) > 0)
        {
            stop("feature subset not supported for processing via\n",
                 "        kernel matrix\n")
        }

        ## calculate kernel
        km <- getKernelMatrix(kernel=model@svmInfo@selKernel, x=x, selx=sel)

        if (cross != 0)
        {
            ## subset immediately
            if (length(sel) > 0)
                x <- x[sel]

            if (showCVTimes)
            {
                timeCV <- system.time(
                    cvModel <- performCrossValidation(object=km, x=x, y=y,
                                sel=integer(0), model=model, cross=cross,
                                noCross=noCross, groupBy=groupBy,
                                perfParameters=perfParameters, verbose=verbose))
                cat("\nCross Validation Time:\n")
                print(timeCV)
                cat("\n")
                return(cvModel)
            }
            else
            {
                return(performCrossValidation(object=km, x=x, y=y,
                                sel=integer(0), model=model, cross=cross,
                                noCross=noCross, groupBy=groupBy,
                                perfParameters=perfParameters, verbose=verbose))
            }
        }
        else
        {
            ## invoke selected svm
            allArgs <- model@svmInfo@selSVMPar
            allArgs[["x"]] <- km
            allArgs[["y"]] <- y
            allArgs[["svmInfo"]] <- model@svmInfo
            allArgs[["verbose"]] <- verbose
            model@svmModel <- do.call(trainSVM, allArgs)

            if (is.null(model@svmModel))
            {
                stop(paste(model@svmInfo@selPackage, "-",
                           model@svmInfo@selSVM, "did not return a model\n"))
            }

            if (model@ctlInfo@classification == TRUE)
            {
                model@numClasses <- getSVMSlotValue("numClasses", model)
                model@classNames <- getSVMSlotValue("classNames", model)

                if (model@numClasses > 2)
                    model@ctlInfo@multiclassType <- getMulticlassType(model)

                if (model@svmInfo@probModel == TRUE)
                {
                    model@probA <- getSVMSlotValue("probA", model)
                    model@probB <- getSVMSlotValue("probB", model)
                }
            }
            else
            {
                ## copy param for regression with probability model
                if (model@svmInfo@probModel == TRUE)
                    model@sigma <- getSVMSlotValue("sigma", model)
            }

            ind <- getSVMSlotValue("svIndex", model)

            if (length(ind) > 0)
            {
                if (length(sel) > 0)
                    indOrig <- sel[ind]
                else
                    indOrig <- ind

                model@SV <- x[indOrig]
                model@svIndex <- indOrig
                model@alphaIndex <- getSVMSlotValue("alphaIndex", model)
            }

            if (featureWeights != "no")
            {
                if (length(ind) > 0)
                {
                    ## generate compact ex rep of SVs
                    if (supportsExplicitRep(model@svmInfo@selKernel))
                    {
                        exRepSV <- getExRep(x=model@SV,
                                            kernel=model@svmInfo@selKernel,
                                            sparse=model@ctlInfo@sparse)

                        if (ncol(exRepSV) == 0)
                            stop("explicit representaton is not available\n")
                    }
                    else
                        exRepSV <- NULL

                    model@featureWeights <-
                        getFeatureWeights(model=model,
                                          exrep=exRepSV,
                                          weightLimit=model@svmInfo@weightLimit)

                    model@b <- getSVMSlotValue("b", model)
                }

                if (length(model@featureWeights) < 1)
                    stop("no feature weights available\n")
            }

            return(model)
        }
    }
    else if (model@ctlInfo@selMethod == "explicitRep")
    {
        if (!supportsExplicitRep(kernel))
            stop("selected kernel does not support explicit representation\n")

        exRep <- getExRep(x=x, selx=sel, kernel=kernel,
                          sparse=model@ctlInfo@sparse,
                          features=features)

        if (ncol(exRep) == 0)
            stop("explicit representaton is not available\n")

        if (model@svmInfo@reqFeatureType == "quadratic" &&
            model@svmInfo@selPackage == "LiblineaR" && cross == 0)
        {
            ## convert to explicit representation for quadratic kernel
            ## as quadratic kernel is not available in LiblineaR
            ## performCrossValidation is called with linear exrep
            exRepLin <- exRep
            exRep <- getExRepQuadratic(exRep)
        }

        ## store feature names in model, needed for native prediction
        ## via ex rep in svm
        model@trainingFeatures <- colnames(exRep)

        ## remove orig kernel and set up ER kernel


        if (cross != 0)
        {
            ## subset immediately
            if (length(sel) > 0)
                x <- x[sel]

            if (showCVTimes)
            {
                timeCV <- system.time(
                    cvModel <- performCrossValidation(object=exRep, x=x, y=y,
                                sel=integer(0), model=model, cross=cross,
                                noCross=noCross, groupBy=groupBy,
                                perfParameters=perfParameters, verbose=verbose))
                cat("\nCross Validation Time:\n")
                print(timeCV)
                cat("\n")
                return(cvModel)
            }
            else
            {
                return(performCrossValidation(object=exRep, x=x, y=y,
                                sel=integer(0), model=model, cross=cross,
                                noCross=noCross, groupBy=groupBy,
                                perfParameters=perfParameters, verbose=verbose))
            }
        }
        else
        {
            ## invoke selected svm
            allArgs <- model@svmInfo@selSVMPar
            allArgs[["x"]] <- exRep
            allArgs[["y"]] <- y
            allArgs[["svmInfo"]] <- model@svmInfo
            allArgs[["verbose"]] <- verbose
            model@svmModel <- do.call(trainSVM, allArgs)

            if (is.null(model@svmModel))
            {
                stop(paste(model@svmInfo@selPackage, "-",
                           model@svmInfo@selSVM, "did not return a model\n"))
            }

            if (model@ctlInfo@classification == TRUE)
            {
                model@numClasses <- getSVMSlotValue("numClasses", model)
                model@classNames <- getSVMSlotValue("classNames", model)

                if (model@numClasses > 2)
                    model@ctlInfo@multiclassType <- getMulticlassType(model)

                if (model@svmInfo@probModel == TRUE)
                {
                    model@probA <- getSVMSlotValue("probA", model)
                    model@probB <- getSVMSlotValue("probB", model)
                }
            }
            else
            {
                ## copy param for regression with probability model
                if (model@svmInfo@probModel == TRUE)
                    model@sigma <- getSVMSlotValue("sigma", model)
            }

            ind <- getSVMSlotValue("svIndex", model)
            exRepSV <- NULL

            if (length(ind) > 0)
            {
                if (length(sel) > 0)
                    indOrig <- sel[ind]
                else
                    indOrig <- ind

                model@SV <- x[indOrig]
                model@svIndex <- indOrig
                model@alphaIndex <- getSVMSlotValue("alphaIndex", model)

                if (model@svmInfo@reqFeatureType == "quadratic" &&
                    model@svmInfo@selPackage == "LiblineaR")
                {
                    exRepSV <- exRepLin[ind,]
                }
                else
                    exRepSV <- exRep[ind,]
            }

            if (featureWeights != "no")
            {
                model@featureWeights <-
                    getFeatureWeights(model=model,
                                      exrep=exRepSV,
                                      features=NULL,
                                      weightLimit=model@svmInfo@weightLimit)

                model@b <- getSVMSlotValue("b", model)

                if (length(model@featureWeights) < 1)
                    stop("no feature weights available\n")
            }

            return(model)
        }
    }
}


#' @rdname kbsvm-methods

#' @title KeBABS Training Methods
#'
#' @description Train an SVM-model with a sequence kernel on biological
#' sequences
#'
#' @param x multiple biological sequences in the form of a
#' \code{\linkS4class{DNAStringSet}}, \code{\linkS4class{RNAStringSet}},
#' \code{\linkS4class{AAStringSet}} (or as \code{\link{BioVector}}). Also
#' a precomputed kernel matrix (see \code{\link{getKernelMatrix}} or a
#' precomputed explicit representation (see \code{\link{getExRep}} can be used
#' instead. If they were precomputed with a sequence kernel this kernel should
#' be specified in the parameter \code{kernel} in this case.
#'
#' @param sel subset of indices into \code{x}. When this parameter is present
#' the training is performed for the specified subset of samples only.
#' Default=\code{integer(0)}
#'
#' @param y response vector which contains one value for each sample in 'x'.
#' For classification tasks this can be either a character vector, a factor or
#' a numeric vector, for regression tasks it must be a numeric vector. For 
#' numeric labels in binary classification the positive class must have the 
#' larger value, for factor or character based labels the positive label must 
#' be at the first position when sorting the labels in descendent order 
#' according to the C locale. If the parameter sel is used to perform training 
#' with a sample subset the response vector must have the same length as 'sel'.
#'
#' @param kernel a sequence kernel object or a string kernel from package
#' \link[kernlab:stringdot]{kernlab}. In case of grid search or model selection
#' a list of sequence kernel objects can be passed to training.
#'
#' @param pkg name of package which contains the SVM implementation to be used
#' for training, e.g. \code{kernlab}, \code{e1071} or \code{LiblineaR}. For
#' gridSearch or model selection multiple packages can be passed as character
#' vector. (see also parameter \code{svm} below). Default="auto"
#'
#' @param svm name of the SVM used for the classification or regression task,
#' e.g. "C-svc". For gridSearch or model selection multiple SVMs can be passed
#' as character vector. For each entry in this character vector a corresponding
#' entry in the character vector for parameter \code{pkg} is required, if
#' multiple SVMs are used in one cross validation or model selection run.
#'
#' @param explicit this parameter controls whether training should be performed
#' with the kernel matrix (see \code{\link{getKernelMatrix}}) or explicit
#' representation (see \code{\link{getExRep}}). When the parameter is set to
#' "no" the kernel matrix is used, for "yes" the model is trained from the
#' explicit representation. When set to "auto" KeBABS automatically selects a
#' variant based on runtime heuristics. Default="auto"
#'
#' @param explicitType this parameter is only relevant when parameter
#' 'explicit' is different from "no". The values "sparse" and "dense"
#' indicate whether a sparse or dense explicit representation should
#' be used. When the parameter is set to "auto" KeBABS selects a variant.
#' Default="auto"
#'
#' @param featureType when the parameter is set to "linear" single features
#' areused in the analysis (with a linear kernel matrix or a linear kernel
#' applied to the linear explicit representation). When set to "quadratic"
#' the analysis is based on feature pairs. For an SVM from
#' \link[LiblineaR:LiblineaR]{LiblineaR} (which does not support kernels)
#' KeBABS generates a quadratic explicit representation. For the other SVMs
#' a polynomial kernel of degree 2 is used for learning via explicit
#' representation. In the case of learning via kernel matrix a quadratic
#' kernel matrix (quadratic here in the sense of linear kernel matrix
#' with each element taken to power 2) is generated. Default="linear"
#'
#' @param featureWeights with the values "no" and "yes" the user can control
#' whether feature weights are calulated as part of the training. When the
#' parameter is set to "auto" KeBABS selects a variant (see below).
#' Default="auto"
#'
#' @param weightLimit the feature weight limit is a single numeric value and
#' allows pruning of feature weights. All feature weights with an absolute
#' value below this limit are set to 0 and are not considered in the model and
#' for further predictions. This parameter is only relevant when featureWeights
#' are calculated in KeBABS during training.
#' Default=.Machine$double.eps
#'
#' @param classWeights a numeric named vector of weights for the different
#' classes, used for asymmetric class sizes. Each element of the vector must
#' have one of the class names but not all class names must be present.
#' Default=1
#
#' @param cross an integer value K > 0 indicates that k-fold cross validation
#' should be performed. A value -1 is used for Leave-One-Out (LOO) cross
#' validation. (see above) Default=0
#'
#' @param noCross an integer value larger than 0 is used to specify the number
#' of repetitions for cross validation. This parameter is only relevant if
#' 'cross' is different from 0. Default=1
#
#' @param groupBy allows a grouping of samples during cross validation. The
#' parameter is only relevant when 'cross' is larger than 1. It is an integer
#' vector or factor with the same length as the number of samples used for
#' training and specifies for each sample to which group it belongs. Samples
#' from the same group are never spread over more than one fold. (see
#' \code{crossValidation}). Grouped cross validation can also be used in
#' grid search for each grid point. Default=NULL
#
#' @param perfParameters a character vector with one or several values from
#' the set "ACC" , "BACC", "MCC" and "ALL". "ACC" stands for accuracy, "BACC"
#' for balanced accuracy, "MCC" for Matthews Correlation Coefficient and "ALL"
#' for all three. This parameter defines which performance parameters are
#' collected in cross validation, grid search and model selection for display
#' purpose. Default=NULL
#'
#' @param perfObjective a singe character string from the set "ACC", "BACC"
#' and "MCC" (see previous parameter). The parameter is only relevant in
#' grid search and model selection and defines which performance measure is
#' used to determine the best performing parameter set. Default="ACC"
#'
#' @param showProgress when setting this boolean parameter to TRUE the
#' progress of a cross validation is displayed. The parameter is only relevant
#' for cross validation. Default=FALSE
#'
#' @param showCVTimes when setting this boolean parameter to TRUE the runtimes
#' of the cross validation runs are shown after the cross validation is
#' finished. The parameter is only relevant for cross validation.
#' Default=FALSE
#'
#' @param nestedCross in integer value K > 0 indicates that a model selection
#' with nested cross validation should be performed with a k-fold outer cross
#' validation. The inner cross validation is defined with the 'cross'
#' parameter (see below), Default=0
#
#' @param noNestedCross an integer value larger than 0 is used to specify the
#' number of repetitions for the nested cross validation. This parameter is
#' only relevant if 'nestedCross' is larger than 0. Default=1
#
#' @param features feature subset of the specified kernel in the form of a
#' character vector. When a feature subset is passed to the function all other
#' features in the feature space are not considered for training (see below).
#' A feature subset can only be used when a single kernel object is specified
#' in the 'kernel' parameter. Default=NULL
#'
#' @param probModel when setting this boolean parameter to TRUE a probability
#' model is determined as part of the training (see below). Default=FALSE
#'
#' @param runtimeWarning when setting this boolean parameter to FALSE a
#' warning for long runtimes will not be shown in case of large feature
#' space dimension or large number of samples. Default=TRUE
#'
#' @param verbose boolean value that indicates whether KeBABS should print
#' additional messages showing the internal processing logic in a verbose
#' manner. The default value depends on the R session verbosity option.
#' Default=getOption("verbose")
#'
#' @param ... additional parameters which are passed to SVM training
#' transparently.
#'
#' @details
#'
#' Overview\cr
#'
#' The kernel-related functionality provided in this package is specifically
#' centered around biological sequences, i.e. DNA-, RNA- or AA-sequences (see
#' also \code{\link{DNAStringSet}}, \code{\link{RNAStringSet}} and
#' \code{\link{AAStringSet}}) and Support Vector Machine (SVM) based methods.
#' Apart from the implementation of the most relevant kernels for sequence
#' analysis (see \code{\link{spectrumKernel}}, \code{\link{mismatchKernel}},
#' \code{\link{gappyPairKernel}} and \code{\link{motifKernel}}) KeBABS also
#' provides a framework which allows easy interworking with existing SVM
#' implementations in other R packages. In the current implementation the SVMs
#' provided in the packages \code{\link[kernlab:ksvm]{kernlab}},
#' \code{\link[e1071:svm]{e1071}} and
#' \code{\link[LiblineaR:LiblineaR]{LiblineaR}} are in focus.\cr\cr
#' This framework can be considered like a "meta-SVM", which provides
#' a simple and unified user interface to these SVMs for classification (binary
#' and multiclass) and regression tasks. The user calls the "meta-SVM" in a
#' classical SVM-like manner by passing sequence data, a sequence kernel with
#' kernel parameters and the SVM which should be used for the learning task
#' togehter with SVM parameters. KeBABS internally generates the relevant
#' representations (see \code{\link{getKernelMatrix}} or \code{\link{getExRep}})
#' from the sequence data using the specified kernel, adapts parameters and
#' formats to the selected SVM and internally calls the actual SVM
#' implementation in the requested package. KeBABS unifies the
#' result returned from the invoked SVM and returns a unified data structure,
#' the KeBABS model, which also contains the SVM-specific model (see
#' \code{\link{svmModel}}.\cr\cr
#' The KeBABS model is used in prediction (see \code{\link{predict}}) to
#' predict the response for new sequence data. On user request the feature
#' weights are computed and stored in the Kebabs model during training (see
#' below). The feature weights are used for the generation of prediction
#' profiles (see \code{\link{getPredictionProfile}}) which show the importance
#' of sequence positions for a specfic learning task.\cr\cr
#'
#' Training of biological sequences with a sequence kernel\cr
#'
#' Training is performed via the method \code{kbsvm} for classification and
#' regression tasks. The user passes sequence data, the response vector, a
#' sequence kernel object and the requested SVM along with SVM parameters
#' to \code{kbsvm} and receives the training results in the form of a
#' KeBABS model object of class \code{\linkS4class{KBModel}}. The accessor
#' \code{svmModel} allows to retrieve the SVM specific model from the KeBABS
#' model object. However, for regular operation a detailed look into the SVM
#' specific model is usually not necessary.
#'
#' The standard data format for sequences in KeBABS are the
#' \code{XStringSet}-derived classes \code{\linkS4class{DNAStringSet}},
#' \code{\linkS4class{RNAStringSet}} and \code{\linkS4class{AAStringSet}}.
#' (When repeat regions are coded as lowercase characters and should be
#' excluded from the analysis the sequence data can be passed as
#' \code{\linkS4class{BioVector}} which also supports lowercase characters
#' instead of \code{\linkS4class{XStringSet}} format. Please note that the
#' classes derived from \code{\linkS4class{XStringSet}} are much more
#' powerful than the \code{\linkS4class{BioVector}} derived classes and
#' should be used in all cases where lowercase characters are not needed).
#'
#' Instead of sequences also a precomputed explicit representation or
#' a precomputed kernel matrix can be used for training. Examples for
#' training with kernel matrix and explicit representation can be found on
#' the help page for the prediction method \code{\link{predict}}.
#'
#' Apart from SVM training \code{kbsvm} can be also used for cross
#' validation (see \link{crossValidation} and parameters \code{cross} and
#' \code{noCross}), grid search for SVM- and kernel-parameter values (see
#' \link{gridSearch}) and model selection (see \link{modelSelection} and
#' parameters \code{nestedCross} and \code{noNestedCross}).\cr\cr
#'
#' Package and SVM selection\cr
#'
#' The user specifies the SVM implementation to be used for a learning task by
#' selecting the package with the \code{pkg} parameter and the SVM method in
#' the package with the \code{SVM} parameter. Currently the packages
#' code{\link[kernlab:ksvm]{kernlab}}, \code{\link[e1071:svm]{e1071}} and
#' \code{\link[LiblineaR:LiblineaR]{LiblineaR}} are supported. The names for
#' SVM methods vary from package to package and KeBABS provide following
#' unified names which can be selected across packages. The following table
#' shows the available SVM methods:\cr
#'
#' \tabular{ll}{
#'   SVM name \tab description\cr
#'   ----------------------- \tab -----------------------------------------
#' ---------\cr
#'   C-svc:      \tab C classification (with L2 regularization and L1 loss)\cr
#'   l2rl2l-svc: \tab classif. with L2 regularization and L2 loss (dual)\cr
#'   l2rl2lp-svc:\tab classif. with L2 regularization and L2 loss (primal)\cr
#'   l1rl2l-svc: \tab classification with L1 regularization and L2 loss\cr
#'   nu-svc:     \tab nu classification\cr
#'   C-bsvc:     \tab bound-constraint SVM classification\cr
#'   mc-natC:    \tab Crammer, Singer native multiclass\cr
#'   mc-natW:    \tab Weston, Watkins native multiclass\cr
#'   one-svc:    \tab one class classification\cr
#'   eps-svr:    \tab epsilon regression\cr
#'   nu-svr:     \tab nu regression\cr
#'   eps-bsvr:   \tab bound-constraint svm regression\cr
#' }
#'
#' Pairwise multiclass can be selected for \code{C-svc} and \code{nu-svc} if
#' the label vector contains more than two classes. For
#' \code{\link[LiblineaR:LiblineaR]{LiblineaR}} the multiclass implementation
#' is always based on "one against the rest" for all SVMs except for
#' \code{mc-natC} which implements native multiclass according to Crammer and
#' Singer. The following table shows which SVM method is available in which
#' package:\cr
#'
#' \tabular{lccc}{
#'   SVM name   \tab kernlab \tab e1071 \tab LiblineaR \cr
#'   -------------------- \tab -------------- \tab -------------- \tab------
#' --------\cr
#'   C-svc:      \tab x \tab x \tab x \cr
#'   l2rl2l-svc: \tab - \tab - \tab x \cr
#'   l2rl2lp-svc:\tab - \tab - \tab x \cr
#'   l1rl2l-svc: \tab - \tab - \tab x \cr
#'   nu-svc:     \tab x \tab x \tab - \cr
#'   C-bsvc:     \tab x \tab - \tab - \cr
#'   mc-natC:    \tab x \tab - \tab x \cr
#'   mc-natW:    \tab x \tab - \tab - \cr
#'   one-svc:    \tab x \tab x \tab - \cr
#'   eps-svr:    \tab x \tab x \tab - \cr
#'   nu-svr:     \tab x \tab x \tab - \cr
#'   eps-bsvr:   \tab x \tab - \tab - \cr
#' }
#'
#'
#' SVM parameters\cr
#'
#' To avoid unnecessary changes of parameters names when switching between SVM
#' implementation in different packages unified names for identical parameters
#' are available. They are translated by KeBABS to the SVM specific name. The
#' obvious example is the cost parameter for the C-svm. It is named \code{C} in
#' \code{\link[kernlab:ksvm]{kernlab}} and \code{cost} in
#' \code{\link[e1071:svm]{e1071}} and
#' \code{\link[LiblineaR:LiblineaR]{LiblineaR}}. The unified name in KeBABS is
#' cost. If the parameter is passed to \code{kbsvm} in a package specific
#' version it is translated back to the KeBABS name internally. This applies to
#' following parameters - here shown with their unified names:\cr
#'
#' \tabular{ll}{
#'   parameter name \tab description\cr
#'   ----------------------- \tab -----------------------------------------
#' -----------\cr
#'   cost:             \tab cost parameter of C-SVM\cr
#'   nu:               \tab nu parameter of nu-SVM\cr
#'   eps:              \tab epsilon parameter of eps-SVR and nu-SVR\cr
#'   classWeights:     \tab class weights for asymmetrical class size\cr
#'   tolerance:        \tab tolerance as termination crit. for optimization\cr
#'   cross:            \tab number of folds in k-fold cross validation\cr
#' }
#'
#'
#' The following table shows the relevance of the SVM parameters cost, nu and
#' eps for the different SVMs:\cr
#'
#' \tabular{lccc}{
#'   SVM name   \tab cost \tab nu \tab eps \cr
#'   -------------------- \tab -------------- \tab -------------- \tab-----
#' ---------\cr
#'   C-svc:      \tab x \tab - \tab - \cr
#'   l1rl2l-svc: \tab x \tab - \tab - \cr
#'   l1rl2lp-svc:\tab x \tab - \tab - \cr
#'   l1rl2l-svc: \tab x \tab - \tab - \cr
#'   nu-svc:     \tab - \tab x \tab - \cr
#'   C-bsvc:     \tab x \tab - \tab - \cr
#'   mc-natC:    \tab x \tab - \tab - \cr
#'   mc-natW:    \tab x \tab - \tab - \cr
#'   one-svc:    \tab x \tab - \tab - \cr
#'   eps-svr:    \tab - \tab - \tab x \cr
#'   nu-svr:     \tab - \tab x \tab - \cr
#'   eps-bsvr:   \tab - \tab - \tab x \cr
#' }
#'
#'
#' Hint: Please be aware that identical parameter names between different SVMs
#' do not necessarily mean, that their values are also identical between
#' packages but they depend on the actual SVM formulation which could be
#' different. For example the \code{cost} parameter is identical between
#' C-SVMs in packages \code{\link[kernlab:ksvm]{kernlab}},
#' \code{\link[e1071:svm]{e1071}} and
#' \code{\link[LiblineaR:LiblineaR]{LiblineaR}} but is for example different
#' from the \code{cost} parameter in l2rl2l-svc in
#' \code{\link[LiblineaR:LiblineaR]{LiblineaR}} because the C-SVM
#' uses a linear loss but the l2rl2l-svc uses a quadratic loss.\cr\cr
#'
#' Feature weights\cr
#'
#' On user request (see parameter \code{featureWeights}) feature weights are
#' computed amd stored in the model (for a detailed description see
#' \code{\link{getFeatureWeights}}). Pruning of feature weights can be achieved
#' with the parameter \code{weightLimit} which defines the cutoff for
#' small feature weights not stored in the model.
#'
#' Hint: For training with a precomputed kernel matrix feature weights are
#' not available. For multiclass prediction is currently not performed via
#' feature weights but native in the SVM.\cr\cr
#'
#' Cross validation, grid search and model selection\cr
#'
#' Cross validation can be controlled with the parameters \code{cross} and
#' \code{noCross}. For details on cross validation see \link{crossValidation}.
#' Grid search can be performed by passing multiple SVM parameter values as
#' vector instead of a single value to \code{kbsvm}. Also multiple sequence
#' kernel objects and multiple SVMs can be used for grid search. For details
#' see \link{gridSearch}. For model selection nested cross validation is used
#' with the parameters \code{nestedCross} and \code{noNestedCross} for the
#' outer and \code{cross} and \code{noCross} for the inner cross validation.
#' For details see \link{modelSelection}.\cr\cr
#'
#' Training with feature subset\cr
#'
#' After performing feature selection repeating the learning task with a
#' feature subset can easily be achieved by specifying a feature subset
#' with the parameter \code{features} as character vector. The feature subset
#' must be a subset from the feature space of the sequence kernel passed in
#' the parameter \code{kernel}. Grid search and model selection with a feature
#' subset can only be used for a single sequence kernel object in the parameter
#' \code{kernel}.\cr\cr
#' Hint: For normalized kernels all features of the feature space are used for
#' normalization not just the feature subset. For a normalized motif kernel
#' (see \code{\link{motifKernel}}) only the features listed in the motif list
#' are part of the feature space. Therefore the motif kernel defined with the
#' same feature subset leads to a different result in the normalized case.\cr\cr
#'
#' Probability model\cr
#'
#' SVMs from the packages \code{\link[kernlab:ksvm]{kernlab}} and
#' \code{\link[e1071:svm]{e1071}} support the generation of a probability
#' model using Platt scaling (for details see
#' \code{\link[kernlab:ksvm]{kernlab}},
#' \code{\link[kernlab:predict.ksvm]{predict.ksvm}},
#' \code{\link[e1071:svm]{svm}} and
#' \code{\link[e1071:predict.svm]{predict.svm}})
#' allowing the computation of class probabilities during prediction. The
#' parameter \code{probabilityModel} controls the generation of a probability
#' model during training (see also parameter \code{predictionType} in
#' \code{\link{predict}}).
#'
#' @return
#' kbsvm: upon successful completion, the function returns a model of class
#' \code{\linkS4class{KBModel}}. Results for cross validation can be retrieved
#' from this model with the accessor \code{\link{cvResult}}, results for grid
#' search or model selection with \code{\link{modelSelResult}}. In case of 
#' model selection the results of the outer cross validation loop can be 
#' retrieved with with the accessor \code{\link{cvResult}}.
#'
#' @seealso \code{\link{predict}}, \code{\link{getKernelMatrix}},
#' \code{\link{getExRep}}, \code{\link{kernelParameters-method}},
#' \code{\link{spectrumKernel}}, \code{\link{mismatchKernel}},
#' \code{\link{gappyPairKernel}}, \code{\link{motifKernel}},
#' \code{\link{getFeatureWeights}}
#'
#' @examples
#'
#' ## load transcription factor binding site data
#' data(TFBS)
#' enhancerFB
#' ## we use 70 of the samples for training and the rest for test
#' train <- sample(1:length(enhancerFB), length(enhancerFB) * 0.7)
#' test <- c(1:length(enhancerFB))[-train]
#' ## create the kernel object for dimers without normalization
#' specK2 <- spectrumKernel(k=2)
#' ## show details of kernel object
#' specK2
#'
#' ## run training with kernel matrix
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=specK2,
#'                pkg="kernlab", svm="C-svc", C=10, explicit="no")
#'
#' ## show KeBABS model
#' model
#' ## show class of KeBABS model
#' class(model)
#' ## show native SVM model contained in KeBABS model
#' svmModel(model)
#' ## show class of native SVM model
#' class(svmModel(model))
#'
#' \dontrun{
#' ## examples for package and SVM selection
#' ## now run the same samples with the same kernel on e1071 which
#' ## currently only supports an explicit representation
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=specK2,
#'                pkg="e1071", svm="C-svc", C=10, explicit="yes")
#'
#' ## show KeBABS model
#' model
#' ## show native SVM model contained in KeBABS model
#' svmModel(model)
#' ## show class of native SVM model
#' class(svmModel(model))
#'
#' ## run the same samples with the same kernel on e1071 with nu-SVM
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=specK2,
#'                pkg="e1071", svm="nu-svc",nu=0.7, explicit="yes")
#'
#' ## show KeBABS model
#' model
#'
#' ## training with feature weights
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=specK2,
#'                pkg="e1071", svm="C-svc", C=10, explicit="yes",
#'                featureWeights="yes")
#'
#' ## show feature weights
#' dim(featureWeights(model))
#' featureWeights(model)[,1:5]
#'
#' ## training without feature weights
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=specK2,
#'                pkg="e1071", svm="C-svc", C=10, explicit="yes",
#'                featureWeights="no")
#'
#' ## show feature weights
#' featureWeights(model)
#'
#' ## pruning of feature weights
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=specK2,
#'                pkg="e1071", svm="C-svc", C=10, explicit="yes",
#'                featureWeights="yes", weightLimit=0.5)
#'
#' dim(featureWeights(model))
#'
#' ## training with precomputed kernel matrix
#' km <- getKernelMatrix(specK2, x=enhancerFB, selx=train)
#' model <- kbsvm(x=km, y=yFB[train], kernel=specK2,
#'                pkg="kernlab", svm="C-svc", C=10, explicit="no")
#'
#' ## training with precomputed explicit representation
#' exrep <- getExRep(enhancerFB, sel=train, kernel=specK2)
#' model <- kbsvm(x=exrep, y=yFB[train], kernel=specK2,
#'                pkg="e1071", svm="C-svc", C=10, explicit="yes")
#'
#' ## computing of probability model via Platt scaling during training
#' ## in prediction class membership probabilities can be computed
#' ## from this probability model
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=specK2,
#'                pkg="e1071", svm="C-svc", C=10, explicit="yes",
#'                probModel=TRUE)
#'
#' ## show parameters of the fitted probability model which are the parameters
#' ## probA and probB for the fitted sigmoid function in case of classification
#' ## and the value sigma of the fitted Laplacian in case of a regression
#' probabilityModel(model)
#'
#' ## cross validation, grid search and model selection are also performed
#' ## via the kbsvm method. Examples can be found on the respective help pages
#' ## (see Details section)
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}
#' @keywords kbsvm
#' @keywords training
#' @keywords cross validation
#' @keywords grid search
#' @keywords model selection
#' @keywords feature weights
#' @keywords methods

#' @rdname kbsvm-methods
#' @aliases
#' kbsvm
#' kbsvm,BioVector-method
#' @export
#'

## no method with XString - at least 2 seqs needed for training
setMethod("kbsvm",signature(x="BioVector"), kbsvm.seqs)

#' @rdname kbsvm-methods
#' @aliases
#' kbsvm,XStringSet-method
#' @export

setMethod("kbsvm",signature(x="XStringSet"), kbsvm.seqs)

kbsvm.ExplicitRep <- function(x, y, kernel=NULL, pkg="auto", svm="C-svc",
            explicit="auto",  explicitType="auto", featureType="linear",
            featureWeights="auto", weightLimit=.Machine$double.eps,
            classWeights = numeric(0), cross=0, noCross=1, groupBy=NULL,
            nestedCross=0, noNestedCross=1, perfParameters=character(0),
            perfObjective="ACC", probModel=FALSE, sel=integer(0),
            showProgress=FALSE, showCVTimes=FALSE, runtimeWarning=TRUE,
            verbose = getOption("verbose"), ...)
{
    addArgs <- list(...)

    if (length(perfParameters) >= 1 && perfParameters[1] == "ALL")
        perfParameters <- c("ACC", "BACC", "MCC")

    ## for explicit representation subsetting on user level only
    ## allows more efficient handling on user level
    if (length(sel) > 0)
        stop("please subset the explicit representation before the training\n")

    if (x@quadratic)
        stop("please pass linear explicit representation to training\n")

    ## check parameters and allocate model
    model <- checkKBSVMParams(x=x, y=y, kernel=kernel, svm=svm, pkg=pkg,
                    explicit=explicit, explicitType=explicitType,
                    featureType=featureType, featureWeights=featureWeights,
                    weightLimit=weightLimit, classWeights=classWeights,
                    probModel=probModel, sel=sel, cross=cross, noCross=noCross,
                    groupBy=groupBy, nestedCross=nestedCross,
                    noNestedCross=noNestedCross, perfParameters=perfParameters,
                    perfObjective=perfObjective, runtimeWarning=runtimeWarning,
                    addArgs=addArgs)

    model@call <- deparse(sys.call(-1))

    if (!is.null(model@modelSelResult))
    {
        return(performModelSelection(object=x, model=model, y=y, sel=sel,
                    explicit=explicit, featureWeights=featureWeights,
                    weightLimit=weightLimit, groupBy=groupBy,
                    perfParameters=perfParameters, perfObjective=perfObjective,
                    showProgress=showProgress, showCVTimes=showCVTimes,
                    addArgs=addArgs, verbose=verbose))
    }

    ## select SVM and method
    result <- selectSVMMethod(model, kebabsInfo@supportedPkgs, x)

    if (length(result) == 2)
    {
        model@svmInfo <- result[[1]]
        model@ctlInfo <- result[[2]]
    }
    else
        stop("wrong result from selectSVMMethod\n")

    ## generate runtime warnings
    model@ctlInfo <- generateRuntimeMessage(model, x)

    ## convert parameters to target SVM
    ## and remove unneeded parameters
    model@svmInfo <- convertSVMParameters(model)

    if (model@ctlInfo@selMethod == "KernelMatrix")
    {
        ## LiblineaR does not support kernel matrix => quadratic
        ## ER not needed in this case

        km <- linearKernel(x=x, selx=sel)
        sel <- integer(0)

        model <- kbsvm(x=km, y=y, kernel=kernel, svm=svm, pkg=pkg,
                       explicit=explicit, explicitType=explicitType,
                       featureWeights=featureWeights, cross=cross,
                       noCross=noCross, groupBy=groupBy, verbose=verbose,
                       weightLimit=weightLimit, sel=sel, ...)

        return(model)
    }
    else
    {
        if (model@svmInfo@reqFeatureType == "quadratic" &&
            model@svmInfo@selPackage == "LiblineaR" && cross == 0)
        {
            ## convert to explicit representation for quadratic kernel
            ## as quadratic kernel is not available in LiblineaR
            ## performCrossValidation is called with linear exrep
            exRepLin <- x
            x <- getExRepQuadratic(x)
        }

        ## store feature names in model, needed for direct prediction
        ## via ex rep in svm without feature weights
        model@trainingFeatures <- colnames(x)

        ## convert parameters to target SVM

        ## remove orig kernel and set up ER kernel


        if (cross != 0)
        {
            if (showCVTimes)
            {
                timeCV <- system.time(
                    cvModel <- performCrossValidation(object=x, x=NULL, y=y,
                                    sel=sel, model=model, cross=cross,
                                    noCross=noCross, groupBy=groupBy,
                                    perfParameters=perfParameters,
                                    verbose=verbose))

                cat("\nCross Validation Time:\n")
                print(timeCV)
                cat("\n")
                return(cvModel)
            }
            else
            {
                return(performCrossValidation(object=x, x=NULL, y=y,
                                    sel=sel, model=model, cross=cross,
                                    noCross=noCross, groupBy=groupBy,
                                    perfParameters=perfParameters,
                                    verbose=verbose))
            }
        }
        else
        {
            ## invoke selected svm
            allArgs <- model@svmInfo@selSVMPar
            allArgs[["x"]] <- x
            allArgs[["y"]] <- y
            allArgs[["svmInfo"]] <- model@svmInfo
            allArgs[["verbose"]] <- verbose
            model@svmModel <- do.call(trainSVM, allArgs)

            if (is.null(model@svmModel))
            {
                stop(paste(model@svmInfo@selPackage, "-",
                           model@svmInfo@selSVM, "did not return a model\n"))
            }

            if (model@ctlInfo@classification == TRUE)
            {
                model@numClasses <- getSVMSlotValue("numClasses", model)
                model@classNames <- getSVMSlotValue("classNames", model)

                if (model@numClasses > 2)
                    model@ctlInfo@multiclassType <- getMulticlassType(model)

                if (model@svmInfo@probModel == TRUE)
                {
                    model@probA <- getSVMSlotValue("probA", model)
                    model@probB <- getSVMSlotValue("probB", model)
                }
            }
            else
            {
                ## copy param for regression with probability model
                if (model@svmInfo@probModel == TRUE)
                    model@sigma <- getSVMSlotValue("sigma", model)
            }

            ind <- getSVMSlotValue("svIndex", model)
            exRepSV <- NULL

            if (length(ind) > 0)
            {
                model@SV <- x[ind]
                model@svIndex <- ind
                model@alphaIndex <- getSVMSlotValue("alphaIndex", model)

                if (model@svmInfo@reqFeatureType == "quadratic" &&
                    model@svmInfo@selPackage == "LiblineaR")
                    exRepSV <- exRepLin[ind,]
                else
                    exRepSV <- x[ind,]
            }

            if (featureWeights != "no")
            {
                model@featureWeights <-
                    getFeatureWeights(model=model,
                                      exrep=exRepSV,
                                      weightLimit=model@svmInfo@weightLimit)

                model@b <- getSVMSlotValue("b", model)

                if (length(model@featureWeights) < 1)
                    stop("no feature weights available\n")
            }

            return(model)
        }
    }
}

#' @rdname kbsvm-methods
#' @aliases
#' kbsvm,ExplicitRepresentation-method
#' @export
#'

setMethod("kbsvm",signature(x="ExplicitRepresentation"), kbsvm.ExplicitRep )

kbsvm.KernelMatrix <- function(x, y, kernel=NULL, pkg="auto", svm="C-svc",
            explicit="no",  explicitType="auto", featureType="linear",
            featureWeights="no", classWeights = numeric(0), cross=0, noCross=1,
            groupBy=NULL, nestedCross=0, noNestedCross=1,
            perfParameters=character(0), perfObjective="ACC", probModel=FALSE,
            sel=integer(0), showProgress=FALSE, showCVTimes=FALSE,
            runtimeWarning=TRUE, verbose = getOption("verbose"), ...)
{
    addArgs <- list(...)

    if (length(perfParameters) >= 1 && perfParameters[1] == "ALL")
        perfParameters <- c("ACC", "BACC", "MCC")

    ## for kernel matrix subsetting on user level only
    ## allows more efficient handling on user level
    if (length(sel) > 0)
        stop("please subset the kernel matrix before the training\n")

    ## parameters explicit and featureWeights included to avoid
    ## getting them into the transparent parameters

    if (!(explicit %in% c("no", "auto")))
    {
        stop("processing via explicit representation not possible\n",
             "        for call with kernel matrix\n")
    }

    if (featureWeights != "no")
    {
        stop("prediction via feature weights is not possible\n",
             "        for call with kernel matrix\n")

    }

    weightLimit <- .Machine$double.eps

    ## check parameters and allocate model
    model <- checkKBSVMParams(x=x, y=y, kernel=kernel, svm=svm, pkg=pkg,
                    explicit=explicit, explicitType=explicitType,
                    featureType=featureType, featureWeights=featureWeights,
                    weightLimit=weightLimit, classWeights=classWeights,
                    probModel=probModel, sel=sel, cross=cross, noCross=noCross,
                    groupBy=groupBy, nestedCross=nestedCross,
                    perfParameters=perfParameters, perfObjective=perfObjective,
                    noNestedCross=noNestedCross, runtimeWarning=TRUE,
                    addArgs=addArgs)

    model@call <- deparse(sys.call(-1))

    if (!is.null(model@modelSelResult))
    {
        return(performModelSelection(object=x, model=model, y=y, sel=sel,
                    explicit=explicit, featureWeights=featureWeights,
                    weightLimit=weightLimit, groupBy=groupBy,
                    perfParameters=perfParameters, perfObjective=perfObjective,
                    showProgress=showProgress, showCVTimes=showCVTimes,
                    addArgs=addArgs, verbose=verbose))
    }

    ## select SVM and method
    result <- selectSVMMethod(model, kebabsInfo@supportedPkgs, x)

    if (length(result) == 2)
    {
        model@svmInfo <- result[[1]]
        model@ctlInfo <- result[[2]]
    }
    else
        stop("wrong result from selectSVMMethod\n")

    ## generate runtime warnings
    model@ctlInfo <- generateRuntimeMessage(model, x)

    ## convert parameters to target SVM
    model@svmInfo <- convertSVMParameters(model)

    # change parameters if necessary and remove unneeded parameters
    ## $$$ TODO modify dotArgs
    ## finally do.call(<function name>, named pars, dotArgs)

    if (model@ctlInfo@selMethod == "KernelMatrix")
    {
        if (cross != 0)
        {
            if (showCVTimes)
            {
                timeCV <- system.time(
                    cvModel <- performCrossValidation(object=x, x=NULL, y=y,
                                    sel=sel, model=model, cross=cross,
                                    noCross=noCross, groupBy=groupBy,
                                    perfParameters=perfParameters,
                                    verbose=verbose))

                cat("\nCross Validation Time:\n")
                print(timeCV)
                cat("\n")
                return(cvModel)
            }
            else
            {
                return(performCrossValidation(object=x, x=NULL, y=y,
                                    sel=sel, model=model, cross=cross,
                                    noCross=noCross, groupBy=groupBy,
                                    perfParameters=perfParameters,
                                    verbose=verbose))
        }
        }
        else
        {
            ## invoke selected svm
            allArgs <- model@svmInfo@selSVMPar
            allArgs[["x"]] <- x
            allArgs[["y"]] <- y
            allArgs[["svmInfo"]] <- model@svmInfo
            allArgs[["verbose"]] <- verbose
            model@svmModel <- do.call(trainSVM, allArgs)

            if (is.null(model@svmModel))
            {
                stop(paste(model@svmInfo@selPackage, "-",
                           model@svmInfo@selSVM, "did not return a model\n"))
            }

            if (model@ctlInfo@classification == TRUE)
            {
                model@numClasses <- getSVMSlotValue("numClasses", model)
                model@classNames <- getSVMSlotValue("classNames", model)

                if (model@numClasses > 2)
                    model@ctlInfo@multiclassType <- getMulticlassType(model)

                if (model@svmInfo@probModel == TRUE)
                {
                    model@probA <- getSVMSlotValue("probA", model)
                    model@probB <- getSVMSlotValue("probB", model)
                }
            }
            else
            {
                ## copy param for regression with probability model
                if (model@svmInfo@probModel == TRUE)
                    model@sigma <- getSVMSlotValue("sigma", model)
            }

            ind <- getSVMSlotValue("svIndex", model)

            ## just store the sv indices, because the data vectors are missing
            if (length(ind) > 0)
            {
                model@svIndex <- ind
                model@alphaIndex <- getSVMSlotValue("alphaIndex", model)
            }
            return(model)
        }
    }
    else
        stop("explict representation not available for kernel matrix\n")
}

#' @rdname kbsvm-methods
#' @aliases
#' kbsvm,KernelMatrix-method
#' @export
#'

setMethod("kbsvm",signature(x="KernelMatrix"), kbsvm.KernelMatrix)

kbsvm.kernelMatrix <- function(x, y, kernel=NULL, svm="C-svc", pkg="auto",
            explicit="no", explicitType="auto", featureType="linear",
            featureWeights="no", classWeights=numeric(0), cross=0,
            noCross=1, groupBy=NULL, nestedCross=0, noNestedCross=1,
            perfParameters="ACC", perfObjective="ACC", probModel=FALSE,
            weightLimit=.Machine$double.eps, sel=integer(0),
            showProgress=FALSE, showCVTimes=FALSE,
            verbose = getOption("verbose"), ...)
{
    kbsvm(x=as.KernelMatrix(x), y=y, kernel=kernel, svm=svm, pkg=pkg,
          explicit=explicit, explicitType=explicitType,
          featureType=featureType, featureWeights=featureWeights,
          classWeights=classWeights, cross=cross, noCross=noCross,
          groupBy=groupBy, nestedCross=nestedCross, noNestedCross=noNestedCross,
          perfParameters=perfParameters, perfObjective=perfObjective,
          probModel=probModel, weightLimit=weightLimit, sel=sel,
          showProgress=showProgress, showCVTimes=showCVTimes,
          verbose=verbose, ...)
}

#setMethod("kbsvm",signature(x="kernelMatrix"), kbsvm.kernelMatrix)
