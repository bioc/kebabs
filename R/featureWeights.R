##345678901234567890123456789012345678901234567890123456789012345678901234567890

convertSymmetricMatrixToVector <- function(x)
{
    if (!is(x, "matrix") || !isSymmetric(x))
    stop("'x' must be a symmetric matrix\n")

    x <- as(x, "matrix")
    numRows <- nrow(x)
    result <- rep(0, (numRows * (numRows + 1) / 2))
    sqrt2 <- sqrt(2)

    j <- 1
    for (i in 1:numRows)
    {
        result[j] <- x[i,i]

        if (i < numRows)
            result[(j+1):(j + numRows - i)] <- sqrt2 * x[i,((i+1):numRows)]

        names(result)[j:(j + numRows - i)] <-
        paste(colnames(x)[i], colnames(x)[i:numRows], sep="_")
        j <- j + numRows - i + 1
    }

    return(result)
}

getFeatureWeightsPosIndep <- function(model, exrep=NULL, svmIndex=1,
               features=NULL, weightLimit=.Machine$double.eps)
{
    ## $$$ TODO
    ## adapt feature weights for PSVM

    ## $$$ TODO adapt for sparse exrep / sparse weights / sparse pruning


    if (is(model@svmModel, "LiblineaR"))
    {
        bias <- getSVMSlotValue("bias", model)
        weights <- getSVMSlotValue("weights", model)
        noWeights <- ncol(weights)

        if (bias)
        {
            noWeights <- noWeights - 1
            weights <- weights[svmIndex, 1:noWeights, drop=FALSE]
        }

        classNames <- getSVMSlotValue("classNames", model)

        ## assign class names for multiclass with one against the rest
        if (length(classNames) > 2)
            rownames(weights) <- classNames[svmIndex]

        if (colnames(model@svmModel$W)[1] != "W1")
        {
            colnames(weights) <-
                colnames(model@svmModel$W)[1:noWeights]
        }
        else
        {
            if (length(model@trainingFeatures) > 0)
                colnames(weights) <- model@trainingFeatures
        }

        ## prune weights
        if (weightLimit > 0)
        {
            weights[which(abs(weights) < weightLimit)] <- 0
            weights <- weights[,which(weights != 0), drop=FALSE]
        }

        return(weights)
    }

    ## for dense ER the function must be called with the compact form
    ## no check to avoid unnecessary performance loss

    ## adaptation of coef and b to unique format across all svms is done
    ## in access routine getSVMSlotValue
    coef <- getSVMSlotValue("coef", model)

    ## subset to requested model
    coef <- coef[, svmIndex]

    if (!is.numeric(coef) || length(coef) != nrow(exrep))
        stop("invalid coefficients returned from model\n")

    weights <- NULL

    if (model@svmInfo@reqFeatureType == "quadratic")
    {
        if (is(exrep, "ExplicitRepresentationSparse"))
            weights <- ersTransposedAsdgCMatrix(exrep) %*% (coef * exrep)
        else
            weights <- t(exrep) %*% (coef * exrep)

        rownames(weights) <- colnames(exrep)
        colnames(weights) <- colnames(exrep)

        tempWeights <- convertSymmetricMatrixToVector(as(weights, "matrix"))
        weights <- matrix(tempWeights, nrow=1)
        colnames(weights) <- names(tempWeights)

        ## prune weights
        if (weightLimit > 0)
            weights <- weights[,which(abs(weights) > weightLimit), drop=FALSE]

        return(weights)
    }

    ## $$$ TODO feature weights for Weston/Watkins
    tempWeights <- t(coef) %*% exrep

    ## prune weights
    if (weightLimit > 0)
    {
        tempWeights[which(abs(tempWeights) < weightLimit)] <- 0
        nonZero <- which(colSums(tempWeights) != 0)
        weights <- tempWeights[, nonZero, drop=FALSE]
    }
    else
        weights <- tempWeights

    return(weights)
}

getFeatureWeightsPosDep <- function(model, svmIndex=1, features=NULL,
                                    weightLimit=.Machine$double.eps)
{
    if (missing(model) || !is(model, "KBModel"))
        ("model' must be a model object of class \"KBModel\"\n")

    if (is.null(model@SV))
    {
        stop("support vectors are missing in model for computation of\n",
             "        position specific feature weights\n")
    }

    if (inherits(model@SV, "ExplicitRepresentationresentation"))
    {
        stop("position specific feature weights cannot be generated from the\n",
             "         explicit representation\n")
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
        k <- 0
        motifs <- kernelParameters(model@svmInfo@selKernel)$motifs
        motifLengths <- kernelParameters(model@svmInfo@selKernel)$motifLengths
        maxMotifLength <- max(motifLengths)
        maxPatternLength <- max(nchar(motifs))

        ## rough limit for no of nodes in motif tree from no of
        ## chars and no of substitution groups, add one for root
        nodeLimit <- sum(motifLengths) + 1 +
        sum(sapply(gregexpr("[", motifs, fixed=TRUE),
                   function(x) length(unlist(x))))
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

    if (!class(model@svmInfo@selKernel) %in% names(kernelType))
        stop("wrong kernel class\n")

    distWeight <- kernelParameters(model@svmInfo@selKernel)$distWeight

    if (isTRUE(all.equal(distWeight, c(1, rep(0, length(distWeight)-1)))))
        posSpec <- TRUE
    else
        posSpec <- FALSE

    normalized <- kernelParameters(model@svmInfo@selKernel)$normalized
    ignoreLower <- kernelParameters(model@svmInfo@selKernel)$ignoreLower
    presence <- kernelParameters(model@svmInfo@selKernel)$presence
    reverseComplement <- kernelParameters(model@svmInfo@selKernel)$revComplement
    bioCharset <- getBioCharset(model@SV, TRUE)
    unmapped <- is(model@SV, "DNAStringSet") || is(model@SV, "RNAStringSet")
    isXStringSet <- is(model@SV, "XStringSet")
    maxSeqLength <- max(width(model@SV))
    offsetSV <- mcols(model@SV)[["offset"]]

    if (is.null(offsetSV))
    {
        minPos <- 1
        maxPos <- max(width(model@SV))
        offsetSV <- integer(0)
    }
    else
    {
        minPos <- min(offsetSV)
        maxPos <- max(offsetSV + width(model@SV) - 1)
        offsetSV <- offsetSV[model@svIndex]
    }

    coefs <- getSVMSlotValue("coef", model)

    if (!is.list(model@alphaIndex))
        svIndices <- which(model@svIndex %in% model@alphaIndex)
    else
        svIndices <- which(model@svIndex %in% model@alphaIndex[[svmIndex]])

    on.exit(.Call("freeHeapFeatureWeightsC"))

    ## distanceWeight is not really relevant here because actual distance
    ## weighting is done at prediction / generation of prediction profile

    ## get position specific feature weights
    posDepFeatureWeights <- .Call("getFeatureWeightsPosDepC",
            model@SV, svIndices - 1, offsetSV,
            as.logical(isXStringSet), as.integer(maxSeqLength),
            as.integer(svmIndex), as.double(weightLimit),
            as.integer(kernelType[class(model@svmInfo@selKernel)]),
            as.integer(k), as.integer(m), as.integer(bioCharset[[2]]), motifs,
            motifLengths, as.integer(maxMotifLength),
            as.integer(maxPatternLength), as.integer(nodeLimit), coefs,
            as.logical(reverseComplement), as.logical(posSpec),
            as.integer(minPos), as.integer(maxPos), as.logical(normalized),
            as.logical(!ignoreLower), as.logical(unmapped))

    return(posDepFeatureWeights)
}

#' @rdname featureWeights
#' @title Feature Weights
#'
#' @description Compute Feature Weights for KeBABS Model
#'
#' @param model model object of class \code{\linkS4class{KBModel}} created
#' by \code{\link{kbsvm}}.
#'
#' @param exrep optional explicit representation of the support vectors from
#' which the feature weights should be computed. If no explicit representation
#' is passed to the function the explicit representation is generated internally
#' from the support vectors stored in the model. default=\code{NULL}
#'
#' @param features feature subset of the specified kernel in the form of a
#' character vector. When a feature subset is passed to the function all other
#' features in the feature space are not considered for the explicit
#' representation. (see below) default=\code{NULL}
#'
#' @param weightLimit the feature weight limit is a single numeric value and
#' allows pruning of feature weights. All feature weights with an absolute
#' value below this limit are set to 0 and are not considered in the feature
#' weights. Default=.Machine$double.eps
#'
#' @details
#' Overview\cr\cr
#' Feature weights represent the contribution to the decision value for a
#' single occurance of the feature in the sequence. In this way they give a
#' hint concerning the importance of the individual features for a given
#' classification or regression task. Please consider that for a pattern length
#' larger than 1 patterns at neighboring sequence positions overlap and are no
#' longer independent from each other. Apart from the obvious overlapping
#' possibility of patterns for e.g. gappy pair kernel, motif kernel or mixture
#' kernels multiple patterns can be relevant for a single position. Therefore
#' feature weights do not describe the relevance for individual features
#' exactly.\cr\cr
#' Computation of feature weights\cr\cr
#' Feature weights can be computed automatically as part of the training (see
#' parameter \code{featureWeights} in method \code{\link{kbsvm}}. In this case
#' the function getFeatureWeights is called during training automatically.
#' When this parameter is not set during training computation of feature weights
#' after training is possible with the function getFeatureWeights. The function
#' also supports pruning of feature weights (see parameter \code{weightLimit}
#' allowing to test different prunings without retraining.\cr\cr
#' Usage of feature weights\cr\cr
#' Feature weights are used during prediction to speed up the prediction
#' process. Prediction via feature weights is performed in KeBABS when feature
#' weights are available in the model (see \code{\link{featureWeights}}).
#' When feature weights are not available or for multiclass prediction KeBABS
#' defaults to the native prediction in the SVM used during training.\cr\cr
#' Feature weights are also used during generation of prediction profiles
#' (see \code{\link{getPredictionProfile}}). In the feature weights the general
#' relevance of features is reflected. When generating prediction profiles
#' for a given set of sequences from the feature weights the relevance of
#' single sequence positions is shown for the individual sequences according
#' to the given learning task.\cr\cr
#' Feature weights for position dependent kernels\cr\cr
#' For position dependent kernels the generation of feature weights is not
#' possible during training. In this case the featureWeights slot in the model
#' contains a data representation that allows simple computation of feature
#' weights during prediction or during generation of prediction profiles.
#' @return
#' Upon successful completion, the function returns the feature weights as
#' numeric vector. For quadratic kernels a matrix of feature weights is returned
#' giving the feature weights for pairs of features. In case of multiclass
#' the function returns the feature weights for the pairwise SVMs as list of
#' numeric vectors (or matrices for quadratic kernels).
#'
#' @seealso \code{\link{kbsvm}}, \code{\link{predict}},
#' \code{\link{getPredictionProfile}} \code{\link{featureWeights}},
#' \code{\linkS4class{KBModel}}
#'
#' @examples
#'
#' ## standard method to create feature weights automatically during training
#' ## model <- kbsvm( .... , featureWeights="yes", .....)
#' ## this example describes the case where feature weights were not created
#' ## during training but should be added later to the model
#'
#' ## load example sequences and select a small set of sequences
#' ## to speed up training for demonstration purpose
#' data(TFBS)
#' ## create sample indices of training and test subset
#' train <- sample(1:length(yFB), 200)
#' test <- c(1:length(yFB))[-train]
#' ## determin all labels
#' allLables <- unique(yFB)
#'
#' ## create a kernel object
#' gappyK1M4 <- gappyPairKernel(k=1, m=4)
#'
#' ## model is trainded with creation of feature weights
#' model <- kbsvm(enhancerFB[train], yFB[train], gappyK1M4,
#'                pkg="LiblineaR", svm="C-svc", cost=20)
#'
#' ## feature weights included in model
#' featureWeights(model)
#'
#' \dontrun{
#' ## model is originally trainded without creation of feature weights
#' model <- kbsvm(enhancerFB[train], yFB[train], gappyK1M4,
#'                pkg="LiblineaR", svm="C-svc", cost=20, featureWeights="no")
#'
#' ## no feature weights included in model
#' featureWeights(model)
#'
#' ## later after training add feature weights and model offset of model to
#' ## KeBABS model
#' featureWeights(model) <- getFeatureWeights(model)
#' modelOffset(model) <- getSVMSlotValue("b", model)
#'
#' ## show a part of the feature weights and the model offset
#' featureWeights(model)[1:7]
#' modelOffset(model)
#'
#' ## another scenario for getFeatureWeights is to test the performance
#' ## behavior of different prunings of the feature weights
#'
#' ## show histogram of full feature weights
#' hist(featureWeights(model), breaks=30)
#'
#' ## show number of features
#' length(featureWeights(model))
#'
#' ## first predict with full feature weights to see how performance
#' ## when feature weights are included in the model prediction is always
#' ## performed with the feature weights
#' ## changes through pruning
#' pred <- predict(model, enhancerFB[test])
#' evaluatePrediction(pred, yFB[test], allLabels=allLables)
#'
#' ## add feature weights with pruning to absolute values larger than 0.6
#' ## model offset was assigned above and is not impacted by pruning
#' featureWeights(model) <- getFeatureWeights(model, weightLimit=0.6)
#'
#' ## show histogram of full feature weights
#' hist(featureWeights(model), breaks=30)
#'
#' ## show reduced number of features
#' length(featureWeights(model))
#'
#' ## now predict with pruned feature weights
#' pred <- predict(model, enhancerFB, sel=test)
#' evaluatePrediction(pred, yFB[test], allLabels=allLables)
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}
#' @keywords feature weights
#' @keywords methods
#' @export

getFeatureWeights <- function(model, exrep=NULL, features=NULL,
                              weightLimit=.Machine$double.eps)
{
    if (missing(model) || !is(model, "KBModel"))
        stop("model' must be a model object of class \"KBModel\"\n")

    if (!isSingleNumber(weightLimit))
        stop("'weightLimit' must be a number\n")

    if (!is.null(exrep) && !inherits(exrep, "ExplicitRepresentation"))
        stop("'exrep' must be an explicit representation\n")

    if (!is.null(exrep) && (exrep@quadratic != FALSE))
        stop("'exrep' must be a linear explicit representation\n")

    if (length(kernelParameters(model@svmInfo@selKernel)$distWeight) > 0)
    {
        if (length(model@ctlInfo@multiclassType) > 0 &&
            model@ctlInfo@multiclassType %in% c("pairwise", "oneAgainstRest",
                                                "CrammerSinger"))
        {
            featureWeights <- list()

            for (i in 1:choose(model@numClasses, 2))
            {
                if (model@svmInfo@selPackage != "LiblineaR")
                    svmNames <- colnames(getSVMSlotValue("coef", model))
                else
                {
                    svmNames <- getSVMSlotValue("classNames", model, raw=TRUE)

                    if (model@ctlInfo@multiclassType == "oneAgainstRest")
                        svmNames <- paste(svmNames, "/R", sep="")
                }

                featureWeights[[svmNames[i]]] <-
                    getFeatureWeightsPosDep(model=model, svmIndex=i,
                                weightLimit=weightLimit, features=features)

                ## assign name of svm
                rownames(featureWeights[[svmNames[i]]]) <- svmNames[i]
            }

            return(featureWeights)
        }
        else
        {
            return(getFeatureWeightsPosDep(model=model, weightLimit=weightLimit,
                                           features=features))
        }
    }

    if (is.null(exrep) && !is(model@svmModel, "LiblineaR"))
    {
        if (length(model@SV) < 1)
            stop("missing support vectors for feature weight computation\n")

        if (inherits(model@SV, "ExplicitRepresentation"))
            exrep <- model@SV
        else if (is(model@SV, "BioVector") || is(model@SV, "XStringSet"))
        {
            exrep <- getExRep(x=model@SV, kernel=model@svmInfo@selKernel,
                              sparse=model@ctlInfo@sparse,
                              features=features)
        }
        else
            stop("Feature weights cannot be computed from kernel matrix\n")
    }

    if (length(model@ctlInfo@multiclassType) > 0 &&
        model@ctlInfo@multiclassType %in% c("pairwise", "oneAgainstRest",
                                            "CrammerSinger"))
    {
        featureWeights <- list()

        for (i in 1:choose(model@numClasses, 2))
        {
            if (model@svmInfo@selPackage != "LiblineaR")
                svmNames <- colnames(getSVMSlotValue("coef", model))
            else
            {
                svmNames <- getSVMSlotValue("classNames", model, raw=TRUE)

                if (model@ctlInfo@multiclassType == "oneAgainstRest")
                   svmNames <- paste(svmNames, "/R", sep="")
            }

            featureWeights[[svmNames[i]]] <-
                getFeatureWeightsPosIndep(model=model, exrep=exrep,
                                          svmIndex=i, features=features,
                                          weightLimit=weightLimit)
            ## assign name of svm
            rownames(featureWeights[[svmNames[i]]]) <- svmNames[i]
        }

        return(featureWeights)
    }
    else
    {
        return(getFeatureWeightsPosIndep(model=model, exrep=exrep,
                    features=features, weightLimit=weightLimit))
    }
}
