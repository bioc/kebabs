##345678901234567890123456789012345678901234567890123456789012345678901234567890
getPredictionProfile.BioVectorOrXSS <- function(object, kernel, featureWeights,
                      b, svmIndex=1, sel=NULL, weightLimit=.Machine$double.eps)
{
    ## $$$ TODO clarify if prediction profiles are relevant for quadratic kernel
    if (!(inherits(object, "XStringSet") || inherits(object, "BioVector")))
        stop("Prediction profiles can only be computed for sequence data\n")

    if (length(object) < 1)
        stop("missing sequence data\n")

    if (!is(kernel, "SequenceKernel"))
        stop("'kernel' must be a sequence kernel\n")

    if (isUserDefined(kernel))
        stop("prediction profiles not supported for user-defined kernel\n")
    
    if (is.null(featureWeights))
        stop("'featureWeights' are missing\n")

    if (is.vector(featureWeights, mode="numeric"))
    {
        if (length(names(featureWeights)) == length(featureWeights))
        {
            featureWeights <- matrix(featureWeights, nrow=1,
                                     dimnames=list(NULL, names(featureWeights)))
        }
        else
            stop("missing feature names in 'featureWeights'\n")
    }

    if (!(is.matrix(featureWeights) || is(featureWeights, "dgCMatrix") ||
          is.list(featureWeights)))
        stop("'featureWeights' must be a matrix or a list of matrices\n")

    if (!is.numeric(b))
        stop("'b' must be a single numeric\n")

    if (!is.null(sel))
    {
        if (!is.numeric(sel) || length(sel) > length(object) ||
            any(sel < 0) || any (sel > length(object)))
            stop("'sel' must be an integer vector with indices into 'object'\n")
        sel <- as.integer(sel)
    }
    else
        sel <- 1:length(object)

    seqLength <- width(object[sel])[1]
    maxSeqLength <- max(width(object[sel]))

    predProf <- new("PredictionProfile")
    predProf@sequences <- object[sel]
    predProf@kernel <- kernel
    predProf@baselines <- - b[svmIndex] / width(object[sel])

    kernelType <- 1:6
    names(kernelType) <- c("SpectrumKernel",
                           "MixedSpectrumKernel",
                           "MismatchKernel",
                           "MotifKernel",
                           "WeightedDegreeKernel",
                           "GappyPairKernel")

    if (!class(kernel) %in% names(kernelType))
        stop("wrong kernel class\n")

    if (!is(kernel, "MismatchKernel"))
        annSpec <- kernelParameters(kernel)$annSpec
    else
        annSpec <- FALSE

    offsetX <- mcols(object)[["offset"]]

    if (is.null(offsetX))
    {
        minPos <- 1
        maxPos <- max(width(object)[sel])
        offsetX <- integer(0)
    }
    else
    {
        startPosX <- -offsetX[sel] + 1
        minPos <- min(startPosX)
        maxPos <- max(startPosX + width(object)[sel]) - 1
    }

    if (!is(kernel, "MotifKernel"))
    {
        k <- kernelParameters(kernel)$k
        motifs <- character(0)
        motifLengths <- integer(0)
        maxMotifLength <- 0
        maxPatternLength <- 0
        fwMotifs <- character(0)
        fwMotifLengths <- integer(0)
        fwMaxMotifLength <- 0
        fwMaxPatternLength <- 0
        nodeLimit <- 0
    }
    else
    {
        ## for annotation specific and position dependent kernel
        ## use original motifs otherwise create tree from motifs
        ## in featureWeights

        motifs <- kernelParameters(kernel)$motifs
        motifLengths <- kernelParameters(kernel)$motifLengths

        maxMotifLength <- max(motifLengths)
        maxPatternLength <- max(nchar(motifs))

        fwMotifs <- colnames(featureWeights)
        fwMotifLengths <- nchar(fwMotifs)

        result <- .Call("validateMotifsC", fwMotifs, fwMotifLengths)
        fwMaxMotifLength <- max(fwMotifLengths)
        fwMaxPatternLength <- max(nchar(fwMotifs))

        if (length(kernelParameters(kernel)$distWeight) == 0 &&
            !kernelParameters(kernel)$annSpec)
        {
            ## rough limit for no of nodes in motif tree from no of
            ## chars and no of substitution groups, add one for root
            nodeLimit <- sum(fwMotifLengths) + 1 +
                         sum(sapply(gregexpr("[", fwMotifs, fixed=TRUE),
                                    function(x) length(unlist(x))))
        }
        else
        {
            ## rough limit for no of nodes in motif tree from no of
            ## chars and no of substitution groups, add one for root
            nodeLimit <- sum(motifLengths) + 1 +
                             sum(sapply(gregexpr("[", motifs, fixed=TRUE),
                                        function(x) length(unlist(x))))
        }

        k <- 0
    }

    if (is(kernel, "MismatchKernel") || is(kernel, "GappyPairKernel"))
        m <- kernelParameters(kernel)$m
    else
        m <- 0

    distWeight <- kernelParameters(kernel)$distWeight

    if (length(distWeight) > 0)
    {
        if (is(kernel, "SpectrumKernel"))
            minFeatureLength <- k
        else if (is(kernel, "GappyPairKernel"))
            minFeatureLength <- k
        else if (is(kernel, "MotifKernel"))
            minFeatureLength <- min(motifLengths)

        if (is.list(featureWeights))
        {
            posNames <- lapply(featureWeights, colnames)
            minPosSV <- min(unlist(lapply(posNames, "[", 1)))
            maxPosSV <- max(unlist(lapply(posNames, "[",
            unlist(lapply(posNames, length)))))
        }
        else
        {
            posNames <- colnames(featureWeights)
            minPosSV <- as.numeric(posNames[1])
            maxPosSV <- as.numeric(posNames[length(posNames)])
        }

        maxDist <- max(maxPos, maxPosSV) - min(minPos, minPosSV)

        if (is.function(distWeight))
        {
            ## precompute distance weight vector
            ## terminate on stop and warning
            ## assuming that all distances are partially overlapping
            distWeight <- tryCatch(distWeight(0:(maxDist - minFeatureLength + 1)),
                                   warning=function(w) {stop(w)},
                                   error=function(e) {stop(e)})

            if (!(is.numeric(distWeight) && length(distWeight) ==
                  maxDist - minFeatureLength + 2))
            {
                stop("distWeight function did not return a numeric vector\n",
                     "       of correct length\n")
            }

            ## limit to values larger than .Machine$double.eps
            ## for non-monotonic decreasing functions search from end
            for (i in (maxDist - minFeatureLength + 2):1)
            {
                if (distWeight[i] > .Machine$double.eps)
                    break
            }

            distWeight <- distWeight[1:i]
        }

        if (isTRUE(all.equal(distWeight, c(1, rep(0, length(distWeight)-1)))))
        {
            posSpec <- TRUE
            distWeight <- numeric(0)
        }
        else
            posSpec <- FALSE
    }

    normalized <- kernelParameters(kernel)$normalized
    presence <- kernelParameters(kernel)$presence
    ignoreLower <- kernelParameters(kernel)$ignoreLower
    reverseComplement <- kernelParameters(kernel)$revComplement
    bioCharset <- getBioCharset(object, TRUE)
    unmapped <- is(object, "DNAStringSet") || is(object, "RNAStringSet")
    isXStringSet <- is(object, "XStringSet")
    annot <- NULL
    annCharset <- NULL

    if (annSpec)
    {
        annCharset <- metadata(object)$annotationCharset

        if (is.null(annCharset))
            stop("missing annotation characterset metadata in 'object'\n")

        annot <- mcols(object)[["annotation"]]

        if (is.null(annot))
            stop("missing annotation information in 'object'\n")
    }

    if (length(kernelParameters(kernel)$distWeight) == 0)
    {
        ## set exit hook for C heap cleanup
        on.exit(.Call("freeHeapCallocsC",
                      as.integer(kernelType[class(kernel)])))

        predProf@profiles <- .Call("generatePredictionProfilesC", object,
                as.logical(!isXStringSet), sel - 1, as.integer(length(sel)),
                annCharset, annot, as.integer(maxSeqLength),
                as.logical(unmapped), as.logical(reverseComplement),
                as.integer(kernelType[class(kernel)]),
                as.integer(k), as.integer(m), as.integer(bioCharset[[2]]),
                featureWeights, as.integer(svmIndex - 1), motifs,
                motifLengths, as.integer(maxMotifLength),
                as.integer(maxPatternLength), fwMotifs, fwMotifLengths,
                as.integer(fwMaxMotifLength), as.integer(fwMaxPatternLength),
                as.integer(nodeLimit), as.logical(!ignoreLower),
                as.logical(normalized), as.logical(presence))
    }
    else
    {
        pos1 <- as.numeric(colnames(featureWeights)[1])

        predProf@profiles <- .Call("getPosDepPredOrProfC", featureWeights,
                as.double(weightLimit), b, object, as.logical(isXStringSet),
                sel - 1, as.integer(length(sel)), offsetX,
                as.integer(maxSeqLength), as.integer(bioCharset[[2]]),
                as.integer(kernelType[class(kernel)]), as.integer(k),
                as.integer(m), motifs, motifLengths, as.integer(maxMotifLength),
                as.integer(maxPatternLength), as.integer(nodeLimit),
                as.logical(posSpec), distWeight, as.logical(ignoreLower),
                as.logical(unmapped), as.logical(reverseComplement),
                as.logical(normalized), as.logical(TRUE), as.numeric(pos1),
                as.integer(minPos), as.integer(maxPos))
    }

    if (length(names(object)) > 0)
        rownames(predProf@profiles) <- names(object)[sel]

    colnames(predProf@profiles) <- paste("Pos", seq.int(from=minPos, to=maxPos))

    return(predProf)
}

## Prediction profiles can be generated
## during prediction (see parameter
## \code{predProfiles} in \code{\link{predict}} or specifically
## for a given set of sequences with \code{\link{getPredictionProfile}}.

#' @rdname getPredictionProfile-methods
#' @aliases
#' getPredictionProfile
#'
#' @title Calculation Of Predicition Profiles
#'
#' @description compute prediction profiles for a given set of biological
#' sequences from a model trained with /code{kbsvm}
#'
#' @param object a single biological sequence in the form of an
#' \code{\linkS4class{DNAString}}, \code{\linkS4class{RNAString}} or
#' \code{\linkS4class{AAString}} or multiple biological sequences as
#' \code{\linkS4class{DNAStringSet}}, \code{\linkS4class{RNAStringSet}},
#' \code{\linkS4class{AAStringSet}} (or as \code{\linkS4class{BioVector}}).
#'
#' @param kernel a sequence kernel object of class
#' \code{\linkS4class{SequenceKernel}}.
#'
#' @param featureWeights a feature weights matrix retrieved from a KeBABS model
#' with the accessor \code{\link{featureWeights}}.
#'
#' @param b model intercept from a KeBABS model.
#'
#' @param svmIndex integer value selecting one of the pairwise SVMs in case of
#' pairwise multiclass classification. Default=1
#'
#' @param sel subset of indices into \code{x} as integer vector. When this
#' parameter is present the prediction profiles are computed for the specified
#' subset of samples only. Default=\code{integer(0)}
#'
#' @param weightLimit the feature weight limit is a single numeric value and
#' allows pruning of feature weights. All feature weights with an absolute
#' value below this limit are set to 0 and are not considered for the
#' prediction profile computation. This parameter is only relevant when
#' feature weights are calculated in KeBABS during training.
#' Default=.Machine$double.eps
#'
#' @details
#'
#' With this method prediction profiles can be generated explicitely for a
#' given set of sequences with a given model represented through its feature
#' weights and the model intercept b. A single prediction profile shows for
#' each position of the sequence the contribution of the patterns at this
#' position to the decision value. The prediciion profile also includes the
#' kernel object used for the generation of the profile and the seqence
#' data.\cr\cr
#' A single profile or a pair can be plotted with method \code{\link{plot}}
#' showing the relevance of sequence positions for the prediction. Please
#' consider that patterns occuring at neighboring sequence positions are not
#' statistically independent which means that the relevance of a specific
#' position is not only determined by the patterns at this position but is also
#' influenced by the neighborhood around this position. Prediction profiles can
#' also be generated implicitely during predction for the predicted samples
#' (see parameter \code{predProfiles} in \code{\link{predict}}).\cr\cr
#'
#' @return
#' getPredictionProfile: upon successful completion, the function returns a set
#' of prediction profiles for the sequences as class
#' \code{\linkS4class{PredictionProfile}}.
#'
#' @seealso \code{\linkS4class{PredictionProfile}}, \code{\link{predict}},
#' \code{\link{plot}}, \code{\link{featureWeights}},
#' \code{\link{getPredProfMixture}}
#'
#'
#'
#' @examples
#'
#'
#' ## set random generator seed to make the results of this example
#' ## reproducable
#' set.seed(123)
#'
#' ## load coiled coil data
#' data(CCoil)
#' gappya <- gappyPairKernel(k=1,m=11, annSpec=TRUE)
#' model <- kbsvm(x=ccseq, y=as.numeric(yCC), kernel=gappya, 
#'                pkg="e1071", svm="C-svc", cost=15)
#'
#' ## show feature weights
#' featureWeights(model)[,1:5]
#' 
#' ## define two new sequences to be predicted
#' GCN4 <- AAStringSet(c("MKQLEDKVEELLSKNYHLENEVARLKKLV",
#'                       "MKQLEDKVEELLSKYYHTENEVARLKKLV"))
#' names(GCN4) <- c("GCN4wt", "GCN_N16Y,L19T")
#' ## assign annotation metadata
#' annCharset <- annotationCharset(ccseq)
#' annot <- c("abcdefgabcdefgabcdefgabcdefga",
#'            "abcdefgabcdefgabcdefgabcdefga")
#' annotationMetadata(GCN4, annCharset=annCharset) <- annot
#'
#' ## compute prediction profiles
#' predProf <- getPredictionProfile(GCN4, gappya, 
#'            featureWeights(model), modelOffset(model))
#'
#' ## show prediction profiles
#' predProf
#'
#' ## plot prediction profile of first aa sequence
#' plot(predProf, sel=1, ylim=c(-0.4, 0.2), heptads=TRUE, annotate=TRUE)
#'
#' ## plot prediction profile of both aa sequences
#' plot(predProf, sel=c(1,2), ylim=c(-0.4, 0.2), heptads=TRUE, annotate=TRUE)
#'
#' ## prediction profiles can also be generated during prediction
#' ## when setting the parameter predProf to TRUE
#' ## plotting longer sequences to pdf is shown in the examples for the 
#' ## plot function 
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs/}\cr\cr
#' (Mahrenholz, 2011) -- C.C. Mahrenholz, I.G. Abfalter, U. Bodenhofer,
#' R. Volkmer, and S. Hochreiter. Complex networks govern coiled coil
#' oligomerization - predicting and profiling by means of a machine learning
#' approach.\cr\cr
#' (Bodenhofer, 2009) -- U. Bodenhofer, K. Schwarzbauer, S. Ionescu, and
#' S. Hochreiter. Modeling Position Specificity in Sequence Kernels by
#' Fuzzy Equivalence Relations. \cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics}, 31(15):2574-2576, 2015.
#' DOI: \doi{10.1093/bioinformatics/btv176}.
#' @keywords prediction profile
#' @keywords feature weights
#' @keywords methods
#'

#' @rdname getPredictionProfile-methods
#' @aliases
#' getPredictionProfile,BioVector-method
#' @export

setMethod("getPredictionProfile", signature(object = "BioVector"),
          getPredictionProfile.BioVectorOrXSS)

#' @rdname getPredictionProfile-methods
#' @aliases
#' getPredictionProfile,XStringSet-method
#' @export
#'

setMethod("getPredictionProfile", signature(object = "XStringSet"),
          getPredictionProfile.BioVectorOrXSS)


getPredictionProfile.XString <- function(object, kernel, featureWeights, b,
                       svmIndex=1, sel=NULL, weightLimit=.Machine$double.eps)
{
    return(
        switch(class(object),
            "DNAString" = getPredictionProfile(object=DNAStringSet(object),
                                               kernel=kernel, b=b,
                                               svmIndex=svmIndex,
                                               featureWeights=featureWeights,
                                               weightLimit=weightLimit),
            "RNAString" = getPredictionProfile(object=RNAStringSet(object),
                                               kernel=kernel, b=b,
                                               svmIndex=svmIndex,
                                               featureWeights=featureWeights,
                                               weightLimit=weightLimit),
            "AAString"  = getPredictionProfile(object=AAStringSet(object),
                                               kernel=kernel, b=b,
                                               svmIndex=svmIndex,
                                               featureWeights=featureWeights,
                                               weightLimit=weightLimit),
                  stop("wrong class of x\n")
        )
    )
}

#' @rdname getPredictionProfile-methods
#' @aliases
#' getPredictionProfile,XString-method
#' @export
#'

setMethod("getPredictionProfile", signature(object = "XString"),
          getPredictionProfile.XString)



getPredProfMixture.BioVectorOrXSS <- function(object, trainseqs, mixModel,
                   kernels, mixCoef, svmIndex=1, sel=1:length(object),
                   weightLimit=.Machine$double.eps)
{
    if (missing(object) || is.null(object) || length(object) < 1 ||
        !(class(object) %in% kebabsInfo@allowedSeqSetClasses))
    {
        stop(paste("'object' must be a",
                   paste(kebabsInfo@allowedSeqClasses, collapse=", "), "\n"))
    }

    if (!(class(trainseqs) %in% kebabsInfo@allowedSeqSetClasses))
    {
        stop(paste("'trainseqs' must be a",
                   paste(kebabsInfo@allowedSeqClasses, collapse=", "), "\n"))
    }

    if (!is(mixModel, "KBModel"))
        stop("'mixModel' must be a mixed kernel model of class KBModel\n")

    if (length(kernels) == 1 && is(kernels, "SequenceKernel"))
        kernels <- list(kernels)

    if (!is.list(kernels) || length(kernels) < 1 ||
        any(unlist(lapply(kernels, function(x) !is(x, "SequenceKernel")))))
        stop("'kernels' must be a list of SequenceKernels\n")

    if (!is.numeric(mixCoef) || length(mixCoef) != length(kernels))
        stop("'mixCoef' does not have the same length as 'kernels'\n")

    ## store support vectors in model for model trained with mixed km
    ind <- getSVMSlotValue("svIndex", mixModel)

    if (length(ind) > 0)
    {
        mixModel@SV <- trainseqs[ind]
        mixModel@svIndex <- ind
        mixModel@alphaIndex <- getSVMSlotValue("alphaIndex", mixModel)
    }
    else
        stop("no support vectors in model\n")

    ## compute predprof for each of the kernels and mix
    b <- getSVMSlotValue("b", mixModel)[svmIndex]
    model1 <- mixModel
    model1@svmInfo@selKernel <- kernels[[1]]
    fw1 <- getFeatureWeights(model=model1,
                             weightLimit=model1@svmInfo@weightLimit)

    prof <- getPredictionProfile(object, kernels[[1]], fw1, b,
                svmIndex=svmIndex, sel=sel, weightLimit=weightLimit)

    prof@profiles <- mixCoef[1] * prof@profiles
    prof@kernel <- kernels
    prof@baselines <- - b / width(object[sel])

    if (length(kernels) == 1)
    {
        prof@kernel <- kernels[[1]]
        return(prof)
    }

    for (i in 2:length(kernels))
    {
        model2 <- mixModel
        model2@svmInfo@selKernel <- kernels[[i]]
        fw2 <- getFeatureWeights(model=model2,
                                 weightLimit=model2@svmInfo@weightLimit)

        prof2 <- getPredictionProfile(object, kernels[[i]], fw2, b,
                     svmIndex=svmIndex, sel=sel, weightLimit=weightLimit)

        prof@profiles <- prof@profiles + mixCoef[i] * prof2@profiles
    }

    return(prof)
}

#' @rdname getPredProfMixture-methods
#' @aliases
#' getPredProfMixture
#'
#' @title Calculation Of Predicition Profiles for Mixture Kernels
#'
#' @description compute prediction profiles for a given set of biological
#' sequences from a model trained with mixture kernels
#'
#' @param object a single biological sequence in the form of an
#' \code{\linkS4class{DNAString}}, \code{\linkS4class{RNAString}} or
#' \code{\linkS4class{AAString}} or multiple biological sequences as
#' \code{\linkS4class{DNAStringSet}}, \code{\linkS4class{RNAStringSet}},
#' \code{\linkS4class{AAStringSet}} (or as \code{\linkS4class{BioVector}}).
#'
#' @param trainseqs training sequences on which the mixture model was
#' trained as
#' \code{\linkS4class{DNAStringSet}}, \code{\linkS4class{RNAStringSet}},
#' \code{\linkS4class{AAStringSet}} (or as \code{\linkS4class{BioVector}}).
#'
#' @param mixModel model object of class \code{\linkS4class{KBModel}}
#' trained with kernel mixture.
#'
#' @param kernels a list of sequence kernel objects of class
#' \code{\linkS4class{SequenceKernel}}. The same kernels must be used as in
#' training.
#'
#' @param mixCoef mixing coefficients for the kernel mixture. The same mixing
#' coefficient values must be used as in training.
#'
#' @param svmIndex integer value selecting one of the pairwise SVMs in case of
#' pairwise multiclass classification. Default=1
#'
#' @param sel subset of indices into \code{x} as integer vector. When this
#' parameter is present the prediction profiles are computed for the specified
#' subset of samples only. Default=\code{integer(0)}
#'
#' @param weightLimit the feature weight limit is a single numeric value and
#' allows pruning of feature weights. All feature weights with an absolute
#' value below this limit are set to 0 and are not considered for the
#' prediction profile computation. This parameter is only relevant when
#' feature weights are calculated in KeBABS during training.
#' Default=.Machine$double.eps
#'
#' @details
#'
#' With this method prediction profiles can be generated explicitely for a
#' given set of sequences with a model trained on a precomputed kernel matrix
#' as mixture of multiple kernels.\cr\cr
#'
#' @return
#' upon successful completion, the function returns a set
#' of prediction profiles for the sequences as class
#' \code{\linkS4class{PredictionProfile}}.
#'
#' @seealso \code{\linkS4class{PredictionProfile}}, \code{\link{predict}},
#' \code{\link{plot}}, \code{\link{featureWeights}},
#' \code{\link{getPredictionProfile}}
#'
#'
#'
#' @examples
#'
#'
#' ## set random generator seed to make the results of this example
#' ## reproducable
#' set.seed(123)
#'
#' ## load coiled coil data
#' data(CCoil)
#' gappya1 <- gappyPairKernel(k=1,m=11, annSpec=TRUE)
#' gappya2 <- gappyPairKernel(k=2,m=9, annSpec=TRUE)
#' kernels <- list(gappya1, gappya2)
#' mixCoef <- c(0.7,0.3)
#'
#' ## precompute mixed kernel matrix
#' km <- as.KernelMatrix(mixCoef[1]*gappya1(ccseq) +
#'                       mixCoef[2]*gappya2(ccseq))
#' mixModel <- kbsvm(x=km, y=as.numeric(yCC),
#'                pkg="e1071", svm="C-svc", cost=15)
#'
#' ## define two new sequences to be predicted
#' GCN4 <- AAStringSet(c("MKQLEDKVEELLSKNYHLENEVARLKKLV",
#'                       "MKQLEDKVEELLSKYYHTENEVARLKKLV"))
#' names(GCN4) <- c("GCN4wt", "GCN_N16Y,L19T")
#' ## assign annotation metadata
#' annCharset <- annotationCharset(ccseq)
#' annot <- c("abcdefgabcdefgabcdefgabcdefga",
#'            "abcdefgabcdefgabcdefgabcdefga")
#' annotationMetadata(GCN4, annCharset=annCharset) <- annot
#'
#' ## compute prediction profiles
#' predProf <- getPredProfMixture(GCN4, ccseq, mixModel,
#'                                kernels, mixCoef)
#'
#' ## show prediction profiles
#' predProf
#'
#' ## plot prediction profile of both aa sequences
#' plot(predProf, sel=c(1,2), ylim=c(-0.4, 0.2), heptads=TRUE, annotate=TRUE)
#'
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs/}\cr\cr
#' (Mahrenholz, 2011) -- C.C. Mahrenholz, I.G. Abfalter, U. Bodenhofer,
#' R. Volkmer, and S. Hochreiter. Complex networks govern coiled coil
#' oligomerization - predicting and profiling by means of a machine learning
#' approach.\cr\cr
#' (Bodenhofer, 2009) -- U. Bodenhofer, K. Schwarzbauer, S. Ionescu, and
#' S. Hochreiter. Modeling Position Specificity in Sequence Kernels by
#' Fuzzy Equivalence Relations. \cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics}, 31(15):2574-2576, 2015.
#' DOI: \doi{10.1093/bioinformatics/btv176}.
#' @keywords prediction profile
#' @keywords feature weights
#' @keywords methods
#'

#' @rdname getPredProfMixture-methods
#' @aliases
#' getPredProfMixture,BioVector-method
#' @export

setMethod("getPredProfMixture", signature(object = "BioVector"),
          getPredProfMixture.BioVectorOrXSS)

#' @rdname getPredProfMixture-methods
#' @aliases
#' getPredProfMixture,XStringSet-method
#' @export
#'

setMethod("getPredProfMixture", signature(object = "XStringSet"),
          getPredProfMixture.BioVectorOrXSS)


getPredProfMixture.XString <- function(object, trainseqs, mixModel,
                                  kernels, mixCoef, svmIndex=1, sel=1,
                                  weightLimit=.Machine$double.eps)
{
    return(
        switch(class(object),
            "DNAString" = getPredProfMixture(object=DNAStringSet(object),
                                trainseqs=trainseqs, kernels=kernels,
                                mixModel=mixModel, mixCoef=mixCoef,
                                svmIndex=svmIndex, sel=sel,
                                weightLimit=weightLimit),
            "RNAString" = getPredProfMixture(object=RNAStringSet(object),
                                trainseqs=trainseqs, kernels=kernels,
                                mixModel=mixModel, mixCoef=mixCoef,
                                svmIndex=svmIndex, sel=sel,
                                weightLimit=weightLimit),
            "AAString" = getPredProfMixture(object=AAStringSet(object),
                                trainseqs=trainseqs, kernels=kernels,
                                mixModel=mixModel, mixCoef=mixCoef,
                                svmIndex=svmIndex, sel=sel,
                                weightLimit=weightLimit),
            stop("wrong class of x\n")
        )
    )
}

#' @rdname getPredProfMixture-methods
#' @aliases
#' getPredProfMixture,XString-method
#' @export
#'

setMethod("getPredProfMixture", signature(object = "XString"),
          getPredProfMixture.XString)

