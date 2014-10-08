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
        maxPos <- max(width(object[sel]))
        offsetX <- integer(0)
    }
    else
    {
        minPos <- min(offsetX)
        maxPos <- max(offsetX + width(object[sel]) - 1)
    }

    if (!is(kernel, "MotifKernel"))
    {
        k <- kernelParameters(kernel)$k
        motifs <- NULL
        motifLengths <- NULL
        maxMotifLength <- 0
        maxPatternLength <- 0
        fwMotifs <- NULL
        fwMotifLengths <- NULL
        fwMaxMotifLength <- 0
        fwMaxPatternLength <- 0
        nodeLimit <- 0
    }
    else
    {
        ## for annotation specific kernel use original motifs
        ## otherwise create tree from motifs in featureWeights

        motifs <- kernelParameters(kernel)$motifs
        motifLengths <- kernelParameters(kernel)$motifLengths

        fwMotifs <- colnames(featureWeights)
        fwMotifLengths <- nchar(motifs)

        result <- .Call("validateMotifsC", motifs, motifLengths)
        maxMotifLength <- max(motifLengths)
        maxPatternLength <- max(nchar(motifs))

        result <- .Call("validateMotifsC", fwMotifs, fwMotifLengths)
        fwMaxMotifLength <- max(fwMotifLengths)
        fwMaxPatternLength <- max(nchar(fwMotifs))

        ## rough limit for no of nodes in motif tree from no of
        ## chars and no of substitution groups, add one for root
        nodeLimit <- sum(motifLengths) + 1 +
        sum(sapply(gregexpr("[", motifs, fixed=TRUE),
                   function(x) length(unlist(x))))
        k <- 0
    }


    if (is(kernel, "MismatchKernel") || is(kernel, "GappyPairKernel"))
        m <- kernelParameters(kernel)$m
    else
        m <- 0

    distWeight <- kernelParameters(kernel)$distWeight

    if (length(distWeight) > 0)
    {
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
#' featureWeights are calculated in KeBABS during training.
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
#' \code{\link{plot}}, \code{\link{featureWeights}{KBModel}}
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
#' (Mahrenholz, 2011) -- Mahrenholz, C.C., Abfalter, I.G., Bodenhofer, U.,
#' Volkmer, R., and Hochreiter, S.; Complex networks govern coiled coil
#' oligomerization - predicting and profiling by means of a machine learning
#' approach.
#'
#' (Bodenhofer, 2009) -- Bodenhofer, U., Schwarzbauer, K.,  Ionescu, S. and
#' Hochreiter,S., Modeling Position Specificity in Sequence Kernels by
#' Fuzzy Equivalence Relations. \cr\cr
#' \url{http://www.bioinf.jku.at/software/kebabs}
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
                                               featureWeights=featureWeights,
                                               weightLimit=weightLimit),
            "RNAString" = getPredictionProfile(object=RNAStringSet(object),
                                               kernel=kernel, b=b,
                                               featureWeights=featureWeights,
                                               weightLimit=weightLimit),
            "AAString"  = getPredictionProfile(object=AAStringSet(object),
                                               kernel=kernel, b=b,
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
