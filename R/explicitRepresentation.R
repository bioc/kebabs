##345678901234567890123456789012345678901234567890123456789012345678901234567890
#' @rdname explicitRepresentation
#' @title Explict Representation
#'
#' @description Create an explicit representation
#'
#' @param x one or multiple biological sequences in the form of a
#' \code{\linkS4class{DNAStringSet}}, \code{\linkS4class{RNAStringSet}},
#' \code{\linkS4class{AAStringSet}} (or as \code{\linkS4class{BioVector}})
#'
#' @param selx subset of indices into \code{x}. When this parameter is present
#' the explicit representation is generated for the specified subset of samples
#' only. default=\code{NULL}
#'
#' @param kernel a sequence kernel object. The feature map of this kernel
#' object is used to generate the explicit representation.
#'
#' @param sparse boolean that indicates whether a sparse or dense explicit
#' representation should be generated. Default=TRUE
#'
#' @param zeroFeatures indicates whether columns with zero feature counts
#' across all samples should be included in the explicit representation.
#' (see below) Default=FALSE
#'
#' @param features feature subset of the specified kernel in the form of a
#' character vector. When a feature subset is passed to the function all other
#' features in the feature space are not considered for the explicit
#' representation. (see below)
#'
#' @param useRowNames if this parameter is set the sample names will be set
#' as row names if available in the provided sequence set.
#' Default=TRUE
#'
#' @param useColNames if this parameter is set the features will be set
#' as column names in the explicit representation.
#' Default=TRUE
#'
#' @details
#' Creation of an explicit representation\cr\cr
#' The function 'getExRep' creates an explicit representation of the given
#' sequence set using the feature map of the specified kernel. It contains
#' the feature counts in a matrix format. The rows of the matrix represent
#' the samples, the columns the features. For a dense explicit representation
#' of class \code{\linkS4class{ExplicitRepresentationDense}} the count data
#' is stored in a dense matrix. To allow efficient storage all features that
#' do not occur in the sequence set are removed from the explicit
#' representation by default. When the parameter \code{zeroFeatures} is set
#' to \code{TRUE} these features are also included resulting an explicit
#' representation which contains the full feature space. For feature spaces
#' larger than one million features the inclusion of zero features is not
#' possible. \cr\cr
#' In case of large feature spaces a sparse explicit representation of class
#' \code{\linkS4class{ExplicitRepresentationSparse}} is much more efficient
#' by storing the count data as \code{dgRMatrix} from package \bold{Matrix}).
#' The class \code{\linkS4class{ExplicitRepresentationSparse}}
#' is derived from \code{dgRMatrix}. As zero features are not stored in a
#' sparse matrix the flag \code{zeroFeatures} only controls whether the
#' column names of features not occuring in the sequences are included or
#' not.\cr\cr
#' Both the dense and the sparse explicit representation also contain the
#' kernel object which was used for it's creation. For an explicit
#' representation without zero features column names are mandatory.
#' An explicit representation can be created for position independent and
#' annotation specific kernel variants (for details see
#' \link{annotationMetadata}). In annotation specific kernels the
#' annotation characters are included as postfix in the features. For kernels
#' with normalization the explicit representation is normalized resulting in
#' row vectors normalized to the unit sphere. For feature subsets used with
#' normalized kernels all features of the feature space are used in the
#' normalization.\cr\cr
#' Usage of explicit representations\cr\cr
#' Learning with linear SVMs (e.g. \code{\link[kernlab]{ksvm}}in package
#' \bold{kernlab} or \code{\link[e1071]{svm}} in package \bold{e1071}) can be
#' performed either through passing a kernel matrix of similarity values or
#' an explicit representation and a linear kernel to the SVM. The SVMs in
#' package \bold{kernlab} support a dense explicit representation or kernel
#' matrix as data representations. The SVMs in packages \bold{e1071}) and
#' \bold{LiblineaR} support dense or sparse explicit representations.
#' In many cases there can be considerable performance differences between
#' the two variants of passing data to the SVM. And especially for larger
#' feature spaces the sparse explicit representation not only brings higher
#' memory efficiency but also leads to drastically improved runtimes during
#' training and prediction. Starting with kebabs version 1.2.0 kernel matrix
#' support is also available for package \code{e1071} via the dense LIBSVM
#' implementation integrated in package \code{kebabs}.\cr\cr
#' In general all of the complexity of converting the sequences with a specific
#' kernel to an explicit representation or a kernel matrix and adapting the
#' formats and parameters to the specific SVM is hidden within the KeBABS
#' training and predict methods (see \code{\link{kbsvm}},
#' \code{\link{predict}}) and the user can concentrate on the actual data
#' analysis task. During training via \code{\link{kbsvm}} the parameter
#' \code{explicit} controls the training via kernel matrix or explicit
#' representation and the parameter \code{explicitType} determines whether a
#' dense or sparse explicit representation is used. Manual generation of
#' explicit representations is only necessary for usage with other learners
#' or analysis methods not supported by KeBABS.\cr\cr
#' Quadratic explicit representation\cr\cr
#' The package \bold{LiblineaR} only provides linear SVMs which are tuned for
#' efficient processing of larger feature spaces and sample numbers. To allow
#' the use of a quadratic kernel on these SVMs a quadratic explicit
#' representation can be generated from the linear explicit representation.
#' It contains counts for feature pairs and the features combined to one pair
#' are separated by '_' in the column names of the quadratic explicit
#' representation. Please be aware that the dimensionality for a quadratic
#' explicit representation increases considerably compared to the linear one.
#' In the other SVMs a linear explicit representation together with a quadratic
#' kernel is used instead. In training via \code{\link{kbsvm}} the use
#' of a linear representation with a quadratic kernel or a quadratic explicit
#' representation instead is indicated through setting the parameter
#' \code{featureType} to the value \code{"quadratic"}.
#' @return
#' getExRep: upon successful completion, dependent on the flag \code{sparse}
#' the function returns either a dense explicit representation of class
#' \code{\linkS4class{ExplicitRepresentationDense}} or a sparse explicit
#' representation of class \code{\linkS4class{ExplicitRepresentationSparse}}.
#'
#' @seealso \code{\linkS4class{ExplicitRepresentationDense}},
#' \code{\linkS4class{ExplicitRepresentationSparse}},
#' \code{\link{getKernelMatrix}}, \code{\link{kernelParameters-method}},
#' \code{\linkS4class{SpectrumKernel}}, \code{\link{mismatchKernel}},
#' \code{\link{gappyPairKernel}}, \code{\link{motifKernel}}
#'
#' @examples
#'
#' ## instead of user provided sequences in XStringSet format
#' ## for this example a set of DNA sequences is created
#' ## RNA- or AA-sequences can be used as well with the spectrum kernel
#' dnaseqs <- DNAStringSet(c("AGACTTAAGGGACCTGGACACCACACTCAGCTAGGGGGACTGGGAGC",
#'                           "ATAAAGGGAGCAGACATCATGACCTTTTTGACCCTAATTATTTCAGC",
#'                           "CAGGAATCAGCACAGGCAGGGGCACTGCATCCCAAGACATCTGGGCC",
#'                           "GGACATATACCCACCCTTACCTGCCATACAGGATAGGGCCACTGCCC",
#'                           "ATAAAGGATGCAGACATCATGGCCTTTTTGACCCTAATTATTTCAGC"))
#' names(dnaseqs) <- paste("S", 1:length(dnaseqs), sep="")
#'
#' ## create the kernel object for dimers with normalization
#' speck <- spectrumKernel(k=2)
#' ## show details of kernel object
#' speck
#'
#' ## generate the dense explicit representation for the kernel
#' erd <- getExRep(dnaseqs, speck, sparse=FALSE)
#' dim(erd)
#' erd[1:5,]
#'
#' ## generate the dense explicit representation with zero features
#' erd <- getExRep(dnaseqs, speck, sparse=FALSE, zeroFeatures=TRUE)
#' dim(erd)
#' erd[1:5,]
#'
#' ## generate the sparse explicit representation for the kernel
#' ers <- getExRep(dnaseqs, speck)
#' dim(ers)
#' ers[1:5,]
#'
#' ## generate the sparse explicit representation with zero features
#' ers <- getExRep(dnaseqs, speck, zeroFeatures=TRUE)
#' dim(ers)
#' ers[1:5,]
#'
#' ## generate the quadratic explicit representation
#' erdq <- getExRepQuadratic(erd)
#' dim(erdq)
#' erdq[1:5,1:15]
#'
#' \dontrun{
#' ## run taining and prediction with dense linear explicit representation
#' data(TFBS)
#' enhancerFB
#' train <- sample(1:length(enhancerFB), length(enhancerFB) * 0.7)
#' test <- c(1:length(enhancerFB))[-train]
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=speck,
#'                pkg="LiblineaR", svm="C-svc", cost=10, explicit="yes",
#'                explicitType="dense")
#' pred <- predict(model, x=enhancerFB[test])
#' evaluatePrediction(pred, yFB[test], allLabels=unique(yFB))
#'
#' ## run taining and prediction with sparse linear explicit representation
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=speck,
#'                pkg="LiblineaR", svm="C-svc", cost=10, explicit="yes",
#'                explicitType="sparse")
#' pred <- predict(model, x=enhancerFB[test])
#' evaluatePrediction(pred, yFB[test], allLabels=unique(yFB))
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs/}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics}, 31(15):2574-2576, 2015.
#' DOI: \doi{10.1093/bioinformatics/btv176}.
#' @keywords explicit representation
#' @keywords methods
#' @export
#'

getExRep <- function(x, kernel=spectrumKernel(), sparse=TRUE,
                     zeroFeatures=FALSE, features=NULL, useRowNames=TRUE,
                     useColNames=TRUE, selx=NULL)
{
    if (missing(x) || is.null(x))
    {
        stop(paste("'x' must be a",
                   paste(kebabsInfo@allowedSeqClasses, collapse=", "), "\n"))
    }

    if (class(x) %in% c("DNAString", "RNAString", "AAString"))
    {
        ## avoid copying of x for non XString
        x <- switch(class(x),
                    "DNAString" = DNAStringSet(x),
                    "RNAString" = RNAStringSet(x),
                    "AAString"  = AAStringSet(x)
                    )
    }

    if (!(class(x) %in% kebabsInfo@allowedSeqSetClasses))
    {
        stop(paste("'x' must be a",
                   paste(kebabsInfo@allowedSeqClasses, collapse=", ")))
    }

    if (length(selx) > 0)
    {
        if (!is.numeric(selx) || length(selx) > length(x))
            stop("'selx' must be a numeric vector with indices into 'x'\n")

        selx <- as.integer(selx)
    }
    else
        selx <- 1:length(x)

    if (!isTRUEorFALSE(sparse))
        stop("'sparse' must be TRUE or FALSE\n")

    if (!isTRUEorFALSE(useRowNames))
        stop("'useRowNames' must be TRUE or FALSE\n")

    if (!isTRUEorFALSE(useColNames))
        stop("'useColNames' must be TRUE or FALSE\n")

    if (!isTRUEorFALSE(zeroFeatures))
        stop("'zeroFeatures' must be TRUE or FALSE\n")

    if (!useColNames && zeroFeatures == FALSE)
        stop("column names are required when zero features are removed\n")

    if (!is(kernel, "SequenceKernel"))
        stop("'kernel' must be of class SequenceKernel\n")

    if (isUserDefined(kernel))
    {
        stop("explicit representation not supported for user-defined ",
             "kernel\n")
    }

    if (is(kernel, "SymmetricPairKernel"))
    {
        stop("explicit representation for symmetric pair kernel is\n",
             "        not supported\n")
    }

    if (!is(kernel, "MotifKernel"))
    {
        k <- kernelParameters(kernel)$k
        motifs <- NULL
        motifLengths <- NULL
        maxMotifLength <- 0
        maxPatternLength <- 0
        nodeLimit <- 0

        # limit k to 64 bit feature space
        if (class(x) %in% c("AAVector", "AAStringSet"))
        {
            if (k > 14)
                stop("'k' must be smaller than or equal to 14\n")
        }
        else
        {
            ## only exact charset supported for explicit representation
            if (k > 31)
                stop("for exact charset 'k' must be smaller than or ",
                     "equal to 31\n")
        }
    }
    else
    {
        if (!is.null(features))
            stop("motif kernel does not support feature subsets\n")

        motifs <- kernelParameters(kernel)$motifs
        motifLengths <- kernelParameters(kernel)$motifLengths
        maxMotifLength <- max(motifLengths)
        maxPatternLength <- max(nchar(motifs))

        ## rough limit for no of nodes in motif tree from no of
        ## chars and no of substitution groups, add one for root
        nodeLimit <- sum(motifLengths) + 1 +
                     sum(sapply(gregexpr("[", motifs, fixed=TRUE),
                                function(x) length(unlist(x))))

        k <- 0
    }

    if (length(kernelParameters(kernel)$distWeight) > 0)
    {
        stop("Explicit representation for position specific or \n",
             "      distance weighted kernel not supported\n")
    }

    r <- kernelParameters(kernel)$r

    if (is(kernel, "MismatchKernel") || is(kernel, "GappyPairKernel"))
        m <- kernelParameters(kernel)$m
    else
        m <- 0

    if (r != 1)
        stop("only linear explicit representation supported\n")

    if (is(kernel, "SpectrumKernel") || is(kernel, "GappyPairKernel"))
        mixCoef <- kernelParameters(kernel)$mixCoef
    else
        mixCoef <- numeric(0)

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

    normalized <- kernelParameters(kernel)$normalized
    ignoreLower <- kernelParameters(kernel)$ignoreLower
    presence <- kernelParameters(kernel)$presence
    reverseComplement <- kernelParameters(kernel)$revComplement

    if (is(kernel, "GappyPairKernel"))
        minSeqLength <- 2 * k + m
    else if (is(kernel, "MotifKernel"))
        minSeqLength <- 2 * k + m
    else
        minSeqLength <- maxMotifLength

    if (!is.null(features) && !(is(kernel, "MotifKernel")))
    {
        if (!is.character(features))
            stop("'features' must be a character vector\n")

        if (is(kernel, "GappyPairKernel"))
        {
            minKmerLength <- 2*k
            maxKmerLength <- minKmerLength + m

            if (annSpec)
            {
                minKmerLength <- 2 * minKmerLength
                maxKmerLength <- 2 * maxKmerLength
            }

            if (any(nchar(features) < minKmerLength) ||
                any(nchar(features) > maxKmerLength))
                stop("length of given features is not valid for kernel\n")
        }
        else
        {
            kmerlength <- k

            if (annSpec)
                kmerlength <- 2 * kmerlength

            if ((any(nchar(features) < kmerlength)) ||
                (any(nchar(features) > kmerlength+m)))
                stop("length of given features is not valid for kernel\n")
        }
    }

    bioCharset <- getBioCharset(x, TRUE)

    seqLength <- width(x)[selx]

    ## $$$ TODO limit k

    maxSeqLength <- max(seqLength)

    annX <- NULL
    annCharset <- NULL

    if (annSpec)
    {
        if (is(kernel, "MotifKernel") && zeroFeatures == TRUE)
        {
            stop("motif kernel with annotation does not support explicit \n",
                 "        representation with zero features\n")
        }

        annCharset <- metadata(x)$annotationCharset

        if (is.null(annCharset))
            stop("missing annotation characterset metadata in 'x'\n")

        annX <- mcols(x)[["annotation"]]

        if (is.null(annX))
            stop("missing annotation information in 'x'\n")
    }

    ## set exit hook for C heap cleanup
    on.exit(.Call("freeHeapCallocsC", as.integer(kernelType[class(kernel)])))

    unmapped <- is(x, "DNAStringSet") || is(x, "RNAStringSet")
    isXStringSet <- is(x, "XStringSet")

    if (sparse)
    {
        expRep <- .Call("genExplRepC", x, as.logical(isXStringSet), (selx - 1L),
                        annCharset, annX, as.integer(maxSeqLength),
                        as.integer(kernelType[class(kernel)]),
                        as.integer(k), as.integer(m),
                        as.integer(bioCharset[[2]]), features, motifs,
                        motifLengths, as.integer(maxMotifLength),
                        as.integer(maxPatternLength), as.integer(nodeLimit),
                        as.logical(presence), as.logical(reverseComplement),
                        as.logical(normalized), as.logical(!ignoreLower),
                        as.logical(unmapped), as.logical(useRowNames),
                        as.logical(useColNames), as.logical(zeroFeatures),
                        as.logical(sparse))

        expRep@usedKernel <- kernel
        expRep@quadratic <- FALSE

        ## $$$ TODO implement mixed spectrum sparse er on C level

    }
    else
    {
        expRep <- new("ExplicitRepresentationDense")
        expRep@usedKernel <- kernel
        expRep@quadratic <- FALSE

        if (length(mixCoef) > 0)
        {
            currK <- 0
            for (i in 1:k)
            {
                currK <- currK + 1

                if (mixCoef[i] > 0)
                    break
            }
        }
        else
            currK <- k

        expRep@.Data <- .Call("genExplRepC", x, as.logical(isXStringSet),
                        (selx - 1L), annCharset, annX, as.integer(maxSeqLength),
                        as.integer(kernelType[class(kernel)]),
                        as.integer(currK), as.integer(m),
                        as.integer(bioCharset[[2]]), features, motifs,
                        motifLengths, as.integer(maxMotifLength),
                        as.integer(maxPatternLength), as.integer(nodeLimit),
                        as.logical(presence), as.logical(reverseComplement),
                        as.logical(normalized), as.logical(!ignoreLower),
                        as.logical(unmapped), as.logical(useRowNames),
                        as.logical(useColNames), as.logical(zeroFeatures),
                        as.logical(sparse))

        ## $$$ TODO realize mixed spectrum er dense on C level
        if (length(mixCoef) > 0)
            expRep@.Data <- sqrt(mixCoef[currK]) * expRep@.Data

        if (length(mixCoef) > 0 && currK < k)
        {
            for (i in (currK+1):k)
            {
                if (mixCoef[i] > 0)
                {
                    expRep@.Data <- cbind(expRep@.Data, sqrt(mixCoef[i]) *
                        .Call("genExplRepC", x, as.logical(isXStringSet),
                        (selx - 1L), annCharset, annX, as.integer(maxSeqLength),
                        as.integer(kernelType[class(kernel)]),
                        as.integer(i), as.integer(m),
                        as.integer(bioCharset[[2]]), features, motifs,
                        motifLengths, as.integer(maxMotifLength),
                        as.integer(maxPatternLength), as.integer(nodeLimit),
                        as.logical(presence), as.logical(reverseComplement),
                        as.logical(normalized), as.logical(!ignoreLower),
                        as.logical(unmapped), as.logical(useRowNames),
                        as.logical(useColNames), as.logical(zeroFeatures),
                        as.logical(sparse)))
                }
            }
        }
    }

    if (useRowNames && length(names(x)) > 0 && length(selx) == nrow(expRep))
        rownames(expRep) <- names(x)[selx]

    return(expRep)
}

getSingleFeaturesFromQuadratic <- function(quadFeatures)
{
    if (!is.character(quadFeatures))
        stop("'quadFeatures' must be a character vector")

    if (length(grep("_", quadFeatures)) != length(quadFeatures))
        stop("'quadFeatures' must contain feature pairs separated by '_'")

    return(unique(unlist(strsplit(quadFeatures,"_"))))
}

mapFeatureVector <- function(x, v1, v2, nofea)
{
    factor <- rep(sqrt(2), nofea*(nofea+1)/2)
    diagIndex <- which(sequence(nofea:1) == 1)
    factor[diagIndex] <- 1
    res <- factor * x[v1] * x[v2]
}

#' @rdname explicitRepresentation
#'
#' @param exRepLin a linear explicit representation
#'
#' @return
#' getExRepQuadratic: upon successful completion, the function returns a
#' quadratic explicit representation
#' @export


getExRepQuadratic <- function(exRepLin, useRowNames=TRUE,
                              useColNames=TRUE, zeroFeatures=FALSE)
{
    if (is(exRepLin, "ExplicitRepresentationSparse"))
    {
        return(getExRepQuadraticSparse(exRepLin, useRowNames,
                                       useColNames, zeroFeatures))
    }

    if (!is(exRepLin, "ExplicitRepresentationDense"))
        stop("'exRepLin' must be of class ExplicitRepresentationDense\n")

    if (exRepLin@quadratic)
        stop("function called with quadratic explicit representation\n")

    nofea <- ncol(exRepLin)
    v1 <- c(unlist(mapply(rep, 1:nofea, nofea:1)))
    v2 <- c(unlist(sapply(1:nofea,
                                   function(x){x:(nofea)})))

    erq <- t(apply(exRepLin, 1,
                   function(x) mapFeatureVector(x=x, v1=v1, v2=v2,
                                                nofea=nofea)))

    ## assign names as pairs
    if (useRowNames && length(rownames(exRepLin)) > 0)
        rownames(erq) <- rownames(exRepLin)
    else
        rownames(erq) <- NULL

    if (useColNames && length(colnames(exRepLin)) > 0)
    {
        linNames <- colnames(exRepLin)
        colnames(erq) <- paste(linNames[v1], linNames[v2], sep="_")
    }
    else
        colnames(erq) <- NULL

    ## remove zero features
#    if (!zeroFeatures)
#        erq <- erq[, which(colSums(erq) != 0)]

    exRepQuad <- new("ExplicitRepresentationDense")
    exRepQuad@quadratic <- TRUE
    exRepQuad@usedKernel <- exRepLin@usedKernel
    exRepQuad@.Data <- erq

    return(exRepQuad)
}

## $$$ TODO the direct solution without dense matrix has a strange memory
##          behavior always varying real mem by more than 1GB in a few second
##          rythm and runtime is unbelievable high, interrupted after 1100 sec
##          looks like this should be done in C

#getExRepQuadraticSparse <- function(exRepLin, useRowNames=TRUE,
#                                    useColNames=TRUE, zeroFeatures=FALSE)
#{
#    if (!is(exRepLin, "ExplicitRepresentationSparse"))
#        stop("'exRepLin' must be of class ExplicitRepresentationSparse\n")
#
#    ## allocate first for better gc
#    exRepQuad <- new("ExplicitRepresentationSparse")
#    exRepQuad@quadratic <- TRUE
#    exRepQuad@usedKernel <- exRepLin@usedKernel
#    exRepQuad@factors <- list()
#
#    nofea <- ncol(exRepLin)
#    v1 <- c(unlist(mapply(rep, 1:nofea, nofea:1)))
#    v2 <- c(unlist(sapply(1:nofea,
#                          function(x){x:(nofea)})))
#
#    ## assemble erqs row by row and avoid heap fragmentation
#   ## in a first run count the number of non-zero elements
#    noOfNonZeroElements <- 0
#
#    for (i in 1:nrow(exRepLin))
#    {
#        noOfNonZeroElements <- noOfNonZeroElements +
#        length(which(mapFeatureVector(x=as.numeric(exRepLin[i,]), v1=v1,
#                                      v2=v2, nofea=nofea) > 0))
#    }
#
#    exRepQuad@p <- rep(0L, nrow(exRepLin) + 1)
#    exRepQuad@j <- rep(0L, noOfNonZeroElements)
#    exRepQuad@x <- rep(0L, noOfNonZeroElements)
#    exRepQuad@Dim <- as.integer(c(exRepLin@Dim[[1]], nofea*(nofea+1)/2))
#
#    for (i in 1:nrow(exRepLin))
#    {
#        qVector <- mapFeatureVector(x=as.numeric(exRepLin[i,]), v1=v1,
#                                    v2=v2, nofea=nofea)
#        jq <- which(qVector > 0)
#        exRepQuad@p[i+1] <- exRepQuad@p[i] + length(jq)
#        exRepQuad@j[(exRepQuad@p[i]+1):exRepQuad@p[i+1]] <- jq
#        exRepQuad@x[(exRepQuad@p[i]+1):exRepQuad@p[i+1]] <- qVector[jq]
#    }
#
#    ## assign names as pairs
#    if (useRowNames && length(rownames(exRepLin)) > 0)
#        rowNames <- rownames(exRepLin)
#    else
#        rowNames <- NULL
#
#    if (useColNames && length(colnames(exRepLin)) > 0)
#    {
#        linNames <- colnames(exRepLin)
#        colNames <- paste(linNames[v1], linNames[v2], sep="_")
#    }
#    else
#        colNames <- NULL
#
#    exRepQuad@Dimnames <- list(rowNames, colNames)
#
#    ## remove zero features
##    if (!zeroFeatures)
##        erq <- erq[, which(colSums(erq) != 0)]
#
#    return(exRepQuad)
#}

## $$$ TODO optimize to generate quadratic ER without going through dense
##          this solution only works up to k=5 with reasonable memory
##          consumption - with 1000 samples a 10000 bases, k=5: 1.4GB real mem
##          final erqs is only 400 MB, input erls is 3MB

getExRepQuadraticSparse <- function(exRepLin, useRowNames=TRUE,
                                    useColNames=TRUE, zeroFeatures=FALSE)
{
    if (!is(exRepLin, "ExplicitRepresentationSparse"))
        stop("'exRepLin' must be of class ExplicitRepresentationSparse\n")

    ## allocate first for better gc
    expRepQuad <- new("ExplicitRepresentationSparse")
    expRepQuad@quadratic <- TRUE
    expRepQuad@usedKernel <- exRepLin@usedKernel
    expRepQuad@factors <- list()

    nofea <- ncol(exRepLin)
    v1 <- c(unlist(mapply(rep, 1:nofea, nofea:1)))
    v2 <- c(unlist(sapply(1:nofea,
                          function(x){x:(nofea)})))

    erq <- t(apply(exRepLin, 1,
                   function(x) mapFeatureVector(x=x, v1=v1, v2=v2,
                                                nofea=nofea)))

    ## assign names as pairs
    if (useRowNames && length(rownames(exRepLin)) > 0)
        rownames(erq) <- rownames(exRepLin)
    else
        rownames(erq) <- NULL

    if (useColNames && length(colnames(exRepLin)) > 0)
    {
        linNames <- colnames(exRepLin)
        colnames(erq) <- paste(linNames[v1], linNames[v2], sep="_")
    }
    else
        colnames(erq) <- NULL

    ## remove zero features
    if (!zeroFeatures)
        erq <- erq[, which(colSums(erq) != 0)]

    ersq <- as(erq, "dgRMatrix")

    expRepQuad@Dimnames <- list(exRepLin@Dimnames[[1]], colnames(ersq))
    expRepQuad@Dim <- c(exRepLin@Dim[[1]], ersq@Dim[[2]])
    expRepQuad@p <- ersq@p
    expRepQuad@j <- ersq@j
    expRepQuad@x <- ersq@x

    rm(erq)
    rm(ersq)
    gc()

    return(expRepQuad)
}
