#2345678901234567890123456789012345678901234567890123456789012345678901234567890

isTRUEorFALSE <- function(x)
{
    !is.null(x) && !is.na(x) && is.logical(x) && length(x) == 1L
}

isSingleInteger <- function(x)
{
    !is.null(x) && !is.na(x) && is.integer(x) && length(x) == 1L
}

isSingleNumber <- function(x)
{
    !is.null(x) && !is.na(x) && is.numeric(x) && length(x) == 1L
}

isSingleString <- function(x)
{
    !is.null(x) && !is.na(x) && is.character(x) && length(x) == 1L
}

isPkgInstalled <- function(pkg)
{
    is.element(pkg, installed.packages()[,1])
}

#' @rdname kebabsCollectInfo
#' @title Collect KeBABS Package Information
#' @aliases
#' kebabsCollectInfo
#'
#' @description Collects and prints general R and package version information.
#' If you have a question related to some KeBABS functionality or observe
#' some unexpected behavior please call this function and send its output
#' together with your information to the contact address specified on the
#' title page of the package vignette. The function by default only outputs
#' the package version of those packages which are directly related to the
#' KeBABS functionality.
#'
#' @param onlyKebabsRelated if set to \code{TRUE} only the packages related
#' to KeBABS are shown, if set to \code{FALSE} all attached packages and
#' all packages loaded via namespace are shown. Default=TRUE
#' @return
#' see above
#' @examples
#'
#' ## collect KeBABS related package information
#' kebabsCollectInfo()
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @export

kebabsCollectInfo <- function(onlyKebabsRelated=TRUE)
{
    sInfo <- sessionInfo()

    if (onlyKebabsRelated == TRUE)
    {
        attachedSubset <- intersect(names(sInfo$otherPkgs),
                          c("kebabs", "Biostrings", "XVector", "IRanges",
                            "S4Vectors", "BiocGenerics", "kernlab", "SparseM"))
        ldViaNamespace <- intersect(names(sInfo$loadedOnly),
                          c("Rcpp","Matrix","e1071","LiblineaR","lattice"))
        
        sInfo$otherPkgs <- sInfo$otherPkgs[attachedSubset]
        sInfo$loadedOnly <- sInfo$loadedOnly[ldViaNamespace]
    }

    print(sInfo)

    if (!isPkgInstalled("SparseM"))
    {
        cat("\nnot installed:\n")
        cat("SparseM\n")
    }
    else
    {
        if(!(("SparseM" %in% names(sInfo$otherPkgs)) ||
             ("SparseM" %in% names(sInfo$loadedOnly))))
        {
            cat("\nnot loaded:\n")
            versionSparseM <- as.character(packageVersion("SparseM"))
            cat("SparseM_", c(versionSparseM), "\n", sep="")
        }
    }
    cat("\n")
}

#' @rdname LinearKernel
#' @title Linear Kernel
#' @aliases
#' linearKernel
#'
#' @description Create a dense or sparse kernel matrix from an explicit
#' representation
#'
#' @param x a dense or sparse explicit representation. \code{x} must be a
#' sparse explicit representation when a sparse kernel matrix should be
#' returned by the function (see parameter \code{sparse}).
#' @param y a dense or sparse explicit representation. If \code{x} is dense,
#' \code{y} must be dense. If \code{x} is sparse, \code{y} must be sparse.
#' @param selx a numeric or character vector for defining a subset of
#' \code{x}. Default=integer(0)
#' @param sely a numeric or character vector for defining a subset of
#' \code{y}. Default=integer(0)
#' @param sparse boolean indicating whether returned kernel matrix
#' should be sparse or dense. For value \code{FALSE} a dense kernel matrix
#' of class \code{\link{KernelMatrix}} is returned. If set to \code{TRUE}
#' the kernel matrix is returned as sparse matrix of class
#' \code{\linkS4class{dgCMatrix}}. In case of a symmetric matrix either the
#' lower triangular part or the full matrix can be returned. Please note that
#' a sparse kernel matrix currently can not be used for SVM based learning
#' in \code{kebabs}. Default=FALSE
#' @param triangular boolean indicating whether just the lower triangular or
#' the full sparse matrix should be returned. This parameter is only relevant
#' for a sparse symmetric kernel matrix. Default=TRUE
#' @param diag boolean indicating whether the diagonal should be included
#' in a sparse triangular matrix. This parameter is only relevant when
#' parameter \code{sparse} and \code{triangular} are set to \code{TRUE}.
#' Default=TRUE
#' @param lowerLimit a numeric value for a similarity threshold. The
#' parameter is relevant for sparse kernel matrices only. If set to a value
#' larger than 0 only similarity values larger than this threshold will
#' be included in the sparse kernel matrix. Default=0
#' @return
#' linearKernel:
#' kernel matrix as class \code{\linkS4class{KernelMatrix}} or sparse
#' kernel matrix of class \code{\linkS4class{dgCMatrix}}
#' dependent on parameter \code{sparse}
#' @examples
#'
#' ## load sequence data and change sample names
#' data(TFBS)
#' names(enhancerFB) <- paste("S", 1:length(enhancerFB), sep="_")
#'
#' ## create the kernel object for dimers with normalization
#' speck <- spectrumKernel(k=5)
#'
#' ## generate sparse explicit representation
#' ers <- getExRep(enhancerFB, speck)
#'
#' ## compute dense kernel matrix (as currently used in SVM based learning)
#' km <- linearKernel(ers)
#' km[1:5, 1:5]
#'
#' ## compute sparse kernel matrix
#' ## because it is symmetric just the lower diagonal
#' ## is computed to save storage
#' km <- linearKernel(ers, sparse=TRUE)
#' km[1:5, 1:5]
#'
#' ## compute full sparse kernel matrix
#' km <- linearKernel(ers, sparse=TRUE, triangular=FALSE)
#' km[1:5, 1:5]
#'
#' ## compute triangular sparse kernel matrix without diagonal
#' km <- linearKernel(ers, sparse=TRUE, triangular=TRUE, diag=FALSE)
#' km[1:5, 1:5]
#'
#' ## plot histogram of similarity values
#' hist(as(km, "numeric"), breaks=30)
#'
#' ## compute sparse kernel matrix with similarities above 0.5 only
#' km <- linearKernel(ers, sparse=TRUE, lowerLimit=0.5)
#' km[1:5, 1:5]
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics}, 31(15):2574-2576, 2015.
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.
#' @keywords kernel, linearKernel
#' @export


linearKernel <- function(x, y=NULL, selx=integer(0), sely=integer(0),
                         sparse=FALSE, triangular=TRUE, diag=TRUE,
                         lowerLimit=0)
{
    if (length(selx) > 0)
    {
        if (is.character(selx))
        {
            if (length(rownames(x)) == 0)
                stop("missing rownames in 'x' for subsetting\n")

            selx <- which(rownames(x) %in% selx)
        }
        else
        {
            if (!is.numeric(selx) || any(selx < 1) || any(selx > nrow(x)))
                stop("'selx' is not valid\n")
        }
    }
    else
        selx <- 1L:nrow(x)

    if (!is.null(y))
    {
        if (length(sely) > 0)
        {
            if (is.character(sely))
            {
                if (length(rownames(y)) == 0)
                    stop("missing rownames in 'y' for subsetting\n")

                sely <- which(rownames(y) %in% sely)
            }
            else
            {
                if (!is.numeric(sely) || any(sely < 1) || any(sely > nrow(y)))
                    stop("'sely' is not valid\n")
            }
        }
        else
            sely <- 1L:nrow(y)
    }
    else
    {
        if (length(sely) > 0)
            stop("'sely' only allowed together with 'y'")

        if (sparse == TRUE && triangular == FALSE)
        {
            y <- x

            if (length(selx) > 0)
                sely <- selx
        }
    }

    if (!is(x, "ExplicitRepresentation"))
        stop("'x' must be a dense or sparse explicit representation\n")

    if (!is.null(y) && !is(y, "ExplicitRepresentation"))
        stop("'y' must be a dense or sparse explicit representation\n")

    if (is(x,"ExplicitRepresentationDense"))
    {
        if (!is.null(y) && !is(y, "ExplicitRepresentationDense"))
            stop("'y' must be a dense explicit representation\n")

        if (sparse)
        {
            stop("for sparse kernel matrix the function linearKernel must be\n",
                 "        called with a sparse explicit representation\n")
        }

        if (!is.null(y) && (length(sely) > 0))
            return(as.KernelMatrix(tcrossprod(x=x[selx], y=y[sely])))
        else
            return(as.KernelMatrix(tcrossprod(x=x[selx], y=y)))
    }
    else
    {
        if (!is.null(y) && !is(y, "ExplicitRepresentationSparse"))
            stop("'y' must be a sparse explicit representation\n")

        if (sparse)
        {
            ## set exit hook for C heap cleanup
            on.exit(.Call("freeHeapLinearKernelC"))

            if (length(rownames(x) > 0))
                rowNames <- rownames(x)[selx]
            else
                rowNames <- character(0)

            if (is.null(y))
            {
                colNames <- rowNames

                return(.Call("linearKernelSparseKMdgRMatrixC",
                             as.integer(nrow(x)), x@p, x@j, x@x, selx-1,
                             as.integer(nrow(x)), integer(0), integer(0),
                             double(0), integer(0), rowNames, colNames,
                             TRUE, as.logical(diag), as.double(lowerLimit)))
            }
            else
            {
                if (length(rownames(y) > 0))
                    colNames <- rownames(y)[sely]
                else
                    colNames <- character(0)

                return(.Call("linearKernelSparseKMdgRMatrixC",
                             as.integer(nrow(x)), x@p, x@j, x@x, selx-1,
                             as.integer(nrow(y)), y@p, y@j, y@x, sely-1,
                             rowNames, colNames, FALSE, as.logical(diag),
                             as.double(lowerLimit)))
            }
        }

        if (is.null(y))
        {
            km <- .Call("linearKerneldgRMatrixC", as.integer(nrow(x)), x@p,
                        x@j, x@x, selx-1, as.integer(nrow(x)), integer(0),
                        integer(0), double(0), integer(0), TRUE)

            if (length(rownames(x)) > 0)
            {
                rownames(km) <- rownames(x)[selx]
                colnames(km) <- rownames(x)[selx]
            }
        }
        else
        {
            km <- .Call("linearKerneldgRMatrixC", as.integer(nrow(x)), x@p,
                        x@j, x@x, selx-1, as.integer(nrow(y)), y@p, y@j, y@x,
                        sely-1, FALSE)

            if (length(rownames(x)) > 0)
                rownames(km) <- rownames(x)[selx]

            if (length(rownames(y)) > 0)
                colnames(km) <- rownames(y)[sely]
        }

        return(as.KernelMatrix(km))
    }
}

## subsetting for BioVector, XStringSet, KernelMatrix and ExplicitRep
subsetSeqRep <- function(x, sel)
{
    if (is.null(sel) || length(sel) == 0)
        stop("'sel' is null in subset operation\n")

    if (is(x, "KernelMatrix") || is(x, "kernelMatrix"))
        noElem <- nrow(x)
    else
        noElem <- length(x)

    if (!is.numeric(sel) || any(sel < 1) || any(sel > noElem))
        stop("'sel' is not valid\n")

    if (is(x, "KernelMatrix") || is(x, "kernelMatrix"))
        return(x[sel, sel])
    else
        return(x[sel])
}

## no of elements for BioVector, XStringSet, KernelMatrix and ExplicitRep
getNoOfElementsOfSeqRep <- function(x)
{
    if (is.null(x))
        return(0)

    if (is(x, "BioVector") || is(x, "XStringSet"))
        return(length(x))
    else
        return(nrow(x))
}

listToString <- function(x)
{
    if (!is.list(x))
        stop("'x' must be a list\n")

    temp <- paste(deparse(x), collapse="")
    endPos <- regexpr(".Names = ", temp) - 4

    if (endPos > 0)
        return(substr(temp, 16, endPos))
    else
        return(NULL)
}

#' @rdname ExplicitRepresentationAccessors
#' @aliases
#' %*%,matrix,dgRMatrix-method
#' @param y in the first case and explicit representation and \code{x} is a
#' \code{matrix}, for the second case a numeric matrix and x is an explicit
#' representation
#' @section Accessor-like methods:
#' \describe{
#'   \item{}{\code{\%*\%}:
#'   this function provides the multiplication of a \code{dgRMatrix} or
#'   a sparse explicit representation (which is derived from \code{dgRMatrix})
#'   with a matrix or a vector. This functionality is not available in package
#'   \bold{Matrix} for a \code{dgRMatrix}.
#'   }
#' }
#' @export
#'

## matrix multiplication of dense matrix with dgRMatrix
## is not implemented in Matrix package
setMethod("%*%", signature(x="matrix", y="dgRMatrix"),
    function(x, y)
    {
        if (ncol(x) == 0 && nrow(y) == 0)
            return(matrix(nrow=0,ncol=0))

        if (ncol(x) != nrow(y))
            stop("non-conformable arguments\n")

        ## dense matrix product is much faster if size allows
        if (nrow(y) * ncol(y) < 10000000)
            return(x %*% as(y, "matrix"))

        result <- .Call("matrixdgRMatrixProductC", x, as.integer(nrow(x)),
                        as.integer(ncol(x)), as.integer(nrow(y)),
                        as.integer(ncol(y)), y@p, y@j, y@x)

        if (length(rownames(x)) > 0)
            rownames(result) <- rownames(x)

        if (length(colnames(y)) > 0)
            colnames(result) <- colnames(y)

        return(result)
    }
)

#' @rdname ExplicitRepresentationAccessors
#' @aliases
#' %*%,dgRMatrix,numeric-method
#' @export
#'

## multiplication of dgRMatrix with vector
setMethod("%*%", signature(x="dgRMatrix", y="numeric"),
    function(x, y)
    {
        if (ncol(x) == 0 && length(y) == 0)
            return(matrix(nrow=0,ncol=0))

        if (ncol(x) != length(y))
            stop("non-conformable arguments\n")

        ## dense matrix product is much faster if size allows
        if (nrow(x) * ncol(x) < 10000000)
            return(as(x, "matrix") %*% y)

        result <- .Call("dgRMatrixNumericVectorProductC", x@p, x@j, x@x,
                        as.integer(nrow(x)), as.integer(ncol(y)), y,
                        as.integer(length(y)))

        if (length(rownames(x)) > 0)
            names(result) <- rownames(x)

        return(result)
    }
)

#' @rdname genRandBioSeqs
#' @title Generate Random Biological Sequences
#'
#' @description Generate biological sequences with uniform random distribution
#' of alphabet characters.
#'
#' @param seqType defines the type of sequence as DNA, RNA or AA and the
#' underlying alphabet. Default="DNA"
#'
#' @param numSequences single numeric value which specifies the number of
#' sequences that should be generated.
#'
#' @param seqLength either a single numeric value or a numeric vector of length
#' 'numSequences' which gives the length of the sequences to be generated.
#'
#' @param biostring if \code{TRUE} the sequences will be generated in XStringSet
#' format otherwise as BioVector derived class. Default=TRUE
#'
#' @param seed when present the random generator will be seeded with the
#' value passed in this parameter
#'
#' @details
#' The function generates a set of sequences with uniform distribution of
#' alphabet characters and returns it as XStringSet or BioVector dependent on
#' the parameter \code{biostring}.
#' @return
#' When the parameter 'biostring' is set to FALSE the function returns a
#' XStringSet derived class otherwise a BioVector derived class.
#'
#' @examples
#' ## generate a set of AA sequences of fixed length as AAStringSet
#' aaseqs <- genRandBioSeqs("AA", 100, 1000, biostring=TRUE)
#'
#' ## show AA sequence set
#' aaseqs
#'
#' \dontrun{
#' ## generate a set of "DNA" sequences as DNAStringSet with uniformly
#' ## distributed lengths between 1500 and 3000 bases
#' seqLength <- runif(300, min=1500, max=3500)
#' dnaseqs <- genRandBioSeqs("DNA", 100, seqLength, biostring=TRUE)
#'
#' ## show DNA sequence set
#' dnaseqs
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics}, 31(15):2574-2576, 2015.
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.
#' @keywords methods
#' @export

genRandBioSeqs <- function(seqType=c("DNA", "RNA", "AA"), numSequences,
                           seqLength, biostring=TRUE, seed)
{
    if (!missing(seed))
        set.seed(seed)

    seqType <- match.arg(seqType)

    if (!(is.numeric(numSequences) && length(numSequences) == 1 &&
          numSequences >= 1))
        stop("'numSequences' must be a integer value larger than 0\n")

    if (!(is.numeric(seqLength) && all(seqLength >= 1)))
        stop("'seqLength' must be a numeric vector with values larger than 0\n")

    numSequences <- as.integer(numSequences)
    seqLength <- as.integer(seqLength)

    if (length(seqLength == 1))
        seqLength <- rep(seqLength, numSequences)

    if (seqType == "DNA")
    {
        seqAlphabet <- c("A", "C", "G", "T")
        seqs <- list()

        for (i in 1:numSequences)
        {
            seqs <- list(seqs, paste(seqAlphabet[round(runif(seqLength[i],
                                                    0.5, 4.5))], collapse=""))
        }

        if (biostring == "TRUE")
            return(DNAStringSet(as.character(unlist(seqs))))
        else
            return(new("DNAVector",(as.character(unlist(seqs)))))
    }
    else if (seqType == "RNA")
    {
        seqAlphabet <- c("A", "C", "G", "U")
        seqs <- list()

        for (i in 1:numSequences)
        {
            seqs <- list(seqs, paste(seqAlphabet[round(runif(seqLength[i],
                                                    0.5, 4.5))], collapse=""))
        }

        if (biostring == "TRUE")
            return(RNAStringSet(as.character(unlist(seqs))))
        else
            return(new("RNAVector",(as.character(unlist(seqs)))))
    }
    else if (seqType == "AA")
    {
        seqAlphabet <- c("A","C","D","E","F","G","H","I","K","L","M","N",
                         "P", "Q","R","S","T","U","V","W","Y")
        seqs <- list()

        for (i in 1:numSequences)
        {
            seqs <- list(seqs, paste(seqAlphabet[round(runif(seqLength[i],
                                                    0.5, 21.5))], collapse=""))
        }

        if (biostring == "TRUE")
            return(AAStringSet(as.character(unlist(seqs))))
        else
            return(new("AAVector",(as.character(unlist(seqs)))))
    }
}

#' @rdname computeROCandAUC
#' @title Compute Receiver Operating Characteristic And Area Under The Curve
#'
#' @description Compute the receiver operating characteristic (ROC) and
#' area under the ROC curve (AUC) as performance measure for binary
#' classification
#'
#' @param prediction prediction results in the form of decision values as
#' returned by \code{\link{predict}} for predictionType="decision".
#'
#' @param labels label vector of same length as parameter 'prediction'.
#'
#' @param allLabels vector containing all occuring labels once. This parameter
#' is required only if the labels parameter is not a factor. Default=NULL
#'
#' @details
#' For binary classfication this function computes the receiver operating
#' curve (ROC) and the area under the ROC curve (AUC).
#'
#' @return
#' On successful completion the function returns an object of class
#' \code{\linkS4class{ROCData}} containing the AUC, a numeric vector of
#' TPR values and a numeric vector containing the FPR values. If the ROC and
#' AUC cannot be computed because of missing positive or negative samples the
#' function returns 3 NA values.
#'
#' @seealso \code{\link{predict}}, \code{\linkS4class{ROCData}}
#' @examples
#'
#' ## load transcription factor binding site data
#' data(TFBS)
#' enhancerFB
#' ## select 70% of the samples for training and the rest for test
#' train <- sample(1:length(enhancerFB), length(enhancerFB) * 0.7)
#' test <- c(1:length(enhancerFB))[-train]
#' ## create the kernel object for gappy pair kernel with normalization
#' gappy <- gappyPairKernel(k=1, m=3)
#' ## show details of kernel object
#' gappy
#'
#' ## run training with explicit representation
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=gappy,
#'                pkg="LiblineaR", svm="C-svc", cost=80, explicit="yes",
#'                featureWeights="no")
#'
#' ## predict the test sequences
#' pred <- predict(model, enhancerFB[test])
#' ## print prediction performance
#' evaluatePrediction(pred, yFB[test], allLabels=unique(yFB))
#'
#' ## compute ROC and AUC
#' preddec <- predict(model, enhancerFB[test], predictionType="decision")
#' rocdata <- computeROCandAUC(preddec, yFB[test], allLabels=unique(yFB))
#'
#' ## show AUC value
#' rocdata
#'
#' \dontrun{
#' ## plot ROC
#' plot(rocdata)
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics}, 31(15):2574-2576, 2015.
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.
#' @keywords prediction performance
#' @keywords methods
#' @export

computeROCandAUC <- function(prediction, labels, allLabels=NULL)
{
    if (length(prediction) != length(labels))
        stop("length of prediction does not match length of label\n")

    if (is.factor(labels))
        allLabels <- levels(labels)
    else
    {
        if (!is.null(allLabels))
            allLabels <- unique(allLabels)
        else
            stop("please specify all labels via parameter 'allLabels'\n")
    }

    if (length(allLabels) > 2)
        stop("AUC is only supported for binary classification\n")

    allLabels <- sortWith_LC_Collate_C(allLabels, decreasing=TRUE)
    numSamples <- length(prediction)
    posSamples <- sum(labels == allLabels[1])
    negSamples <- sum(labels == allLabels[2])

    if ((posSamples + negSamples) != numSamples)
        stop("AUC is only supported for binary classification\n")

    if (posSamples == 0 || negSamples == 0)
        return(new("ROCData", AUC=NA, TPR=numeric(0), FPR=numeric(0)))

    ## contribution per sample
    TPC <- 1 / posSamples
    FPC <- 1 / negSamples
    ind <- order(prediction, decreasing=TRUE)
    TPR <- rep(0, length(prediction) + 1)
    FPR <- rep(0, length(prediction) + 1)
    predOld <- -Inf
    FPOld<- 0
    TPOld <- 0
    TP <- 0
    FP <- 0
    AUC <- 0

    for (i in 1:numSamples)
    {
        if (prediction[i] != predOld)
        {
            AUC <- AUC + (FP - FPOld) * (TP + TPOld) / 2
            predOld <- prediction[i]
            FPOld <- FP
            TPOld <- TP
        }

        if (labels[ind[i]] == allLabels[1])
        {
            TP <- TP + 1
            FPR[i+1] <- FPR[i]
            TPR[i+1] <- TPR[i] + TPC
        }
        else
        {
            FP <- FP + 1
            TPR[i+1] <- TPR[i]
            FPR[i+1] <- FPR[i] + FPC
        }
    }

    AUC <- (AUC + (negSamples - FPOld) * (posSamples + TPOld) / 2) /
           (posSamples * negSamples)

    return(new("ROCData", AUC=AUC, TPR=TPR, FPR=FPR))
}

#' @rdname evaluatePrediction
#' @title Evaluate Prediction
#'
#' @description Evaluate performance results of prediction on a testset based
#' on given labels for binary classification
#'
#' @param prediction prediction results as returned by \code{\link{predict}}
#' for predictionType="response".
#'
#' @param label label vector of same length as parameter 'prediction'.
#'
#' @param allLabels vector containing all occuring labels once. This parameter
#' is required only if the label vector is numeric. Default=NULL
#'
#' @param decValues numeric vector containing decision values for the
#' predictions as returned by the \code{\link{predict}} method with
#' \code{predictionType} set to \code{decision}. This parameter is needed for
#' the determination of the AUC value which is currently only supported for
#' binary classification. Default=NULL
#'
#' @param print This parameter indicates whether performance values should be
#' printed or returned as data frame without printing (for details see below).
#' Default=TRUE
#'
#' @param confmatrix When set to TRUE a confusion matrix is printed. The rows
#' correspond to predictions, the columns to the true labels. Default=TRUE
#'
#' @param numPrecision minimum number of digits to the right of the decimal
#' point. Values between 0 and 20 are allowed. Default=3
#'
#' @param numPosNegTrainSamples optional integer vector with two values giving
#' the number of positive and negative training samples. When this parameter
#' is set the balancedness of the training set is reported. Default=numeric(0)
#'
#' @details
#' For binary classfication this function computes the performance measures
#' accuracy, balanced accuracy, sensitivity, specificity, precision and the
#' Matthews Correlation Coefficient(MCC). If decision values are passed in the
#' parameter \code{decValues} the function additionally determines the AUC.
#' When the number of positive and negative training samples is passed to
#' the function it also shows the balancedness of the training set. The
#' performance results are either printed by the routine directly or returned
#' in a data frame. The columns of the data frame are:
#'
#' \tabular{ll}{
#'   column name \tab performance measure\cr
#'   -------------------- \tab -------------- \cr
#'   TP          \tab true positive\cr
#'   FP          \tab false positive\cr
#'   FN          \tab false negative\cr
#'   TN          \tab true negative\cr
#'   ACC         \tab accuracy \cr
#'   BAL_ACC     \tab balanced accuracy\cr
#'   SENS        \tab sensitivity\cr
#'   SPEC        \tab specificity\cr
#'   PREC        \tab precision\cr
#'   MAT_CC      \tab Matthews correlation coefficient\cr
#'   AUC         \tab area under ROC curve\cr
#'   PBAL        \tab prediction balancedness (fraction of positive samples)\cr
#'   TBAL        \tab training balancedness (fraction of positive samples)\cr
#' }
#'
#'
#' @return
#' When the parameter 'print' is set to FALSE the function returns a data frame
#' containing the prediction performance values (for details see above).
#' @seealso \code{\link{predict}}, \code{\link{kbsvm}}
#' @examples
#'
#' ## set seed for random generator, included here only to make results
#' ## reproducable for this example
#' set.seed(456)
#' ## load transcription factor binding site data
#' data(TFBS)
#' enhancerFB
#' ## select 70% of the samples for training and the rest for test
#' train <- sample(1:length(enhancerFB), length(enhancerFB) * 0.7)
#' test <- c(1:length(enhancerFB))[-train]
#' ## create the kernel object for gappy pair kernel with normalization
#' gappy <- gappyPairKernel(k=1, m=3)
#' ## show details of kernel object
#' gappy
#'
#' ## run training with explicit representation
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=gappy,
#'                pkg="LiblineaR", svm="C-svc", cost=80, explicit="yes",
#'                featureWeights="no")
#'
#' ## predict the test sequences
#' pred <- predict(model, enhancerFB[test])
#'
#' ## print prediction performance
#' evaluatePrediction(pred, yFB[test], allLabels=unique(yFB))
#'
#' \dontrun{
#' ## print prediction performance including AUC
#' ## additionally determine decision values
#' preddec <- predict(model, enhancerFB[test], predictionType="decision")
#' evaluatePrediction(pred, yFB[test], allLabels=unique(yFB),
#'                    decValues=preddec)
#'
#' ## print prediction performance including training set balance
#' trainPosNeg <- c(length(which(yFB[train] == 1)),
#'                  length(which(yFB[train] == -1)))
#' evaluatePrediction(pred, yFB[test], allLabels=unique(yFB),
#'                    numPosNegTrainSamples=trainPosNeg)
#'
#' ## or get prediction performance as data frame
#' perf <- evaluatePrediction(pred, yFB[test], allLabels=unique(yFB),
#'                            print=FALSE)
#'
#' ## show performance values in data frame
#' perf
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics}, 31(15):2574-2576, 2015.
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.
#' @keywords prediction performance
#' @keywords methods
#' @export

evaluatePrediction <- function(prediction, label, allLabels=NULL,
                            decValues=NULL, print=TRUE, confmatrix=TRUE,
                            numPrecision=3, numPosNegTrainSamples=numeric(0))
{
    ## for binary classification only - first label must be the positive class
    ## optional input pos and neg training samples for training balance

    if (length(prediction) != length(label))
        stop("length of 'prediction' and 'label' do not match\n")
    
    computeAUC <- FALSE

    if (is.factor(label))
        allLabels <- levels(label)
    else
    {
        if (!is.null(allLabels))
            allLabels <- unique(allLabels)
        else
            stop("please specify all labels via parameter 'allLabels'\n")
    }

    if (length(allLabels) > 2)
        allLabels <- sortWith_LC_Collate_C(allLabels)
    else
    {
        allLabels <- sortWith_LC_Collate_C(allLabels, decreasing=TRUE)
        
        ## AUC currently only for binary classification
        if (!is.null(decValues) && !(is.na(decValues)))
        {
            if (length(decValues) != length(label))
                stop("length of 'decValues' and 'label' do not match\n")
            
            computeAUC <- TRUE
        }
    }
    x <- matrix(0, nrow=length(allLabels), ncol=length(allLabels))
    rownames(x) <- allLabels
    colnames(x) <- allLabels

    for (i in 1:length(prediction))
    {
        predChar <- as.character(prediction[i])
        labelChar <- as.character(label[i])
        x[predChar, labelChar] <- x[predChar, labelChar] + 1
    }

    if(confmatrix && print)
    {
        print(x)
        cat("\n")
    }

    auc <- NA

    if (length(allLabels) == 2)
    {
        tp <- x[1, 1]
        tn <- x[2, 2]
        fp <- x[1, 2]
        fn <- x[2, 1]
        l <- sum(x)
        bal = sum(x[,1]) / l
        accuracy = 100 * (tp + tn) / l

        if (tp + fn > 0 && tn + fp > 0)
        {
            balAccuracy <- 50 * (tp / (tp + fn) + tn / (tn + fp))
            sensitivity <- 100 * tp / (tp + fn)
            specificity <- 100 * tn / (tn + fp)
        }
        else if (tp + fn > 0)
        {
            balAccuracy <- NA
            sensitivity <- 100 * tp / (tp + fn)
            specificity <- NA
        }
        else
        {
            balAccuracy <- NA
            sensitivity <- NA
            specificity <- 100 * tn / (tn + fp)
        }

        precision <- ifelse(tp + fp > 0, 100 * tp / (tp + fp), NA)
        denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fn) * (tn + fp))
        mcc <- ifelse(denominator == 0, 0, (tp * tn - fp * fn) / denominator)

        if (length(numPosNegTrainSamples) > 0)
        {
            trainBalance <- 100 * numPosNegTrainSamples[1] /
            (numPosNegTrainSamples[1] + numPosNegTrainSamples[2])
        }
        else
            trainBalance <- NA
            
        if ((computeAUC == TRUE) && sum(x[,1]) > 0 && sum(x[,2]) > 0)
            auc <- auc(computeROCandAUC(decValues, label, allLabels))
    }
    else
    {
        ## Multiclass
        bal <- NA
        trainBalance <- NA
        l <- sum(x)
        corrClassified <- diag(x)
        accuracy <- 100 * sum(corrClassified) / l
        balAccuracy <- 0

        for (i in 1:length(allLabels))
            balAccuracy <- balAccuracy + corrClassified[i] / sum(x[,i])

        balAccuracy <- 100 * balAccuracy / length(allLabels)

        ## Multiclass MCC according to Gorodkin J.,
        ## Comparing two K-category assignments by K-category correlation
        ## coefficient, Comput. Biol. Chem., 2004 Dec, 28(5-6) 367-74.
        ## according to fomula 8 in paper
        ##
        mcc <- (l * sum(diag(x)) - sum(x %*% x)) / sqrt(
                   (l^2 - sum(tcrossprod(x))) * (l^2 - sum(crossprod(x))))

        sensitivity <- NA
        specificity <- NA
        precision <- NA
        tp <- NA
        tn <- NA
        fp <- NA
        fn <- NA
    }

    if (print)
    {
        if (length(allLabels) > 2)
        {
            cat("Accuracy:            ",
                format(accuracy, width=8, digits=3, nsmall=numPrecision,
                       trim=FALSE), "% (", sum(corrClassified), " of ", l, 
                       ")\n", sep="")
            cat("Balanced accuracy:   ",
                format(balAccuracy, width=8, digits=3, nsmall=numPrecision,
                       trim=FALSE), "%\n", sep="")
            cat("Matthews CC:         ",
                format(mcc, width=8, digits=3, nsmall=numPrecision, trim=FALSE),
                "\n", sep="")

            for (i in 1:length(allLabels))
            {
                cat("\nClass ", allLabels[i],":\n", sep="")
                cat("Sensitivity:         ",
                    format(100 * corrClassified[i] / sum(x[,i]), width=8,
                           digits=3, nsmall=numPrecision, trim=FALSE),
                    "% (", corrClassified[i] , " of ", sum(x[,i]), ")\n", 
                    sep="")
                cat("Specificity:         ",
                    format(100 * sum(corrClassified[-i]) / sum(x[,-i]), width=8,
                           digits=3, nsmall=numPrecision, trim=FALSE),
                    "% (", sum(corrClassified[-i]), " of ", sum(x[,-i]), ")\n",
                    sep="")
                cat("Precision:           ",
                    format(100 * corrClassified[i] / sum(x[i,]), width=8,
                           digits=3, nsmall=numPrecision, trim=FALSE),
                    "% (", corrClassified[i], " of ", sum(x[i,]), ")\n", sep="")
            }
        }
        else
        {
            if (length(numPosNegTrainSamples) > 0)
            {
                cat("Training set balance:",
                    format(trainBalance, width=8, digits=3, 
                           nsmall=numPrecision),
                    "% (", numPosNegTrainSamples[1], " positive of ",
                    numPosNegTrainSamples[1] + numPosNegTrainSamples[1],
                    ")\n", sep="")
                cat("\n")
            }

            cat("Accuracy:            ",
                format(accuracy, width=8, digits=3, nsmall=numPrecision,
                       trim=FALSE), "% (", tp + tn, " of ", l, ")\n", sep="")
            cat("Balanced accuracy:   ",
                format(balAccuracy, width=8, digits=3, nsmall=numPrecision,
                       trim=FALSE), "% (", tp, " of ", tp + fn, " and ",
                tn, " of ", tn + fp, ")\n", sep="")
            cat("Matthews CC:         ",
                format(mcc, width=8, digits=3, nsmall=numPrecision, trim=FALSE),
                "\n", sep="")
                
            if (!is.na(auc))
            {
                cat("AUC:                 ",
                    format(auc, width=8, digits=3, nsmall=numPrecision,
                    trim=FALSE),
                    "\n", sep="")
            }
            cat("\n")

            cat("Sensitivity:         ",
                format(sensitivity, width=8, digits=3, nsmall=numPrecision,
                       trim=FALSE), "% (", tp, " of ", tp + fn, ")\n", sep="")
            cat("Specificity:         ",
                format(specificity, width=8, digits=3, nsmall=numPrecision,
                       trim=FALSE), "% (", tn, " of ", tn + fp, ")\n", sep="")
            cat("Precision:           ",
                format(precision, width=8, digits=3, nsmall=numPrecision,
                       trim=FALSE), "% (", tp, " of ", tp + fp, ")\n", sep="")
        }
        return(invisible(NULL))
    }
    else
    {
        if (length(numPosNegTrainSamples) > 0)
        {
            return(data.frame(NUM=l, PBAL=bal, TP=tp, FP=fp, FN=fn, TN=tn,
                              ACC=accuracy, BAL_ACC=balAccuracy,
                              SENS=sensitivity, SPEC=specificity,
                              PREC=precision, MAT_CC=mcc, AUC=auc,
                              TBAL=trainBalance,
                              PT=numPosNegTrainSamples[1],
                              NT=numPosNegTrainSamples[2]))
        }
        else
        {
            return(data.frame(NUM=l, PBAL=bal, TP=tp, FP=fp, FN=fn, TN=tn,
                              ACC=accuracy, BAL_ACC=balAccuracy,
                              SENS=sensitivity, SPEC=specificity,
                              PREC=precision, MAT_CC=mcc, AUC=auc))
        }
    }
}
