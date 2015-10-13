#2345678901234567890123456789012345678901234567890123456789012345678901234567890
#' @rdname mismatchKernel
#' @title Mismatch Kernel
#'
#' @description Create a mismatch kernel object and the kernel matrix
#'
#' @param k length of the substrings also called kmers; this parameter defines
#' the size of the feature space, i.e. the total number of features considered
#' in this kernel is |A|^k, with |A| as the size of the alphabet (4 for DNA
#' and RNA sequences and 21 for amino acid sequences). Default=3
#'
#' @param m number of maximal mismatch per kmer. The allowed value range is
#' between 1 and k-1. The processing effort for this kernel is highly dependent
#' on the value of m and only small values will allow efficient processing.
#' Default=1
#'
#' @param r exponent which must be > 0 (see details section in
#' \link{spectrumKernel}). Default=1
#'
#' @param normalized a kernel matrix or explicit representation generated with
#' this kernel will be normalized(details see below). Default=TRUE
#'
#' @param exact use exact character set for the evaluation (details see below).
#' Default=TRUE
#'
#' @param ignoreLower ignore lower case characters in the sequence. If the
#' parameter is not set lower case characters are treated like uppercase.
#' Default=TRUE
#'
#' @param presence if this parameter is set only the presence of a kmers will
#' be considered, otherwise the number of occurances of the kmer is used.
#' Default=FALSE
#' @details
#' Creation of kernel object\cr\cr
#' The function 'mismatchKernel' creates a kernel object for the mismatch
#' kernel. This kernel object can then be used with a set of DNA-, RNA- or
#' AA-sequences to generate a kernel matrix or an explicit representation for
#' this kernel. For values different from 1 (=default value) parameter
#' \code{r} leads to a transfomation of similarities by taking each element of
#' the similarity matrix to the power of r. If \code{normalized=TRUE}, the
#' feature vectors are scaled to the unit sphere before computing the
#' similarity value for the kernel matrix. For two samples with the feature
#' vectors \code{x} and \code{y} the similarity is computed as:
#' \deqn{s=\frac{\vec{x}^T\vec{y}}{\|\vec{x}\|\|\vec{y}\|}}{s=(x^T y)/(|x| |y|)}
#' For an explicit representation generated with the feature map of a
#' normalized kernel the rows are normalized by dividing them through their
#' Euclidean norm. For parameter \code{exact=TRUE} the sequence characters
#' are interpreted according to an exact character set. If the flag is not
#' set ambigous characters from the IUPAC characterset are also evaluated.
#' The annotation specific variant (for details see \link{positionMetadata})
#' and the position dependent variant (for details see
#' \link{annotationMetadata}) are not available for this kernel.\cr\cr
#' Creation of kernel matrix\cr\cr
#' The kernel matrix is created with the function \code{\link{getKernelMatrix}}
#' or via a direct call with the kernel object as shown in the examples below.
#' @return
#' mismatchKernel: upon successful completion, the function returns a kernel
#' object of class \code{\linkS4class{MismatchKernel}}.
#' @seealso \code{\link{kernelParameters}}, \code{\link{getKernelMatrix}},
#' \code{\link{getExRep}}, \code{\link{spectrumKernel}},
#' \code{\link{gappyPairKernel}}, \code{\link{motifKernel}},
#' \code{\linkS4class{MismatchKernel}}
#' @examples
#'
#' ## instead of user provided sequences in XStringSet format
#' ## for this example a set of DNA sequences is created
#' ## RNA- or AA-sequences can be used as well with the mismatch kernel
#' dnaseqs <- DNAStringSet(c("AGACTTAAGGGACCTGGTCACCACGCTCGGTGAGGGGGACGGGGTGT",
#'                           "ATAAAGGTTGCAGACATCATGTCCTTTTTGTCCCTAATTATTTCAGC",
#'                           "CAGGAATCAGCACAGGCAGGGGCACGGCATCCCAAGACATCTGGGCC",
#'                           "GGACATATACCCACCGTTACGTGTCATACAGGATAGTTCCACTGCCC",
#'                           "ATAAAGGTTGCAGACATCATGTCCTTTTTGTCCCTAATTATTTCAGC"))
#' names(dnaseqs) <- paste("S", 1:length(dnaseqs), sep="")
#'
#' ## create the kernel object with one mismatch per kmer
#' mm <- mismatchKernel(k=2, m=1, normalized=FALSE)
#' ## show details of kernel object
#' mm
#'
#' ## generate the kernel matrix with the kernel object
#' km <- mm(dnaseqs)
#' dim(km)
#' km[1:5, 1:5]
#'
#' ## alternative way to generate the kernel matrix
#' km <- getKernelMatrix(mm, dnaseqs)
#' km[1:5,1:5]
#'
#' \dontrun{
#' ## plot heatmap of the kernel matrix
#' heatmap(km, symm=TRUE)
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' (Leslie, 2002) -- C. Leslie, E. Eskin, J. Weston and W.S. Noble. 
#' Mismatch String Kernels for SVM Protein Classification. \cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics}, 31(15):2574-2576, 2015.
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.
#' @keywords kernel
#' @keywords mismatchKernel, mismatch
#' @keywords methods
#' @export


mismatchKernel <- function(k=3, m=1, r=1, normalized=TRUE, exact=TRUE,
                     ignoreLower=TRUE, presence=FALSE)
{
    ## check data independent kernel parameters and create closure

    if (!is.numeric(k) || any(k < 1))
        stop("k must be an integer larger than 0\n")

    if (!is.numeric(m) || any(m < 1) || any(sapply(k, function(ki) any(m >= ki))))
        stop("m must be an integer larger than 0 and smaller than k\n")

    if (!isSingleNumber(r) || r <= 0)
        stop("r must be a number larger than 0\n")

    if (!isTRUEorFALSE(normalized))
        stop("normalized must be TRUE or FALSE\n")

    if (!isTRUEorFALSE(exact))
        stop("exact must be TRUE or FALSE\n")

    if (!isTRUEorFALSE(ignoreLower))
        stop("ignoreLower must be TRUE or FALSE\n")

    if (!isTRUEorFALSE(presence))
        stop("presence must be TRUE or FALSE\n")

    if (length(k) == 1 && length(m) == 1)
    {
        rval<- function(x, y = NULL, selx = NULL, sely = NULL, self=NULL)
        {
            return(mismatchProcessing(x=x, y=y, selx=selx, sely=sely, k=k, m=m,
                                      r=r, normalized=normalized, exact=exact,
                                      ignoreLower=ignoreLower,
                                      presence=presence, self=self))
        }

        return(new("MismatchKernel", .Data=rval, .userDefKernel=FALSE, k=k,
                   m=m, r=r, normalized=normalized, annSpec=FALSE,
                   distWeight=numeric(0), exact=exact, ignoreLower=ignoreLower,
                   presence=presence))
    }
    else
    {
        kmPairs <- as.matrix(expand.grid(m,k))
        colnames(kmPairs) <- NULL

        ## return list of kernel objects
        kernels <- mapply(mismatchKernel, k=kmPairs[,2], m=kmPairs[,1],
                          MoreArgs=list(r=r, normalized=normalized,
                          exact=exact, ignoreLower=ignoreLower,
                          presence=presence))
        return(kernels)
    }
}

mismatchProcessing <- function(x, y, selx, sely, k, m, r, normalized,
                               exact, ignoreLower, presence, self=NULL)
{
    if (!is.null(self))
    {
        ## retrieval of kernel parameters
        return(list(k=self@k, m = self@m, r=self@r,
                    normalized=self@normalized,
                    annSpec=FALSE,
                    exact=self@exact,
                    ignoreLower=self@ignoreLower,
                    presence=self@presence,
                    distWeight=numeric(0),
                    revComplement=self@revComplement))
    }

    if (missing(x) || is.null(x))
    {
        stop(paste("x must be a",
                   paste(kebabsInfo@allowedSeqClasses, collapse=", "), "\n"))
    }

    if (length(x) < 1)
        stop("sequence info is missing\n")

    if (missing(y))
        y <- NULL

    if (missing(selx) || is.null(selx))
        selx <- integer(0)

    if (missing(sely) || is.null(sely))
        sely <- integer(0)

    if (class(x) %in% c("DNAString", "RNAString", "AAString"))
    {
        x <- switch(class(x),
                    "DNAString" = DNAStringSet(x),
                    "RNAString" = RNAStringSet(x),
                    "AAString"  = AAStringSet(x)
                    )
    }

    if (!(class(x) %in% kebabsInfo@allowedSeqSetClasses))
    {
        stop(paste("x must be a",
                   paste(kebabsInfo@allowedSeqClasses, collapse=", ")))
    }

    if (length(selx) > 0)
    {
        if (!is.numeric(selx) || length(selx) > length(x))
            stop("selx must be a numeric vector with indices into 'x'\n")

        selx <- as.integer(selx)
    }
    else
        selx <- 1L:length(x)

    selxC <- selx - 1L
    symmetric <- TRUE

    if (length(k) > 1 || length(m) > 1)
    {
        stop("multiple values for 'k' or 'm' are only allowed\n",
             "        in model selection\n")
    }

    if (!is.null(y))
    {
        symmetric <- FALSE

        if (class(y) %in% c("DNAString", "RNAString", "AAString"))
        {
            y <- switch(class(y), "DNAString" = DNAStringSet(y),
                        "RNAString" = RNAStringSet(y),
                        "AAString"  = AAStringSet(y)
                        )
        }

        if (!(class(y) %in% kebabsInfo@allowedSeqSetClasses))
        {
            stop(paste("y must be a",
                       paste(kebabsInfo@allowedSeqClasses, collapse=", "),"\n"))
        }

        if (class(x) != class(y))
            stop("x and y must be of matching classes\n")

        if (length(sely) > 0)
        {
            if (!is.numeric(sely) || length(sely) > length(y))
                stop("sely must be a numeric vector with indices into 'y'\n")

            sely <- as.integer(sely)
        }
        else
            sely <- 1L:length(y)

        selyC <- sely - 1L
    }
    else
    {
        sely <- integer(0)
        selyC <- sely
    }

    bioCharset <- getBioCharset(x, exact)

    ## limit k to 32 bit feature space
    if (names(bioCharset[[1]]) %in% c("AAexact", "AAiupac") && k > 7)
        stop("'k' must be smaller than or equal to 7\n")
    else if (names(bioCharset[[1]]) %in% c("DNAexact", "RNAexact") && k > 15)
        stop("for exact charset 'k' must be smaller than or equal to 15\n")
    else if (names(bioCharset[[1]]) %in% c("DNAiupac", "RNAiupac") &&
             k > 7)
        stop("for iupac charset 'k' must be smaller than or equal to 7\n")

    if (!is.null(y))
        seqLength <- c(width(x)[selx], width(y)[sely])
    else
        seqLength <- width(x)[selx]

    if (any(seqLength < k))
    {
        stop("mismatch kernel does not accept strings shorter\n",
             "       than k\n")
    }

    maxSeqLength <- max(seqLength)
    isXStringSet <- is(x, "XStringSet")
    unmapped <- is(x, "DNAStringSet") || is(x, "RNAStringSet")

    res <- .Call("mismatchKernelMatrixC", x, y, selxC, selyC,
                 as.integer(length(selxC)), as.integer(length(selyC)),
                 as.logical(isXStringSet), as.logical(symmetric),
                 as.integer(bioCharset[[2]]), as.logical(!ignoreLower),
                 as.logical(unmapped), as.integer(maxSeqLength),
                 as.integer(k), as.integer(m), as.logical(normalized),
                 as.logical(presence))

    if (length(names(x)) > 0 && length(selx) == nrow(res))
    {
        rownames(res) <- names(x)[selx]

        if (symmetric)
            colnames(res) <- rownames(res)
    }

    if (!symmetric && length(names(y)) > 0 && length(sely) == ncol(res))
    {
        colnames(res) <- names(y)[sely]
    }

    if (r != 1)
        return(as.KernelMatrix(res^r))
    else
        return(as.KernelMatrix(res))
}

#' @rdname mismatchKernel
#' @aliases
#' getFeatureSpaceDimension,MismatchKernel-method
#'
#' @param kernel a sequence kernel object
#' @param x one or multiple biological sequences in the form of a
#' \code{\linkS4class{DNAStringSet}}, \code{\linkS4class{RNAStringSet}},
#' \code{\linkS4class{AAStringSet}} (or as \code{\linkS4class{BioVector}})
#' @return of getDimFeatureSpace:
#' dimension of the feature space as numeric value
#' @export
#'

## feature space is identical to spectrum kernel
setMethod("getFeatureSpaceDimension",
          signature=signature(kernel="MismatchKernel"),
          kebabs:::getFeatureSpaceDimension.spectrum)
