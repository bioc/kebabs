#2345678901234567890123456789012345678901234567890123456789012345678901234567890
#' @rdname symmetricPairKernel
#' @title Symmetric Pair Kernel
#'
#' @description Create a symmetric pair kernel object
#'
#' @param siKernel kernel for single instances
#'
#' @param kernelType defines the type of pair kernel. It specifies in which
#' way the similarity between two pairs of sequences are computed. Allowed
#' values are "mean", and "TPPK" (see also details section). Default="mean"
#'
#' @param r exponent which must be > 0 (details see below). Default=1
#'
#' @details
#' Creation of kernel object\cr\cr
#' The function 'symmetricPairKernel' creates a kernel object for the symmetric
#' pair kernel. This kernel is an example for multiple instance learning and
#' can be used for learning based on pairs of sequences. The single instance
#' kernel passed to the symmetric pair kernel computes a similarity between
#' two individual sequences giving a similarity for one pair of sequences.
#' The symmetric pair kernel function gets as input two pairs of sequences and
#' computes a similarity value between the two pairs. This similarity is
#' computed dependent on the value of the argument \code{kernelType} from the
#' similarities delivered by the single instance kernel in the following
#' way:\cr\cr
#'
#' \code{mean} (arithmetic mean):\cr\cr
#' k(<a,b>, <c,d>) = 1/4 * (k(a,c) + k(a,d) + k(b,c) + k(b,d))
#'
#' \code{TPKK} (tensor pairwise product kernel):\cr\cr
#' k(<a,b>, <c,d>) = (k(a,c) * k(b,d) + k(a,d) * k(b,c))
#'
#' Every sequence kernel available in KeBABS can be used as single instance
#' kernel for the symmetric pair kernel allowing to create similarity
#' measures between two pairs of sequences based on different similarity
#' measures between individual sequences.
#'
#' The row names and column names of a kernel matrix generated from a symmetric
#' pair kernel object describe the sequence pair with the names of the
#' individual sequences in the pair separated by the underscore character.
#'
#' For values different from 1 (=default value) parameter \code{r}
#' leads to a transfomation of similarities by taking each element of the
#' similarity matrix to the power of r. Only integer values larger than 1
#' should be used for r in context with SVMs requiring positive definite
#' kernels.
#'
#' The symmetricPairKernel can be used in sequence based learning like any
#' single instance kernel. Label values are defined against pairs of sequences
#' in this case. Explicit representation, feature weights and prediction
#' profiles are not available for the symmetric pair kernel. As kernels
#' computed through sums and products of postive definite kernels all variants
#' of this kernel are positive definite.
#'
#'
#' @return
#' symmetricPairKernel: upon successful completion, the function returns a
#' kernel object of class \code{\linkS4class{SymmetricPairKernel}}.
#'
#' @seealso \code{\link{kernelParameters-method}},
#' \code{\link{getKernelMatrix}}, \code{\link{spectrumKernel}},
#' \code{\link{mismatchKernel}}, \code{\link{motifKernel}},
#' \code{\link{gappyPairKernel}}, \code{\linkS4class{SymmetricPairKernel}}
#' @examples
#' ## load sample sequences from transcription factor binding dataset
#' data(TFBS)
#' ## in this example we just use the first 30 sequences and rename samples
#' x <- enhancerFB[1:30]
#' names(x) <- paste("S", 1:length(x), sep="")
#'
#' ## create the single instance kernel object
#' specK5 <- spectrumKernel(k=5)
#' ## show details of single instance kernel object
#' specK5
#'
#' ## create the symmetric pair kernel object for the single instance kernel
#' tppk <- symmetricPairKernel(siKernel=specK5, kernelType="TPPK")
#'
#' ## generate the kernel matrix with the symmetric pair kernel object which
#' ## contains similarity values between two pairs of sequences.
#' ## Hint: The kernel matrix for the single instance kernel is computed
#' ## internally.
#' km <- tppk(x)
#' dim(km)
#' km[1:5,1:5]
#'
#' \dontrun{
#' ## plot heatmap of the kernel matrix
#' heatmap(km, symm=TRUE)
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs/}\cr\cr
#' (Hue, 2002) -- M.Hue and J.-P.Vert. On learning with kernels for unordered
#' pairs. \cr\cr
#' (Ben-Hur, 2005) -- A. Ben-Hur and W.S. Noble. Kernel methods for predicting
#' protein-protein interactions. \cr\cr
#' (Gaertner, 2002) -- T. Gaertner, P.A. Flach, A. Kowalczyk, A.J. Smola.
#' Multi-Instance Kernels. \cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics}, 31(15):2574-2576, 2015.
#' DOI: \doi{10.1093/bioinformatics/btv176}.
#' @keywords kernel
#' @keywords symmetric, pair, symmetricPairKernel, multiple instance learning
#' @keywords methods
#' @export
#'

symmetricPairKernel <- function(siKernel, kernelType=c("mean", "TPPK"),
                                r=1)
{
    ## define function for kernel processing

    rval<- function(x, y = NULL, selx=NULL, sely = NULL)
    {
        return(symmPairProcessing(x=x, y=y, selx=selx, sely=sely,
                                  siKernel=siKernel, kernelType=kernelType,
                                  r=r))
    }

    if ((length(siKernel) != 1) ||
        (!is(siKernel, "SequenceKernel") && !is(siKernel, "stringkernel")))
        stop("'siKernel' must be a single SequenceKernel or stringkernel\n")

    if (isUserDefined(siKernel))
        stop("User defined single instance kernel not supported\n")

    if (!isSingleNumber(r) || r <= 0)
        stop("r must be a number larger than 0\n")

    kernelType <- match.arg(kernelType)

    return(new("SymmetricPairKernel", .Data=rval, .userDefKernel=FALSE,
               siKernel=siKernel, kernelType=kernelType, r=r))
}

symmPairProcessing <- function(x, y, selx, sely, siKernel, kernelType, r,
                               self=NULL)
{
    if (!is.null(self))
    {
        ## retrieval of kernel parameters
        return(list(siKernel=self@siKernel, kernelType=self@kernelType,
                    r=self@r))
    }

    if (missing(x) || is.null(x))
    {
        stop(paste("x must be a",
                   paste(kebabsInfo@allowedSeqClasses, collapse=", "), "\n"))
    }

    if (length(x) < 1)
        stop("sequence info is missing in symmetric pair kernel\n")

    if (missing(y))
        y <- NULL

    if (missing(selx) || is.null(x))
        selx <- integer(0)

    if (missing(sely) || is.null(x))
        sely <- integer(0)

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
                   paste(kebabsInfo@allowedSeqClasses, collapse=", "), "\n"))
    }

    if (!is.null(selx))
    {
        if (!is.matrix(selx) || (ncol(selx) != 2) ||
            (nrow(selx) > length(x)^2))
            stop("'selx' must be a numeric matrix defining pairs in 'x'\n")
    }
    else
    {
        selx <- matrix(c(rep(1:(length(x)-1), (length(x)-1):1),
                         unlist(sapply(2:length(x),
                                       function(i) seq(i,length(x))))),
                       ncol=2)
        colnames(selx) <- c("Inst1","Inst2")
    }

    indX <- unique(as.numeric(selx))

    if (min(indX) < 1 || max(indX) > length(x))
        stop("indices in 'selx' must be between 1 and length of 'x'\n")

    symmetric <- TRUE

    if (!is.null(y))
    {
        symmetric <- FALSE

        if (class(y) %in% c("DNAString", "RNAString", "AAString"))
        {
            y <- switch(class(y),
                        "DNAString" = DNAStringSet(y),
                        "RNAString" = RNAStringSet(y),
                        "AAString"  = AAStringSet(y)
                        )
        }

        if (!(class(y) %in% kebabsInfo@allowedSeqSetClasses))
        {
            stop(paste("'y' must be a",
                       paste(kebabsInfo@allowedSeqClasses, collapse=", ")))
        }

        if (class(x) != class(y))
            stop("'x' and 'y' must be of matching classes\n")

        if (!is.null(sely))
        {
            if (!is.matrix(sely) || (ncol(sely) != 2) ||
                (nrow(sely) > length(y)^2))
            stop("'sely' must be a numeric matrix defining pairs in 'y'\n")
        }
        else
        {
            sely <- matrix(c(rep(1:(length(y)-1), (length(y)-1):1),
                             unlist(sapply(2:length(y),
                                           function(i) seq(i,length(y))))),
                           ncol=2)
            colnames(sely) <- c("Inst1","Inst2")
        }

        indY <- unique(as.numeric(sely))

        if (min(indY) < 1 || max(indY) > length(y))
            stop("indices in 'sely' must be between 1 and length of 'y'\n")
    }

    kernelTypeInt <- which(c("mean", "TPPK") == kernelType)

    if (symmetric)
    {
        kmSI <- siKernel(x, selx=indX)

        selxC <- selx - 1

        km <- .Call("symmetricPairKernelC", kmSI, selxC, selxC,
                    as.integer(nrow(selxC)), as.integer(nrow(selxC)),
                    as.integer(kernelTypeInt), as.logical(symmetric))
    }
    else
    {
        kmSI <- siKernel(x, y, selx=indX, sely=indY)

        selxC <- selx - 1
        selyC <- sely - 1

        km <- .Call("symmetricPairKernelC", kmSI, selxC, selyC,
                    as.integer(nrow(selxC)), as.integer(nrow(selyC)),
                    as.integer(kernelTypeInt), as.logical(symmetric))
    }

    if (length(names(x)) > 0)
    {
        rownames(km) <- sapply(1:nrow(selx),
                               function(i) paste(names(x)[selx[i,1]],
                                                 names(x)[selx[i,2]], sep="_"))

        if (is.null(y))
            colnames(km) <- rownames(km)
    }

    if (length(names(y)) > 0)
    {
        colnames(km) <- sapply(1:nrow(sely),
                               function(i) paste(names(y)[sely[i,1]],
                                                 names(y)[sely[i,2]], sep="_"))
    }

    if (r > 1)
        return(as.KernelMatrix(km^r))
    else
        return(as.KernelMatrix(km))
}

setMethod("kernelParameters", signature=signature(object="SymmetricPairKernel"),
          function(object){object(self=object)})

