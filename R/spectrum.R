#2345678901234567890123456789012345678901234567890123456789012345678901234567890
#' @rdname spectrumKernel
#' @title Spectrum Kernel
#'
#' @description Create a spectrum kernel object
#'
#' @param k length of the substrings (also called kmers). This parameter defines
#' the size of the feature space, i.e. the total number of features considered
#' in this kernel is |A|^k, with |A| as the size of the alphabet (4 for DNA
#' and RNA sequences and 21 for amino acid sequences). When multiple kernels
#' with different k values should be generated e.g. for model selection a
#' range e.g. k=3:5 can be specified. In this case a list of kernel objects
#' with the individual k values from the range is generated as result.
#' Default=3
#'
#' @param r exponent which must be > 0 (details see below). Default=1
#'
#' @param annSpec boolean that indicates whether sequence annotation should
#' be taken into account (details see on help page for
#' \code{\link{annotationMetadata}}). For the annotation specific spectrum
#' kernel the total number of features increases to |A|^k * |a|^k with
#' |A| as the size of the sequence alphabet and |a| as the size of the
#' annotation alphabet. Default=FALSE
#'
#' @param distWeight a numeric distance weight vector or a distance weighting
#' function (details see on help page for \code{\link{gaussWeight}}).
#' Default=NULL
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
#'
#' @param revComplement if this parameter is set a kmer and its reverse
#' complement are treated as the same feature. Default=FALSE
#'
#' @param mixCoef mixing coefficients for the mixture variant of the spectrum
#' kernel. A numeric vector of length k is expected for this parameter with
#' the unused components in the mixture set to 0. Default=numeric(0)
#'
#' @details
#' Creation of kernel object\cr\cr
#' The function 'spectrumKernel' creates a kernel object for the spectrum
#' kernel. This kernel object can then be used with a set of DNA-, RNA- or
#' AA-sequences to generate a kernel matrix or an explicit representation for
#' this kernel. The spectrum kernel uses all subsequences for length k (also
#' called kmers). For sequences shorter than k the self similarity (i.e. the
#' value on the main diagonal in the square kernel matrix) is 0. The explicit
#' representation contains only zeros for such a sample. Dependent on the
#' learning task it might make sense to remove such sequences from the data
#' set as they do not contribute to the model but still influence performance
#' values. \cr\cr
#' For values different from 1 (=default value) parameter \code{r}
#' leads to a transfomation of similarities by taking each element of the
#' similarity matrix to the power of r. Only integer values larger than 1
#' should be used for r in context with SVMs requiring positive definite
#' kernels. If \code{normalized=TRUE}, the feature vectors are scaled to the
#' unit sphere before computing the similarity value for the kernel matrix.
#' For two samples with the feature vectors \code{x} and \code{y} the
#' similarity is computed as:
#' \deqn{s=\frac{\vec{x}^T\vec{y}}{\|\vec{x}\|\|\vec{y}\|}}{s=(x^T y)/(|x| |y|)}
#' For an explicit representation generated with the feature map of a
#' normalized kernel the rows are normalized by dividing them through their
#' Euclidean norm. For parameter \code{exact=TRUE} the sequence characters
#' are interpreted according to an exact character set. If the flag is not
#' set ambigous characters from the IUPAC characterset are also evaluated.
#' For sequences shorter than k the self similarity (i.e. the value on the
#' main diagonal in the square kernel matrix) is 0.
#'
#' The annotation specific variant (for details see
#' \link{annotationMetadata}) and the position dependent variants (for
#' details see \link{positionMetadata}) either in the form of a position
#' specific or a distance weighted kernel are supported for the spectrum
#' kernel. The generation of an explicit representation is not possible for
#' the position dependent variants of this kernel.\cr\cr
#' Creation of kernel matrix\cr\cr
#' The kernel matrix is created with the function \code{\link{getKernelMatrix}}
#' or via a direct call with the kernel object as shown in the examples below.
#' @return
#' spectrumKernel: upon successful completion, the function returns a kernel
#' object of class \code{\linkS4class{SpectrumKernel}}.
#'
#' @seealso \code{\link{kernelParameters-method}},
#' \code{\link{getKernelMatrix}}, \code{\link{getExRep}},
#' \code{\link{mismatchKernel}}, \code{\link{motifKernel}},
#' \code{\link{gappyPairKernel}}, \code{\linkS4class{SpectrumKernel}}
#' @examples
#'
#' ## instead of user provided sequences in XStringSet format
#' ## for this example a set of DNA sequences is created
#' ## RNA- or AA-sequences can be used as well with the spectrum kernel
#' dnaseqs <- DNAStringSet(c("AGACTTAAGGGACCTGGTCACCACGCTCGGTGAGGGGGACGGGGTGT",
#'                           "ATAAAGGTTGCAGACATCATGTCCTTTTTGTCCCTAATTATTTCAGC",
#'                           "CAGGAATCAGCACAGGCAGGGGCACGGCATCCCAAGACATCTGGGCC",
#'                           "GGACATATACCCACCGTTACGTGTCATACAGGATAGTTCCACTGCCC",
#'                           "ATAAAGGTTGCAGACATCATGTCCTTTTTGTCCCTAATTATTTCAGC"))
#' names(dnaseqs) <- paste("S", 1:length(dnaseqs), sep="")
#'
#' ## create the kernel object for dimers without normalization
#' speck <- spectrumKernel(k=2, normalized=FALSE)
#' ## show details of kernel object
#' speck
#'
#' ## generate the kernel matrix with the kernel object
#' km <- speck(dnaseqs)
#' dim(km)
#' km[1:5,1:5]
#'
#' ## alternative way to generate the kernel matrix
#' km <- getKernelMatrix(speck, dnaseqs)
#' km[1:5,1:5]
#'
#' \dontrun{
#' ## plot heatmap of the kernel matrix
#' heatmap(km, symm=TRUE)
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' (Leslie, 2002) -- C. Leslie, E. Eskin and W.S. Noble. The Spectrum Kernel:
#' A String Kernel for SVM Protein Classification. \cr\cr
#' (Bodenhofer, 2009) -- U. Bodenhofer, K. Schwarzbauer, M. Ionescu and
#' S. Hochreiter. Modelling position specificity in sequence kernels by fuzzy
#' equivalence relations. \cr\cr
#' (Mahrenholz, 2011) -- C.C. Mahrenholz, I.G. Abfalter, U. Bodenhofer, R. Volkmer
#' and S. Hochreiter. Complex networks govern coiled-coil oligomerizations -
#' predicting and profiling by means of a machine learning approach.\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics}, 31(15):2574-2576, 2015.
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.
#' @keywords kernel
#' @keywords spectrum, spectrumKernel
#' @keywords methods
#' @export


spectrumKernel <- function(k=3, r=1, annSpec=FALSE, distWeight=numeric(0),
                     normalized=TRUE, exact=TRUE, ignoreLower=TRUE,
                     presence=FALSE, revComplement=FALSE, mixCoef=numeric(0))
{
    ## check data independent kernel parameters and create closure

    if (!is.numeric(k) || any(k < 1))
        stop("'k' must be larger than 0\n")

    if (!isSingleNumber(r) || r <= 0)
        stop("'r' must be a number larger than 0\n")

    if (!isTRUEorFALSE(normalized))
        stop("'normalized' must be TRUE or FALSE\n")

    if (!isTRUEorFALSE(exact))
        stop("'exact' must be TRUE or FALSE\n")

    if (!isTRUEorFALSE(ignoreLower))
        stop("'ignoreLower' must be TRUE or FALSE\n")

    if (!isTRUEorFALSE(presence))
        stop("'presence' must be TRUE or FALSE\n")

    if (!isTRUEorFALSE(revComplement))
        stop("'revComplement' must be TRUE or FALSE\n")

    if (!isTRUEorFALSE(annSpec))
        stop("'annSpec' must be TRUE or FALSE\n")

    if (length(distWeight) > 0)
    {
        if (!(is.numeric(distWeight) || is.function(distWeight)))
            stop("'distWeight' must be a numeric vector or a function\n")
        
        if (is.function(distWeight))
        {
            func <- deparse(distWeight)[2]
            index <- grep("(", strsplit(func, split="")[[1]], fixed=TRUE,
                          value=FALSE)
            
            if (length(index) < 1)
                stop("Missing parentheses in 'distWeight'\n")
        }
    }

    if (presence && length(distWeight) > 0)
    {
        stop("'presence' can only be used with the position independent\n",
             "       spectrum kernel\n")
    }

    if (length(mixCoef) > 0)
    {
        if (!is.numeric(mixCoef) ||
            length(mixCoef) != k)
            stop("'mixCoef' must be a numeric vector of length k\n")

        if (any(mixCoef < 0))
            stop("mixing coefficients must be non-negative\n")
    }

    if (length(k) == 1)
    {
        ## define function for kernel matrix processing
        rval<- function(x, y = NULL, selx = NULL, sely = NULL, self=NULL)
        {
            return(spectrumProcessing(x=x, y=y, selx=selx, sely=sely, k=k, r=r,
                        annSpec=annSpec, distWeight=distWeight,
                        normalized=normalized, exact=exact,
                        ignoreLower=ignoreLower, presence=presence,
                        revComplement=revComplement, mixCoef=mixCoef,
                        self=self))
        }

        return(new("SpectrumKernel", .Data=rval, .userDefKernel=FALSE,
                   k=k, r=r, normalized=normalized, annSpec=annSpec,
                   distWeight=distWeight, exact=exact,
                   ignoreLower=ignoreLower, presence=presence,
                   revComplement=revComplement, mixCoef=mixCoef))
    }
    else
    {
        ## return list of kernel objects
        kernels <- lapply(k, function(kVal) {
            spectrumKernel(k=kVal,r=r, annSpec=annSpec, distWeight=distWeight,
            normalized=normalized, exact=exact, ignoreLower=ignoreLower,
            presence=presence, revComplement=revComplement, mixCoef=mixCoef)})

        return(kernels)
    }
}

spectrumProcessing <- function(x, y, selx, sely, k, r, annSpec, distWeight,
                               normalized, exact, ignoreLower, presence,
                               revComplement, mixCoef, self=NULL)
{
    if (!is.null(self))
    {
        ## retrieval of kernel parameters
        return(list(k=self@k, r=self@r,
                    normalized=self@normalized,
                    annSpec=self@annSpec,
                    exact=self@exact,
                    ignoreLower=self@ignoreLower,
                    presence=self@presence,
                    distWeight=self@distWeight,
                    revComplement=self@revComplement,
                    mixCoef=self@mixCoef))
    }

    if (missing(x) || is.null(x))
    {
        stop(paste("'x' must be a",
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

    if (length(selx) > 0)
    {
        if (!is.numeric(selx) || length(selx) > length(x))
            stop("'selx' must be an integer vector with indices into 'x'\n")

        selx <- as.integer(selx)
    }
    else
        selx <- 1L:length(x)

    selxC <- selx - 1L
    symmetric <- TRUE

    if (length(k) > 1 || length(r) > 1)
    {
        stop("multiple values for kernel parameter are only allowed\n",
             "        in model selection\n")
    }

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
                   paste(kebabsInfo@allowedSeqClasses, collapse=", "), "\n"))
        }

        if (class(x) != class(y))
            stop("'x' and 'y' must be of matching classes\n")

        if (length(sely) > 0)
        {
            if (!is.numeric(sely) || length(sely) > length(y))
                stop("'sely' must be an integer vector with indices into 'y'\n")

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

    if (annSpec)
    {
        ## limit k to 64 bit feature space
        if (names(bioCharset[[1]]) %in% c("AAexact", "AAiupac") && k > 8)
            stop("'k' must be smaller than or equal to 8\n")
        else if (names(bioCharset[[1]]) %in% c("DNAexact", "RNAexact") && 
                 k > 12)
            stop("for exact charset 'k' must be smaller than or equal to 12\n")
        else if (k > 8)
            stop("for iupac charset 'k' must be smaller than or equal to 8\n")
    }
    else
    {
        ## limit k to 64 bit feature space
        if (names(bioCharset[[1]]) %in% c("AAexact", "AAiupac") && k > 14)
            stop("'k' must be smaller than or equal to 14\n")
        else if (names(bioCharset[[1]]) %in% c("DNAexact", "RNAexact") && 
                 k > 30)
            stop("for exact charset 'k' must be smaller than or equal to 30\n")
        else if (names(bioCharset[[1]]) %in% c("DNAiupac", "RNAiupac") &&
                 k > 15)
            stop("for iupac charset 'k' must be smaller than or equal to 15\n")
    }

    if (!is.null(y))
        seqLength <- c(width(x)[selx], width(y)[sely])
    else
        seqLength <- width(x)[selx]

    maxSeqLength <- max(seqLength)
    posSpec <- FALSE
    offsetX <- integer(0)
    offsetY <- integer(0)

    if (length(distWeight) > 0)
    {
        offsetX <- mcols(x)[["offset"]]

        if (!is.null(offsetX))
        {
            if (!is.integer(offsetX))
                stop("position metadata of 'x' must be an integer vector\n")

            maxDist <- max(width(x[selx]) - offsetX[selx]) -
                       min(-offsetX[selx] + 1)
        }
        else
        {
            offsetX <- integer(0)
            maxDist <- maxSeqLength - 1
        }

        if (!is.null(y))
        {
            offsetY <- mcols(y)[["offset"]]

            if (!is.null(offsetY))
            {
                if (!is.integer(offsetY))
                    stop("position metadata of 'y' must be an integer vector\n")

                if (length(offsetX) > 0)
                {
                    maxDist <- max(c(width(x)[selx] - offsetX[selx],
                                     width(y)[sely] - offsetY[sely])) -
                               min(c(-offsetX[selx] + 1, -offsetY[sely] + 1))
                }
                else
                {
                    maxDist <- max(c(width(x)[selx],
                                     width(y)[sely] - offsetY[sely])) -
                               min(c(1, -offsetY[sely] + 1))
                }
            }
            else
            {
                offsetY <- integer(0)

                if (length(offsetX) > 0)
                {
                    maxDist <- max(c(width(x)[selx] - offsetX[selx],
                                     width(y)[sely])) -
                               min(c(-offsetX[selx] + 1, 1))
                }
                else
                    maxDist <- max(c(width(x)[selx], width(y)[sely])) - 1
            }
        }

        if (is.function(distWeight))
        {
            ## precompute distance weight vector
            ## terminate on stop and warning
            ## assuming that all distances are partially overlapping
            distWeight <- tryCatch(distWeight(0:(maxDist - k + 1)),
                                   warning=function(w) {stop(w)},
                                   error=function(e) {stop(e)})

            if (!(is.numeric(distWeight) && length(distWeight) ==
                  (maxDist - k + 2)))
            {
                stop("distance weighting function did not return a numeric\n",
                     "       vector of correct length\n")
            }

            ## limit to values larger than .Machine$double.eps
            ## for non-monotonic decreasing functions search from end
            for (i in (maxDist - k + 2):1)
            {
                if (distWeight[i] > .Machine$double.eps)
                    break
            }

            distWeight <- distWeight[1:i]
        }

        if (length(distWeight) == 0)
            stop("only zero values for distance weights\n")

        if (isTRUE(all.equal(distWeight, c(1, rep(0, length(distWeight)-1)))))
        {
            posSpec <- TRUE
            distWeight <- numeric(0)
        }
    }

    annCharset <- NULL
    annX <- NULL
    annY <- NULL

    if (annSpec)
    {
        annCharset <- metadata(x)$annotationCharset

        if (is.null(annCharset))
            stop("missing annotation characterset metadata in 'x'\n")

        if (length(annCharset) > 1)
            stop("annotation character set must be one character string\n")

        if (nchar(annCharset[1]) < 2)
        {
            stop("at least two annotation characters needed in annotation\n",
                 "        character set\n")
        }

        annX <- mcols(x)[["annotation"]]

        if (!is.null(y))
            annY <- mcols(y)[["annotation"]]

        if (is.null(annX) || (!is.null(y) && is.null(annY)))
            stop("missing annotation information in 'x' and/or 'y'\n")

        if (!is.character(annX))
            stop("annotation metadata of 'x' must be a character vector\n")

        if (!is.null(y) && !is.character(annY))
            stop("annotation metadata of 'y' must be a character vector\n")
    }

    isXStringSet <- inherits(x, "XStringSet")
    unmapped <- is(x, "DNAStringSet") || is(x, "RNAStringSet")

    if (length(mixCoef) == 0)
    {
        res <- .Call("spectrumKernelMatrixC", x, y, selxC, selyC,
                     as.integer(length(selxC)), as.integer(length(selyC)),
                     as.logical(isXStringSet), as.logical(symmetric),
                     offsetX, offsetY, annCharset, annX, annY,
                     as.integer(bioCharset[[2]]), as.logical(ignoreLower),
                     as.logical(unmapped), as.integer(maxSeqLength),
                     as.integer(k), as.logical(posSpec), distWeight,
                     as.logical(normalized), as.logical(presence),
                     as.logical(revComplement))
    }
    else
    {
        currK <- 0

        if (is.null(y))
            res <- matrix(0, length(selxC), length(selxC))
        else
            res <- matrix(0, length(selxC), length(selyC))

        for (i in 1:k)
        {
            currK <- currK + 1

            if (mixCoef[i] != 0)
            {
                res <- res + mixCoef[i] *
                .Call("spectrumKernelMatrixC", x, y, selxC, selyC,
                      as.integer(length(selxC)), as.integer(length(selyC)),
                      as.logical(isXStringSet), as.logical(symmetric),
                      offsetX, offsetY, annCharset, annX, annY,
                      as.integer(bioCharset[[2]]), as.logical(ignoreLower),
                      as.logical(unmapped), as.integer(maxSeqLength),
                      as.integer(currK), as.logical(posSpec), distWeight,
                      as.logical(normalized), as.logical(presence),
                      as.logical(revComplement))
            }
        }
    }

    if (length(names(x)) > 0)
    {
        rownames(res) <- names(x)[selx]

        if (is.null(y))
            colnames(res) <- rownames(res)
    }

    if (length(names(y)) > 0)
        colnames(res) <- names(y)[sely]

    if (r != 1)
        return(as.KernelMatrix(res^r))
    else
        return(as.KernelMatrix(res))
}

getFeatureSpaceDimension.spectrum <- function(kernel, x)
{
    ## if sequences are not available feature space dimension cannot be
    ## determined
    if (is(x, "KernelMatrix") || is(x, "ExplicitRepresentation"))
        return(-1)

    if (!(class(x) %in% kebabsInfo@allowedSeqSetClasses))
    {
        stop(paste("'x' must be a",
                   paste(kebabsInfo@allowedSeqClasses, collapse=", "), "\n"))
    }

    bioCharset <- getBioCharset(x, kernelParameters(kernel)$exact)
    numAlphaChars <- nchar(bioCharset[[1]])

    ## for IUPAC charactersets the characters - and + are not relevant for
    ## kernels
    if (bioCharset[[2]] %in% c(2,4))
        numAlphaChars <- numAlphaChars - 2

    dimFeatureSpace <- numAlphaChars ^ kernelParameters(kernel)$k

    if (kernelParameters(kernel)$annSpec == TRUE)
    {
        numAnnotChars <- nchar(metadata(x)$annotationCharset)

        if (numAnnotChars == 0)
            stop("Missing annotation in 'x'\n")

        dimFeatureSpace <- dimFeatureSpace *
            numAnnotChars ^ kernelParameters(kernel)$k
    }

    dimFeatureSpace
}

#' @rdname spectrumKernel
#' @aliases
#' getFeatureSpaceDimension
#' getFeatureSpaceDimension,SpectrumKernel-method
#' getFeatureSpaceDimension,ANY-method
#'
#' @param kernel a sequence kernel object
#' @param x one or multiple biological sequences in the form of a
#' \code{\linkS4class{DNAStringSet}}, \code{\linkS4class{RNAStringSet}},
#' \code{\linkS4class{AAStringSet}} (or as \code{\linkS4class{BioVector}})
#' @return of getDimFeatureSpace:
#' dimension of the feature space as numeric value
#' @export
#'

setMethod("getFeatureSpaceDimension",
          signature=signature(kernel="SpectrumKernel"),
          getFeatureSpaceDimension.spectrum)
