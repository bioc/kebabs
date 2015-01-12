#2345678901234567890123456789012345678901234567890123456789012345678901234567890
#' @rdname motifKernel
#' @title Motif Kernel
#'
#' @description Create a motif kernel object and the kernel matrix
#'
#' @param motifs a set of motif patterns specified as character vector. The
#' order in which the patterns are passed for creation of the kernel object
#' also determines the order of the features in the explicit representation.
#' Lowercase characters in motifs are always converted to uppercase. For
#' details concerning the definition of motif patterns see below and in the
#' examples section.
#'
#' @param r exponent which must be > 0 (see details section in
#' \link{spectrumKernel}). Default=1
#'
#' @param annSpec boolean that indicates whether sequence annotation should
#' be taken into account (details see on help page for
#' \code{\link{annotationMetadata}}). Default=FALSE
#'
#' @param distWeight a numeric distance weight vector or a distance weighting
#' function (details see on help page for \code{\link{gaussWeight}}).
#' Default=NULL
#'
#' @param normalized generated data from this kernel will be normalized
#' (details see below). Default=TRUE
#'
#' @param exact use exact character set for the evaluation (details see below).
#' Default=TRUE
#'
#' @param ignoreLower ignore lower case characters in the sequence. If the
#' parameter is not set lower case characters are treated like uppercase.
#' default=TRUE
#'
#' @param presence if this parameter is set only the presence of a motif will
#' be considered, otherwise the number of occurances of the motif is used;
#' Default=FALSE
#'
#' @details
#' Creation of kernel object\cr\cr
#' The function 'motif' creates a kernel object for the motif kernel for a set
#' of given DNA-, RNA- or AA-motifs. This kernel object can then be used to
#' generate a kernel matrix or an explicit representation for this kernel.
#' The individual patterns in the set of motifs are built similar to regular
#' expressions through concatination of following elements in arbitrary order:
#' \itemize{
#' \item{a specific character from the used character set - e.g. 'A' or 'G' in
#'   DNA patterns for matching a specific character}
#' \item{the wildcard character '.' which matches any valid character of the
#'   character set except '-'}
#' \item{a substitution group specified by a collection of characters from the
#'   character set enclosed in square brackets - e.g. [AG] - which matches any
#'   of the listed characters; with a leading '^' the character list is
#'   inverted and matching occurs for all characters of the character set
#'   which are not listed except '-'}
#' }
#' For values different from 1 (=default value) parameter \code{r} leads
#' to a transfomation of similarities by taking each element of the similarity
#' matrix to the power of r. For the annotation specific variant of this
#' kernel see \link{annotationMetadata}, for the distance weighted
#' variants see \link{positionMetadata}. If \code{normalized=TRUE}, the
#' feature vectors are scaled to the unit sphere before computing the
#' similarity value for the kernel matrix. For two samples with the feature
#' vectors \code{x} and \code{y} the similarity is computed as:
#' \deqn{s=\frac{\vec{x}^T\vec{y}}{\|\vec{x}\|\|\vec{y}\|}}{s=(x^T y)/(|x| |y|)}
#' For an explicit representation generated with the feature map of a
#' normalized kernel the rows are normalized by dividing them through their
#' Euclidean norm. For parameter \code{exact=TRUE} the sequence characters
#' are interpreted according to an exact character set. If the flag is not
#' set ambigous characters from the IUPAC characterset are also evaluated.
#'
#' The annotation specific variant (for details see
#' \link{annotationMetadata}) and the position dependent variants (for
#' details see \link{positionMetadata}) either in the form of a position
#' specific or a distance weighted kernel are supported for the motif kernel.
#' The generation of an explicit representation is not possible for the
#' position dependent variants of this kernel.
#'
#' Hint: For a normalized motif kernel with a feature subset of a normalized
#' spectrum kernel the explicit representation will not be identical to the
#' subset of an explicit representation for the spectrum kernel because
#' the motif kernel is not aware of the other kmers which are used in the
#' spectrum kernel additionally for normalization.\cr\cr
#' Creation of kernel matrix\cr\cr
#' The kernel matrix is created with the function \code{\link{getKernelMatrix}}
#' or via a direct call with the kernel object as shown in the examples below.
#'
#' @return
#' motif: upon successful completion, the function returns a kernel
#' object of class \code{\linkS4class{MotifKernel}}.
#' @seealso \code{\link{kernelParameters-method}},
#' \code{\link{getKernelMatrix}}, \code{\link{getExRep}},
#' \code{\link{spectrumKernel}}, \code{\link{mismatchKernel}},
#' \code{\link{gappyPairKernel}}
#' @examples
#'
#' ## instead of user provided sequences in XStringSet format
#' ## for this example a set of DNA sequences is created
#' ## RNA- or AA-sequences can be used as well with the motif kernel
#' dnaseqs <- DNAStringSet(c("AGACTTAAGGGACCTGGTCACCACGCTCGGTGAGGGGGACGGGGTGT",
#'                           "ATAAAGGTTGCAGACATCATGTCCTTTTTGTCCCTAATTATTTCAGC",
#'                           "CAGGAATCAGCACAGGCAGGGGCACGGCATCCCAAGACATCTGGGCC",
#'                           "GGACATATACCCACCGTTACGTGTCATACAGGATAGTTCCACTGCCC",
#'                           "ATAAAGGTTGCAGACATCATGTCCTTTTTGTCCCTAATTATTTCAGC"))
#' names(dnaseqs) <- paste("S", 1:length(dnaseqs), sep="")
#'
#' ## create the kernel object with the motif patterns
#' mot <- motifKernel(c("A[CG]T","C.G","G[^A][AT]"), normalized=FALSE)
#' ## show details of kernel object
#' mot
#'
#' ## generate the kernel matrix with the kernel object
#' km <- mot(dnaseqs)
#' dim(km)
#' km
#'
#' ## alternative way to generate the kernel matrix
#' km <- getKernelMatrix(mot, dnaseqs)
#'
#' \dontrun{
#' ## plot heatmap of the kernel matrix
#' heatmap(km, symm=TRUE)
#'
#' ## generate rectangular kernel matrix
#' km <- mot(x=dnaseqs, selx=1:3, y=dnaseqs, sely=4:5)
#' dim(km)
#' km
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' (Ben-Hur, 2003) -- A. Ben-Hur, and D. Brutlag. Remote homology detection:
#' a motif based approach. \cr\cr
#' (Bodenhofer, 2009) -- U. Bodenhofer, K. Schwarzbauer, M. Ionescu and
#' S. Hochreiter. Modelling position specificity in sequence kernels by fuzzy
#' equivalence relations. \cr
#' (Mahrenholz, 2011) -- C. Mahrenholz, I. Abfalter, U. Bodenhofer, R. Volkmer
#' and S. Hochreiter. Complex networks govern coiled-coil oligomerizations -
#' predicting and profiling by means of a machine learning approach.\cr\cr
#' \url{http://www.bioinf.jku.at/software/kebabs}
#' @keywords kernel
#' @keywords motifKernel, motif
#' @keywords methods
#' @export


motifKernel <- function(motifs, r=1, annSpec=FALSE, distWeight=numeric(0),
                  normalized=TRUE, exact=TRUE, ignoreLower=TRUE,
                  presence=FALSE)
{
    ## define function for kernel processing

    rval <- function(x, y = NULL, selx = NULL, sely = NULL, self=NULL)
    {
        return(motifProcessing(x=x, y=y, selx=selx, sely=sely, motifs=motifs,
                               motifLengths=motifLengths, r=r, annSpec=annSpec,
                               distWeight=distWeight, normalized=normalized,
                               exact=exact, ignoreLower=ignoreLower,
                               presence=presence, self=self))
    }

    ## check data independent kernel parameters and create closure

    if (!is.character(motifs) || length(motifs) < 1)
        stop("please specify motifs as character vector\n")

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

    if (!isTRUEorFALSE(annSpec))
        stop("annSpec must be TRUE or FALSE\n")

    if (length(distWeight) > 0 &&
        !(is.numeric(distWeight) || is.function(distWeight)))
        stop("distWeight must be a numeric vector or a function\n")

    if (presence && length(distWeight) > 0)
    {
        stop("presence can only be used with the position independent\n",
             "       motif kernel\n")
    }

    motifs <- toupper(motifs)

    if (length(motifs) != length(unique(motifs)))
        stop("motifs are not unique\n")

    motifLengths <- nchar(motifs)
    result <- .Call("validateMotifsC", motifs, motifLengths)

    if (result == FALSE)
    {
        message("motifs not defined correctly - errors occurred in\n",
                "motif at position: (first 25 errors only)\n")

        wrongMotifs <- 0

        for (i in 1:length(motifs))
        {
            if (motifLengths[i] != -1)
            {
                if (wrongMotifs > 25)
                    break
                else
                {
                    wrongMotifs <- wrongMotifs + 1
                    message(i, "        ", motifLengths[i])
                }
            }
        }

        message("")
        stop("please correct motif definitions\n")
    }

    ## motifs are stored in environment only which is not copied
    return(new("MotifKernel", .Data=rval, r=r,
               normalized=normalized, annSpec=annSpec,
               distWeight=distWeight, exact=exact,
               ignoreLower=ignoreLower, presence=presence))
}

motifProcessing <- function(x, y, selx, sely, motifs, motifLengths, r,
                            annSpec, distWeight, normalized, exact,
                            ignoreLower, presence, self=NULL)
{
    if (!is.null(self))
    {
        ## retrieval of kernel parameters
        return(list(motifs=motifs, r=self@r,
                    motifLengths=motifLengths,
                    normalized=self@normalized,
                    annSpec=self@annSpec,
                    exact=self@exact,
                    ignoreLower=self@ignoreLower,
                    presence=self@presence,
                    distWeight=self@distWeight,
                    revComplement=self@revComplement))
    }

    if (missing(x) || is.null(x))
    {
        stop(paste("x must be a",
                   paste(kebabsInfo@allowedSeqClasses, collapse=", ")))
    }

    if (length(x) < 1)
        stop("sequence info is missing in motif kernel\n")

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
                   paste(kebabsInfo@allowedSeqClasses, collapse=", "),"\n"))
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
                       paste(kebabsInfo@allowedSeqClasses, collapse=", "),"\n"))
        }

        if (class(x) != class(y))
            stop("x and y must be of matching classes\n")

        if (length(sely) > 0)
        {
            if (!is.numeric(sely) || length(sely) > length(y))
                stop("sely must be an integer vector with indices into y\n")

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

    ## feature space is limited by no of motifs

    if (!is.null(y))
        seqLength <- c(width(x)[selx], width(y)[sely])
    else
        seqLength <- width(x)[selx]

    maxSeqLength <- max(seqLength)
    maxMotifLength <- max(motifLengths)
    maxPatternLength <- max(nchar(motifs))

    if (any(seqLength < min(motifLengths)))
    {
        stop("motif kernel does not accept strings shorter\n",
             "       than the minimal motif length\n")
    }

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

            maxDist <- max(abs(offsetX), abs(width(x) - offsetX))
        }
        else
        {
            offsetX <- integer(0)
            maxDist <- maxSeqLength
        }

        if (!is.null(y))
        {
            offsetY <- mcols(y)[["offset"]]

            if (!is.null(offsetY))
            {
                if (!is.integer(offsetY))
                    stop("position metadata of 'y' must be an integer vector\n")

                maxDistY <- max(abs(offsetY), abs(width(y) - offsetY))
                maxDist <- max(maxDist, maxDistY)
            }
            else
                offsetY <- integer(0)
        }

        if (is.function(distWeight))
        {
            ## precompute distance weight vector
            ## terminate on stop and warning
            ## assuming that all distances are partially overlapping
            distWeight <- tryCatch(distWeight(0:(2 * maxDist)),
                                   warning=function(w) {stop(w)},
                                   error=function(e) {stop(e)})

            if (!(is.numeric(distWeight) && length(distWeight) ==
                  2 * maxDist + 1))
            {
                stop("distWeight function did not return a numeric vector\n",
                     "       of correct length\n")
            }

            ## limit to values larger than .Machine$double.eps
            ## for non-monotonic decreasing functions search from end
            for (i in (2 * maxDist):0)
            {
                if (distWeight[i] > .Machine$double.eps)
                    break
            }

            distWeight <- distWeight[0:i]
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
            stop("missing annotation characterset metadata in x\n")

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
            stop("missing annotation information in x and/or y\n")

        if (!is.character(annX))
            stop("annotation metadata of 'x' must be a character vector\n")

        if (!is.null(y) && !is.character(annY))
            stop("annotation metadata of 'y' must be a character vector\n")

        ## set exit hook for C heap cleanup
        on.exit(.Call("freeHeapCallocsC", 4L))
    }

    ## rough limit for no of nodes in motif tree from no of chars and
    ## no of substitution groups, add one for root
    nodeLimit <- sum(motifLengths) + 1 +
                 sum(sapply(gregexpr("[", motifs, fixed=TRUE),
                            function(x) length(unlist(x))))
    isXStringSet <- inherits(x, "XStringSet")
    unmapped <- is(x, "DNAStringSet") || is(x, "RNAStringSet")

    res <- .Call("motifKernelMatrixC", x, y, selxC, selyC,
                 as.integer(length(selxC)), as.integer(length(selyC)),
                 as.logical(isXStringSet), as.logical(symmetric), offsetX,
                 offsetY, annCharset, annX, annY, motifs, motifLengths,
                 as.integer(nodeLimit), as.integer(maxMotifLength),
                 as.integer(maxPatternLength), as.integer(bioCharset[[2]]),
                 as.logical(ignoreLower), as.logical(unmapped),
                 as.integer(max(seqLength)), as.logical(posSpec),
                 distWeight, as.logical(normalized), as.logical(presence))

    if (length(names(x)) > 0)
    {
        rownames(res) <- names(x)[selx]

        if (is.null(y))
            colnames(res) <- names(x)[selx]
    }

    if (length(names(y)) > 0)
    {
        colnames(res) <- names(y)[sely]
    }

    if (r != 1)
        return(as.KernelMatrix(res^r))
    else
        return(as.KernelMatrix(res))
}

getFeatureSpaceDimension.motif <- function(kernel, x)
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

    dimFeatureSpace <- length(kernelParameters(kernel)$motifs)

    if (kernelParameters(kernel)$annSpec == TRUE)
    {
        numAnnotChars <- nchar(metadata(x)$annotationCharset)

        if (numAnnotChars == 0)
            stop("Missing annotation in 'x'\n")

        dimFeatureSpace <- dimFeatureSpace *
            numAnnotChars ^ max(kernelParameters(kernel)$motifLengths)
    }

    dimFeatureSpace
}

#' @rdname motifKernel
#' @aliases
#' getFeatureSpaceDimension,MotifKernel-method
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
          signature=signature(kernel="MotifKernel"),
          getFeatureSpaceDimension.motif)
