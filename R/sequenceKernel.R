#2345678901234567890123456789012345678901234567890123456789012345678901234567890
#' @rdname sequenceKernel
#' @title Sequence Kernel
#' @aliases
#' sequenceKernel
#'
#' @description Create the kernel matrix for a kernel object\cr\cr
#' Retrieve kernel parameters from the kernel object
#'
#' @param kernel one kernel object of class \code{\linkS4class{SequenceKernel}}
#' or one kernlab string kernel (see \code{\link[kernlab:stringdot]{stringdot}}
#' @param x one or multiple biological sequences in the form of a
#' \code{\linkS4class{DNAStringSet}}, \code{\linkS4class{RNAStringSet}},
#' \code{\linkS4class{AAStringSet}} (or as \code{\linkS4class{BioVector}})
#' @param y one or multiple biological sequences in the form of a
#' \code{\linkS4class{DNAStringSet}}, \code{\linkS4class{RNAStringSet}},
#' \code{\linkS4class{AAStringSet}} (or as \code{\linkS4class{BioVector}});
#' if this parameter is specified a rectangular kernel matrix with the
#' samples in \code{x} as rows and the samples in \code{y} as columns is
#' generated otherwise a square kernel matrix with samples in \code{x} as
#' rows and columns is computed; default=\code{NULL}
#' @param selx subset of indices into \code{x}; when this parameter is present
#' the kernel matrix is generated for the specified subset of \code{x} only;
#' default=\code{NULL}
#' @param sely subset of indices into \code{y}; when this parameter is present
#' the kernel matrix is generated for the specified subset of \code{y} only;
#' default=\code{NULL}
#' @details
#' Sequence Kernel\cr\cr
#' A sequence kernel is used for determination of similarity values between
#' biological sequences based on patterns occuring in the sequences. The
#' kernels in this package were specifically written for the biological domain.
#' The corresponding term in the kernlab package is string kernel which is a
#' domain independent implementation of the same functionality which often
#' used in other domains, for example in text classification. For the
#' sequence kernels in this package DNA-, RNA- or AA-acid sequences are used
#' as input with a reduced character set compared to regular text.\cr\cr
#' In string kernels the actual position of a pattern in the sequence/text is
#' irrelevant just the number of occurances of the pattern is important for
#' the similarity consideration. The kernels provided in this package can be
#' created in a position-independent or position-dependent way. Position
#' dependent kernels are using the postion of patterns on the pair of sequences
#' to determine the contribution of a pattern match to the similarity value.
#' For details see help page for \code{\link{positionMetadata}}. As second
#' method of specializing similarity consideration in a kernel is to use
#' annotation information which is placed along the sequences. For details see
#' \code{\link{annotationMetadata}}.
#' Following kernels are available:
#' \itemize{
#' \item{spectrum kernel}
#' \item{mismatch kernel}
#' \item{gappy pair kernel}
#' \item{motif kernel}
#' }
#' These kernels are provided in a position-independent variant. For all
#' kernels except the mismatch also the position-dependent and the
#' annotation-specific variants of the kernel are supported. In addition
#' the spectrum and gappy pair kernel can be created as mixture kernels with
#' the weighted degree kernel and shifted weighted degree kernel being two
#' specific examples of such mixture kernels. The functions described below
#' apply for any kind of kernel in this package.
#' Retrieving kernel paramters from the kernel object\cr\cr
#' The function 'kernelParameters' retrieves the kernel parameters and returns
#' them as list. The function 'seqKernelAsChar' converts a sequnce kernel
#' object into a character string.\cr\cr
#' Generation of kernel matrix\cr\cr
#' The function \code{getKernelMatrix} creates a kernel matrix for the
#' specified kernel and one or two given sets of sequences. It contains
#' similarity values between pairs of samples. If one set of sequences is used
#' the square kernel matrix contains pairwise similarity values for this set.
#' For two sets of sequences the similarities are calculated between these sets
#' resulting in a rectangular kernel matrix. The kernel matrix is always
#' created as dense matrix of the class \code{\linkS4class{KernelMatrix}}.
#' Alternatively the kernel matrix can also be generated via a direct function
#' call with the kernel object. (see examples below)\cr\cr
#' Generation of explicit representation\cr\cr
#' With the function \code{\link{getExRep}} an explicit representation for
#' a specified kernel and a given set of sequences can be generated in sparse
#' or dense form. Applying the linear kernel to the explicit representation
#' with the function \code{\link{linearKernel}} also generates a dense kernel
#' matrix.
#' @return
#' getKernelMatrix: upon successful completion, the function returns a kernel
#' matrix of class \code{\linkS4class{KernelMatrix}} which contains similarity
#' values between pairs of the biological sequences.
#' @seealso \code{\link{as.KernelMatrix}}, \code{\linkS4class{KernelMatrix}},
#' \code{\link{spectrumKernel}}, \code{\link{mismatchKernel}},
#' \code{\link{gappyPairKernel}}, \code{\link{motifKernel}}
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
#' ## create the kernel object with the spectrum kernel
#' spec <- spectrumKernel(k=3, normalized=FALSE)
#'
#' ## generate the kernel matrix
#' km <- getKernelMatrix(spec, dnaseqs)
#' dim(km)
#' km[1:5,1:5]
#'
#' ## alternative way to generate the kernel matrix
#' km <- spec(dnaseqs)
#' km[1:5,1:5]
#'
#' ## generate rectangular kernel matrix
#' km <- getKernelMatrix(spec, x=dnaseqs, selx=1:3, y=dnaseqs, sely=4:5)
#' dim(km)
#' km[1:3,1:2]
#'
#' ## generate a sparse explicit representation
#' er <- getExRep(dnaseqs, spec)
#' er[1:5, 1:8]
#'
#' ## generate kernel matrix from explicit representation
#' km <- linearKernel(er)
#' km[1:5,1:5]
#' 
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}
#' @keywords kernel
#' @keywords methods
#' @export

getKernelMatrix <- function(kernel, x, y, selx, sely)
{
    if (!is(kernel, "stringkernel"))
    {
        return(kernel(x=x, y=y, selx=selx, sely=sely))
    }
    else
    {
        if (missing(y) || is.null(y))
        {
            if (missing(selx) || is.null(selx))
                return(kernelMatrix(kernel=kernel, x=x))
            else
                return(kernelMatrix(kernel=kernel, x=x[selx]))
        }
        else
        {
            if (missing(selx) || is.null(selx))
            {
                if (missing(sely) || is.null(sely))
                    return(kernelMatrix(kernel=kernel, x=x, y=y))
                else
                    return(kernelMatrix(kernel=kernel, x=x, y=y[sely]))
            }
            else
            {
                if (missing(sely) || is.null(sely))
                    return(kernelMatrix(kernel=kernel, x=x[selx], y=y))
                else
                    return(kernelMatrix(kernel=kernel, x=x[selx], y=y[sely]))
            }
        }
    }
}

#' @rdname sequenceKernel
#' @name kernelParameters-method
#' @aliases
#' kernelParameters
#' kernelParameters,SpectrumKernel-method
#' @param object a sequence kernel object
#'
#' @return
#' kernelParameters: the kernel parameters as list
#' @export


setMethod("kernelParameters", signature=signature(object="SpectrumKernel"),
          function(object){object(self=object)})

#' @rdname sequenceKernel
#' @aliases
#' kernelParameters,MismatchKernel-method

setMethod("kernelParameters", signature=signature(object="MismatchKernel"),
          function(object){object(self=object)})

#' @rdname sequenceKernel
#' @aliases
#' kernelParameters,GappyPairKernel-method

setMethod("kernelParameters", signature=signature(object="GappyPairKernel"),
function(object){object(self=object)})

#' @rdname sequenceKernel
#' @aliases
#' kernelParameters,MotifKernel-method

setMethod("kernelParameters", signature=signature(object="MotifKernel"),
          function(object){object(self=object)})

#' @rdname sequenceKernel
#' @aliases
#' kernelParameters,SymmetricPair-method

setMethod("kernelParameters", signature=signature(object="SymmetricPairKernel"),
          function(object){object(self=object)})

## compatibility to kernlab
#setMethod("kernelParameters","kernel", function(object) kpar(object))

getFeatureSpaceDimension.stringkernel <- function(kernel, x)
{
#    library(kernlab)

    if (!(class(x) %in% c("DNAVector", "RNAVector", "AAVector")))
        stop("'x' must be a DNAVector, RNAVector or AAVector\n")

    bioCharset <- getBioCharset(x, TRUE)
    numAlphaChars <- nchar(bioCharset[[1]])

    kernelType <- kpar(kernel)$type
    k <- kpar(kernel)$length
    maxSeqLength <- max(width(x))



    dimFeatureSpace <-
    switch(kernelType,
        "spectrum" = numAlphaChars ^ k,
        "string" = numAlphaChars ^ k,
        "boundrange"  = sum(numAlphaChars ^ (1:k)),
        "exponential"  = sum(numAlphaChars ^ (1:k)),
        "fullstring"  = sum(numAlphaChars ^ (1:k)),
        "constant" = sum(numAlphaChars ^ (1:maxSeqLength))
    )

    return(dimFeatureSpace)
}

## $$$ TODO uncomment
#setMethod("getFeatureSpaceDimension",
#          signature=signature(kernel="stringkernel"),
#          getFeatureSpaceDimension.stringkernel)

getFeatureSpaceDimension.ANY <- function(kernel, x)
{
    ## if kernel is not available dimension cannot be determined
    return(-1)
}

setMethod("getFeatureSpaceDimension",
          signature=signature(kernel="ANY"),
          getFeatureSpaceDimension.ANY)
