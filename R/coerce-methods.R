#2345678901234567890123456789012345678901234567890123456789012345678901234567890
#' @rdname BioVector
#' @param use.names when set to \code{TRUE} the names are preserved
#' @section Coercion methods:
#' In the code snippets below, \code{x} is a \code{BioVector}.
#'
#' \describe{
#'   \item{}{\code{as.character(x, use.names=TRUE)}:
#'   return the sequence set as named or unnamed character vector dependent on
#'   the use.names parameter.
#'   }
#' }
#'

setMethod("as.character", signature(x="BioVector"),
    function(x, use.names=TRUE)
    {
        xc <- x@.Data

        if (use.names && length(x@NAMES) > 0)
            names(xc) <- x@NAMES

        xc
    }
)

asKernelMatrix <- function(x, center=FALSE)
{
    if (center){
        m <- dim(x)[1]
        x <- t(t(x - colSums(x)/m) -  rowSums(x)/m) + sum(x)/m^2
    }

    return(new("KernelMatrix",.Data = x))
}

#' @rdname KernelMatrixAccessors
#' @aliases
#' as.KernelMatrix
#' @param x kernel matrix of class \code{\link{KernelMatrix}}
#' @param center when set to \code{TRUE} the matrix is centered. Default=FALSE
#' @section Coercion methods:
#' In the code snippets below, \code{x} is a kernel matrix.
#'
#' \describe{
#'   \item{}{\code{as.KernelMatrix(x, center=FALSE)}:
#'   center the kernel matrix dependent on the center parameter and coerce
#'   it to class \code{\linkS4class{KernelMatrix}} .
#'   }
#' }

setMethod("as.KernelMatrix", signature(x = "matrix"), asKernelMatrix)

## this method was 33 % faster than the default conversion from Matrix
## for a dgRMatrix to matrix.csr conversion with 15.8 Mio nonzero
## elements in a 10000 x 10000 matrix
## likely has to do with preallocating an empty matrix.csr
## but is commented because of warning in package install for
## character or null from SparseM
##asSparseMCsr.dgRMatrix <- function(from)
##{
##    library(SparseM, warn.conflicts=FALSE)
##
##    ## no dimnames for matrix.csr
##    m1.csr <- as.matrix.csr(matrix(0,1,1))
##    m1.csr@ra <- from@x
##    m1.csr@ja <- as.integer(from@j +1)
##    m1.csr@ia <- as.integer(from@p +1)
##    m1.csr@dimension <- from@Dim
##    return(m1.csr)
##}
#setAs("ExplicitRepresentationSparse", "matrix.csr", asSparseMCsr.dgRMatrix)
#setAs("dgRMatrix", "matrix.csr", asSparseMCsr.dgRMatrix)


asERS.ExplicitRepresentationDense <- function(from)
{
    expRepS <- new("ExplicitRepresentationSparse")
    expRepS@usedKernel <- from@usedKernel
    expRepS@quadratic <- from@quadratic
    smat <- as(as(from@.Data, "RsparseMatrix"), "dgRMatrix")
    expRepS@p <- smat@p
    expRepS@j <- smat@j
    expRepS@Dim <- smat@Dim
    expRepS@Dimnames <- smat@Dimnames
    expRepS@x <- smat@x
    expRepS@factors <- smat@factors
    return(expRepS)
}

setAs("ExplicitRepresentationDense", "ExplicitRepresentationSparse",
      asERS.ExplicitRepresentationDense)

asmatrix.ExplicitRepresentationSparse <- function(from)
{
    smat <- new("dgRMatrix")
    smat@p <- from@p
    smat@j <- from@j
    smat@Dim <- from@Dim
    smat@Dimnames <- from@Dimnames
    smat@x <- from@x
    smat@factors <- from@factors

    return(as(smat, "matrix"))
}

setAs("ExplicitRepresentationSparse", "matrix",
      asmatrix.ExplicitRepresentationSparse)

asERD.ExplicitRepresentationSparse <- function(from)
{
    expRepD <- new("ExplicitRepresentationDense")
    expRepD@usedKernel <- from@usedKernel
    expRepD@quadratic <- from@quadratic
    expRepD@.Data <- as(from, "matrix")
    return(expRepD)
}

setAs("ExplicitRepresentationSparse", "ExplicitRepresentationDense",
      asERD.ExplicitRepresentationSparse)

## fast coercion of ExplicitRepresentationSparse to dgCMatrix together
## with transposition
ersTransposedAsdgCMatrix <- function(from)
{
    if (!is(from, "dgRMatrix"))
        stop("'from' must be a dgRMatrix\n")

    dgC <- new("dgCMatrix")
    dgC@p <- from@p
    dgC@i <- from@j
    dgC@x <- from@x
    dgC@Dim <- as.integer(c(dim(from)[2], dim(from)[1]))
    dgC@Dimnames[1] <- dimnames(from)[2]
    dgC@Dimnames[2] <- dimnames(from)[1]
    return(dgC)
}

distWeightKernelToString <- function(distWeight)
{
    dwLength <- length(distWeight)
    distWeightChar <- ""

    if (is.numeric(distWeight))
    {
        if (dwLength == 1)
            distWeightChar <- distWeight
        else if (dwLength < 7)
        {
            distWeightChar <-
                paste("c(",paste(distWeight,collapse=","),")", sep="")
        }
        else
        {
            distWeightChar <-
                paste("c(", paste(distWeight[1:3],collapse=","), " ... ",
                      paste(distWeight[dwLength - 2:dwLength],collapse=","),
                      ")", sep="")
        }
    }
    else if (typeof(distWeight) == "closure")
    {
        func <- deparse(distWeight)[2]
        index <- grep("(", strsplit(func, split="")[[1]], fixed=TRUE,
                      value=FALSE)

        if (length(index) > 0)
            distWeightChar <- substr(func, 1, index[1] - 1)
        {
            if (distWeightChar %in% c("linWeight", "gaussWeight",
                                      "expWeight"))
            {
                sigma <- get("sigma", envir = environment(distWeight))
                distWeightChar <- paste(distWeightChar, "(sigma=", sigma, ")",
                                        sep="")
            }
            else if (distWeightChar == "swdWeight")
                distWeightChar <- paste(distWeightChar, "()", sep="")
            else
                distWeightChar <- paste(distWeightChar, "( . . . )", sep="")
        }
    }
    distWeightChar
}

#' @rdname sequenceKernel
#' @aliases
#' seqKernelAsChar
#' @param from a sequence kernel object
#'

## coercion is not possible for closure
## in this case call the function directly
seqKernelAsChar <- function(from)
{
    if (!is(from, "SequenceKernel"))
        return(NULL)

    if (isUserDefined(from))
        return(as(from, "character"))

    kChar <- class(from)
    firstPar <- TRUE

    if (!is(from, "MotifKernel"))
    {
        kChar <- paste(kChar, ": k=", kernelParameters(from)$k, sep="")
        firstPar <- FALSE
    }

    if (is(from, "MismatchKernel") || is(from, "GappyPairKernel"))
        kChar <- paste(kChar, ", m=", kernelParameters(from)$m, sep="")

    if (kernelParameters(from)$r != 1)
    {
        if (firstPar)
        {
            kChar <- paste(kChar, ": ", sep="")
            firstPar <- FALSE
        }
        else
            kChar <- paste(kChar, ", ", sep="")

        kChar <- paste(kChar, "r=", kernelParameters(from)$r, sep="")
    }

    if (kernelParameters(from)$annSpec)
    {
        kChar <- paste(kChar, ", annSpec=",
                       kernelParameters(from)$annSpec, sep="")
    }

    if (length(kernelParameters(from)$distWeight) > 0)
    {
        dwString <- distWeightKernelToString(kernelParameters(from)$distWeight)
        kChar <- paste(kChar, ", distWeight=", dwString, sep="")
    }

    if (!kernelParameters(from)$normalized)
    {
        if (firstPar)
        {
            kChar <- paste(kChar, ": ", sep="")
            firstPar <- FALSE
        }
        else
            kChar <- paste(kChar, ", ", sep="")

        kChar <- paste(kChar, "normalized=", kernelParameters(from)$normalized,
                       sep="")
    }

    if (!kernelParameters(from)$exact)
    {
        if (firstPar)
        {
            kChar <- paste(kChar, ": ", sep="")
            firstPar <- FALSE
        }
        else
            kChar <- paste(kChar, ", ", sep="")

        kChar <- paste(kChar, "exact=", kernelParameters(from)$exact, sep="")
    }

    if (!kernelParameters(from)$ignoreLower)
    {
        if (firstPar)
        {
            kChar <- paste(kChar, ": ", sep="")
            firstPar <- FALSE
        }
        else
            kChar <- paste(kChar, ", ", sep="")

        kChar <- paste(kChar, "ignoreLower=",
                       kernelParameters(from)$ignoreLower, sep="")
    }

    if (kernelParameters(from)$presence)
    {
        if (firstPar)
        {
            kChar <- paste(kChar, ": ", sep="")
            firstPar <- FALSE
        }
        else
            kChar <- paste(kChar, ", ", sep="")

        kChar <- paste(kChar, "presence=",
                       kernelParameters(from)$presence, sep="")
    }

    return(kChar)
}

