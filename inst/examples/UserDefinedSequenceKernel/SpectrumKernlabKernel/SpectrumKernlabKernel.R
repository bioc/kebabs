#2345678901234567890123456789012345678901234567890123456789012345678901234567890
## Bioconductor requires an 80 column layout of R sources

## kernel class definition
setClass("SpectrumKernlabKernel",
    representation = representation
    (
        k               = "numeric",
        r               = "numeric",
        normalized      = "logical"
    ),
    prototype = prototype
    (
        k               = 3,
        r               = 1,
        normalized      = TRUE
    ),
    contains="SequenceKernel"
)

## coercion to character
asChar.SpectrumKernlabKernel <- function(from)
{
    kChar <- paste("Spectrum Kernlab Kernel: k=", from@k, sep="")
    if (from@r > 1)
        cat(paste(", r=", object@r, sep=""))
    if (!from@normalized)
        kChar <- paste(kChar, ", normalized=FALSE", sep="")
    return(kChar)
}

## set up coercion method for class SpectrumKernlabKernel
setAs("SpectrumKernlabKernel", "character", asChar.SpectrumKernlabKernel)

## show method to display the kernel object
setMethod("show", signature(object="SpectrumKernlabKernel"),
    function(object)
    {
        cat(as(object, "character"))
        cat("\n")
    }
)

## constructor to create a kernel object
spectrumKernlabKernel <- function(k=3, r=1, normalized=TRUE)
{
    ## check data independent kernel parameters and create closure
    
    if (!is.numeric(k) || any(k < 1))
        stop("'k' must be larger than 0\n")

    if (!isSingleNumber(r) || r <= 0)
        stop("'r' must be a number larger than 0\n")

    if (!isTRUEorFALSE(normalized))
        stop("'normalized' must be TRUE or FALSE\n")

    if (length(k) == 1)
    {
        ## define function for kernel matrix processing
        rval<- function(x, y = NULL, selx = NULL, sely = NULL, self=NULL)
        {
            return(spectrumKernlabProcessing(x=x, y=y, selx=selx, sely=sely,
                                             k=k, r=r, normalized=normalized,
                                             self=self))
        }

        return(new("SpectrumKernlabKernel", .Data=rval, k=k, r=r,
                   normalized=normalized))
    }
    else
    {
        ## return list of kernel objects
        kernels <- lapply(k, function(kVal) {
            spectrumKernlabKernel(k=kVal,r=r, normalized=normalized)})

        return(kernels)
    }
}

## function to compute the kernel matrix for one or two sequence sets
## subsetting of sequence sets happens with parameters selx and sely
## parameter self is only needed for call from method kernelParameters

spectrumKernlabProcessing <- function(x, y, selx, sely, k, r, normalized,
                                      self=NULL)
{
    if (!is.null(self))
    {
        ## retrieval of kernel parameters
        return(list(k=self@k, r=self@r,
        normalized=self@normalized))
    }

    if (missing(x) || is.null(x))
    {
        stop(paste("'x' must be a",
        paste(kebabs:::kebabsInfo@allowedSeqClasses, collapse=", "), "\n"))
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

    if (!(class(x) %in% kebabs:::kebabsInfo@allowedSeqSetClasses))
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

        if (!(class(y) %in% kebabs:::kebabsInfo@allowedSeqSetClasses))
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
    }
    else
        sely <- integer(0)

    ## user-defined logic to compute the kernel matrix starts here
    
    ## If the kernel needs C/C++ code implemented in a shared library
    ## the library must be preloaded with dyn.load() before the R file
    ## with the kernel definition interfacing to the shared library is
    ## sourced. The shared C/C++ library can be built from source on
    ## the commandline with R CMD SHLIB.
    
    ## For usage with the SVMs supported in kebabs the user-defined kernel
    ## must be positive definite. Because a positive definite kernel is
    ## symmetric for a quadratic kernel matrix only the upper triangular
    ## matrix should be computed, the other values can be set by simple
    ## mirroring on the main matrix diagonal.
    
    ## In this example the kernlab stringdot spectrum kernel
    ## is used instead of a user defined-kernel logic.
    
    ## create the kernlab kernel object
    sk <- stringdot(length=k, normalized=normalized)

    ## invoke the kernlab stringdot kernel as user-defined kernel processing
    if (is(x, "BioVector"))
    {
        if (is.null(y))
            res <- kernelMatrix(sk, x[selx])
        else
            res <- kernelMatrix(sk, x[selx], y[sely])
    }
    else
    {
        if (is.null(y))
            res <- kernelMatrix(sk, as.character(x[selx]))
        else
        {
            res <- kernelMatrix(sk, as.character(x[selx]),
                                as.character(y[sely]))
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

## method to get kernel parameters from kernel object
setMethod("kernelParameters", signature=signature(object="SpectrumKernlabKernel"),
          function(object){object(self=object)})
