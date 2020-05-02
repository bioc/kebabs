##345678901234567890123456789012345678901234567890123456789012345678901234567890
#' @rdname annotationSpecificKernel
#' @title Annotation Specific Kernel
#' @aliases
#' annotationSpecificKernel
#' AnnotationSpecificKernel
#'
#' @description Assign annotation metadata to sequences and create a kernel
#' object which evaluates annotation information\cr\cr
#' Show biological sequence together with annotation\cr\cr
#'
#' @details
#' Annotation information for sequences\cr\cr
#' For the annotation specific kernel additional annotation information is
#' added to the sequence data. The annotation for one sequence consist of a
#' character string with a single annotation character per position, i.e.
#' the annotation sequence has the same length as the sequence. The character
#' set used for annotation is defined user specific on XStringSet level
#' with up to 32 different characters. Each biological sequence needs
#' an associated annotation sequence assigned consisting of characters from
#' this character set. The evaluation of annotation information as part of
#' the kernel processing during generation of a kernel matrix or an explict
#' representation can be activated per kernel object. \cr\cr
#' Assignment of annotation information\cr\cr
#' The annotation characterset consists of a character string listing all
#' allowed annotation characters in alphabetical order. Any single byte ASCII
#' character from the decimal range between 32 and 126, except 45, is allowed.
#' The character '-' (ASCII dec. 45) is used for masking sequence parts which
#' should not be evaluated. As it has assigned this special masking function
#' it must not be used in annotation charactersets.\cr\cr
#' The annotation characterset is assigned to the sequence set with the
#' \code{annotationMetadata} function (see below). It is stored in the
#' metadata list as named element \code{annotationCharset} and can be stored
#' along with other metadata assigned to the sequence set. The annotation
#' strings for the individual sequences are represented as a character vector
#' and can be assigned to the XStringSet together with the assignment of the
#' annotation characterset as element related metadata. Element related
#' metadata is stored in a DataFrame and the columns of this data frame
#' represent the different types of metadata that can be assigned in parallel.
#' The column name for the sequence related annotation information is
#' "annotation". (see Example section for an example of annotation metadata
#' assignment) Annotation metadata can be assigned together with position
#' metadata (see \code{\link{positionMetadata}} to a sequence set. \cr\cr
#' Annotation Specific Kernel Processing\cr\cr
#' The annotation specific kernel variant of a kernel, e.g. the spectrum kernel
#' appends the annotation characters corresponding to a specific kmer to this
#' kmer and treats the resulting pattern as one feature - the basic unit for
#' similarity determination. The full feature space of an annotation specific
#' spectrum kernel is the cartesian product of the set of all possible sequence
#' patterns with the set of all possible anntotions patterns. Dependent on the
#' number of characters in the annotation character set the feature space
#' increases drastically compared to the normal spectrum kernel. But through
#' annotation the similarity consideration between two sequences can be split
#' into independent parts considered separately, e.g. coding/non-coding,
#' exon/intron, etc... . For amino acid sequences e.g. a heptad annotation
#' (consisting of a usually periodic pattern of 7 characters (a to g) can be
#' used as annotation like in prediction of coiled coil structures. (see
#' reference Mahrenholz, 2011)\cr\cr
#' The flag \code{annSpec} passed during creation of a kernel object controls
#' whether annotation information is evaluated by the kernel. (see functions
#' \code{\link{spectrumKernel}, \link{gappyPairKernel}, \link{motifKernel}})
#' In this way sequences with annotation can be evaluated annotation specific
#' and without annotation through using two different kernel objects. (see
#' examples below) The annotation specific kernel variant is available for all
#' kernels in this package except for the mismatch kernel.\cr\cr
#' annotationMetadata function \cr\cr
#' With this function annotation metadata can be assigned to sequences defined
#' as XStringSet (or BioVector). The sequence annotation strings are stored
#' as element related information and can be retrieved with the method
#' \code{\link{mcols}}. The characters used for anntation are stored as
#' annotation characterset for the sequence set and can be retrieved
#' with the method \code{\link{metadata}}. For the assignment of annotation
#' metadata to biological sequences this function should be used instead of the
#' lower level functions metadata and mcols. The function
#' \code{annotationMetadata} performs several checks and also takes care
#' that other metadata or element metadata assigned to the object is kept.
#' Annotation metadata are deleted if the parameters \code{annCharset} and
#' \code{annotation} are set to NULL.\cr\cr
#' showAnnotatedSeq function \cr\cr
#' This function displays individual sequences aligned with the annotation
#' string with 50 positions per line. The two header lines show the start
#' postion for each bock of 10 characters.\cr\cr
#'
#' @examples
#'
#' ## create a set of annotated DNA sequences
#' ## instead of user provided sequences in XStringSet format
#' ## for this example a set of DNA sequences is created
#' x <- DNAStringSet(c("AGACTTAAGGGACCTGGTCACCACGCTCGGTGAGGGGGACGGGGTGT",
#'                     "ATAAAGGTTGCAGACATCATGTCCTTTTTGTCCCTAATTATTTCAGC",
#'                     "CAGGAATCAGCACAGGCAGGGGCACGGCATCCCAAGACATCTGGGCC",
#'                     "GGACATATACCCACCGTTACGTGTCATACAGGATAGTTCCACTGCCC",
#'                     "ATAAAGGTTGCAGACATCATGTCCTTTTTGTCCCTAATTATTTCAGC"))
#' names(x) <- paste("S", 1:length(x), sep="")
#' ## define the character set used in annotation
#' ## the masking character '-' is is not part of the character set
#' anncs <- "ei"
#' ## annotation strings for each sequence as character vector
#' ## in the third and fourth sample a part of the sequence is masked
#' annotStrings <- c("eeeeeeeeeeeeiiiiiiiiieeeeeeeeeeeeeeeeiiiiiiiiii",
#'                   "eeeeeeeeeiiiiiiiiiiiiiiiiiiieeeeeeeeeeeeeeeeeee",
#'                   "---------eeeeeeeeeeeeeeeeiiiiiiiiiiiiiiiiiiiiii",
#'                   "eeeeeeeeeeeeeeeeeeeeeeeiiiiiiiiiiiiiiiiiiii----",
#'                   "eeeeeeeeeeeeiiiiiiiiiiiiiiiiiiiiiiieeeeeeeeeeee")
#' ## assign metadata to DNAString object
#' annotationMetadata(x, annCharset=anncs) <- annotStrings
#' ## show annotation
#' annotationMetadata(x)
#' annotationCharset(x)
#'
#' ## show sequence 3 aligned with annotation string
#' showAnnotatedSeq(x, sel=3)
#'
#' ## create annotation specific spectrum kernel
#' speca <- spectrumKernel(k=3, annSpec=TRUE, normalized=FALSE)
#'
#' ## show details of kernel object
#' kernelParameters(speca)
#'
#' ## this kernel object can be now be used in a classification or regression
#' ## task in the usual way or you can use the kernel for example to generate
#' ## the kernel matrix for use with another learning method in another R
#' ## package.
#' kma <- speca(x)
#' kma[1:5,1:5]
#' ## generate a dense explicit representation for annotation-specific kernel
#' era <- getExRep(x, speca, sparse=FALSE)
#' era[1:5,1:8]
#'
#' ## when a standard spectrum kernel is used with annotated
#' ## sequences the anntotation information is not evaluated
#' spec <- spectrumKernel(k=3, normalized=FALSE)
#' km <- spec(x)
#' km[1:5,1:5]
#'
#' ## finally delete annotation metadata if no longer needed
#' annotationMetadata(x) <- NULL
#' ## show empty metadata
#' annotationMetadata(x)
#' annotationCharset(x)
#' 
#' @seealso \code{\link{spectrumKernel}}, \code{\link{gappyPairKernel}},
#' \code{\link{motifKernel}}, \code{\link{positionMetadata}},
#' \code{\link{metadata}}, \code{\link{mcols}}
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' C.C. Mahrenholz, I.G. Abfalter, U. Bodenhofer, R. Volkmer and
#' S. Hochreiter (2011) Complex networks govern coiled coil
#' oligomerization - predicting and profiling by means of a machine
#' learning approach. \emph{Mol. Cell. Proteomics}.
#' DOI: 10.1074/mcp.M110.004994. \cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics}, 31(15):2574-2576, 2015.
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.
#' @keywords kernel
#' @keywords annotation
#' @keywords methods



#' @rdname annotationSpecificKernel
#'
#' @param x biological sequences in the form of a
#' \code{\linkS4class{DNAStringSet}}, \code{\linkS4class{RNAStringSet}},
#' \code{\linkS4class{AAStringSet}} (or as \code{\linkS4class{BioVector}})
#' @param sel single index into x for displaying a specific sequence. Default=1
#' @param ann show annotation information along with the sequence
#' @param pos show position information
#' @param start first postion to be displayed, by default the full sequence is
#' shown
#' @param end last position to be displayed or use parameter 'width'
#' @param width number of positions to be displayed or use parameter 'end'
#' @export


## $$$ TODO position specific handling is still missing
showAnnotatedSeq <- function(x, sel=1, ann=TRUE, pos=TRUE, start=1,
                             end=width(x)[sel], width=NA)
{
    charsPerLine <- 50L

    if (!(is(x, "XStringSet") || is(x, "BioVector")))
        stop("sequence is not an XStringSet or BioVector\n")

    if (!isSingleNumber(sel))
        stop("only one sequence can be shown with annotation and position\n")

    if (!(isSingleNumber(start) && isSingleNumber(end)))
        stop("'start' and 'end' must be single integer values\n")

    if (!missing(width))
    {
        if (!isSingleNumber(width) || width <=0)
            stop("'width' must be a positive integer value\n")

        end <- start + width - 1
    }

    if (end > (start + width(x)[sel] - 1))
        end <- start + width(x)[sel] - 1

    if (start != 1 || end != width(x)[sel])
        seq <- subseq(x[sel], start=start, end=end)
    else
        seq <- x[sel]

    cat("\nSequence", names(x[sel]), "of class", class(x), "\n\n")

    if (ann)
    {
        if (length(metadata(x)) == 0 || is.null(mcols(x)) ||
            is.null(metadata(x)$annotationCharset) ||
            is.null(mcols(x)$annotation))
        {
            cat("no annotation available\n\n")
            ann <- FALSE
        }
        else
        {
            cat("Annotation character set:",
                metadata(x)$annotationCharset, "\n\n")
            annot <- mcols(x)$annotation[sel]
            annot <- substr(annot, start=start, stop=end)
        }
    }

    numChars <- end - start + 1
    widthIndex <- ceiling(log10(numChars)) + 2

    for (i in (1:((numChars + 49) / 50)))
    {
        if (pos==TRUE && ((i-1) %% 10 == 0))
        {
            cat(paste(rep(" ", widthIndex+3), collapse=""),
                paste(0:4, "        ", sep=""), "\n")
            cat(paste(rep(" ", widthIndex+3), collapse=""),
                paste(rep(1,5), "        ", sep=""), "\n")
        }

        cat(format(paste("[", (i-1) * 50 + start, "]", sep=""),
                   width=widthIndex, justify="right"),
            "  ", substring(seq, (i-1) * 50 + 1, i * 50), "\n")

        if (ann)
        {
            cat(paste(rep(" ", widthIndex+3), collapse=""),
                substring(annot, (i-1) * 50 + 1, i * 50), "\n")
        }
    }
    cat("\n")
}

setAnnotationMetadata.seq <- function(x, ... , value)
{
   annCharset <- list(...)$annCharset
   if ((is.null(annCharset) && !is.null(value)) ||
       (!is.null(annCharset) && is.null(value)))
       stop("to reset annotation metadata set 'value' to NULL\n", call.=FALSE)

#    if (is.null(annCharset) && !is.null(value))
#       stop("parameter 'annCharset' is missing\n", call.=FALSE)
#
#    if (!is.null(annCharset) && is.null(value))
#        stop("to reset annotation metadata set 'value' to NULL\n",
#             "        the annotation characterset is reset too\n", 
#             call.=FALSE)

    if (is.null(annCharset) && is.null(value))
    {
        ## reset annotation metadata
        md <- metadata(x)
        md["annotationCharset"] <- NULL
        metadata(x) <- md

        ## reset element metadata
        emd <- mcols(x)

        if (is(emd, "DataFrame_OR_NULL"))
        {
            annCol <- which(colnames(emd) == "annotation")

            if (length(annCol == 1))
            {
                emd <- emd[, -annCol, drop=FALSE]

                if (ncol(emd) == 0)
                    emd <- NULL

                mcols(x) <- emd
            }
        }

        return(x)
    }

    if (!(is.character(annCharset) && length(annCharset) == 1 &&
          nchar(annCharset)[1] > 0))
        stop("'annCharset' must be a single character string\n", call.=FALSE)

    singleChars <- unique(unlist(strsplit(annCharset, split="")))

    if(length(singleChars) > 32)
        stop("only up to 32 annotation characters are allowed\n", call.=FALSE)

    if (length(singleChars) != nchar(annCharset)[1])
    {
        stop("characters in annotation character set are not unique\n",
             call.=FALSE)
    }

    dotPresent <- which(singleChars == ".")

    if (length(dotPresent) > 0)
        stop("'.' is not allowed in annotation character set\n")

    maskingChar <- which(singleChars == "-")

    if (length(maskingChar) > 0)
        singleChars <- singleChars[-maskingChar]

    singleChars <- sort(singleChars)
    annCharset <- paste(singleChars, collapse="")

    if (!is.character(value))
    {
        ## try to convert to character vector
        value <- tryCatch(as.character(value),
                               warning=function(w) {stop(w)},
                               error=function(e) {stop(e)})
    }

    if (length(x) != length(value))
    {
        stop("the length of 'value' must match the length of 'x'\n",
             call.=FALSE)
    }

    if (any(width(x) != nchar(value)))
    {
        stop("length of annotation strings must match length of sequences\n",
             call.=FALSE)
    }

    md <- metadata(x)
    md$annotationCharset <- annCharset

    emd <- mcols(x)

    if (is.null(emd))
    {
        emd <- DataFrame(annotation=value)
        rownames(emd) <- NULL
    }
    else
    {
        if (!is(emd, "DataFrame_OR_NULL"))
        {
            stop("element metadata in 'x' must be a DataFrame\n",
                 "       object or NULL\n", call.=FALSE)
        }

        emd$annotation <- value
    }

    initialize(x, metadata=md, elementMetadata=emd)
}

#' @rdname annotationSpecificKernel
#' @name annotationMetadata<-
#' @aliases
#' annotationMetadata<-
#' annotationMetadata<-,XStringSet-method
#'
#' @usage ## S4 method for signature 'XStringSet'
#' ## annotationMetadata(x, annCharset= ...) <- value
#'
#' ## S4 method for signature 'BioVector'
#' ## annotationMetadata(x, annCharset= ...) <- value
#' @param annCharset character string listing all characters used in annotation
#' sorted ascending according to the C locale, up to 32 characters are possible
#' @param value character vector with annotation strings with same length
#' as the number of sequences. Each anntation string must have the same number
#' of characters as the corresponding sequence. In addition to the characters
#' defined in the annotation character set the character "-" can be used in
#' the annotation strings for masking sequence parts.
#' @details
#' Accessor-like methods\cr\cr
#' The method annotationMetadata<- assigns annotation metadata to a sequence
#' set. In the assignment also the annotation characterset must be specified.
#' Annotation characters which are not listed in the characterset are treated
#' like invalid sequence characters. They interrupt open patterns and lead
#' to a restart of the pattern search at this position.
#' @export
#'

setReplaceMethod("annotationMetadata", "XStringSet", setAnnotationMetadata.seq)

#' @rdname annotationSpecificKernel
#' @name annotationMetadata<-
#' @aliases
#' annotationMetadata<-,BioVector-method
#' @param ... additional parameters which are passed transparently.
#' @export

setReplaceMethod("annotationMetadata", "BioVector", setAnnotationMetadata.seq)

getAnnotationMetadata.seq <- function(x)
{
    mcols(x)$annotation
}

#' @rdname annotationSpecificKernel
#' @name annotationMetadata
#' @aliases
#' annotationMetadata
#' annotationMetadata,XStringSet-method
#' @return \code{annotationMetadata}: a character vector with the annotation 
#' strings
#' @export

setMethod("annotationMetadata", "XStringSet", getAnnotationMetadata.seq)

#' @rdname annotationSpecificKernel
#' @name annotationMetadata
#' @aliases
#' annotationMetadata,BioVector-method
#' @export
#'

setMethod("annotationMetadata", "BioVector", getAnnotationMetadata.seq)

getAnnotationCharset.seq <- function(x)
{
    metadata(x)$annotationCharset
}

#' @rdname annotationSpecificKernel
#' @name annotationCharset
#' @aliases
#' annotationCharset
#' annotationCharset,XStringSet-method
#' @return \code{annotationCharset}: a character vector with the annotation 
#' @export

setMethod("annotationCharset", "XStringSet", getAnnotationCharset.seq)

#' @rdname annotationSpecificKernel
#' @name annotationMetadata
#' @aliases
#' character set 
#' annotationCharset,BioVector-method
#' @export
#'

setMethod("annotationCharset", "BioVector", getAnnotationCharset.seq)

