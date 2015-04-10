##345678901234567890123456789012345678901234567890123456789012345678901234567890
maxAnnotationCharsetLength <- 32

#' @rdname BioVector
#' @title DNAVector, RNAVector, AAVector Objects and BioVector Class
#' @name BioVector
#'
#' @description
#' Create an object containing a set of DNA-, RNA- or amino acid sequences
#'
#' @param x character vector containing a set of sequences as uppercase
#' characters or in mixed uppercase/lowercase form.
#' @details
#' The class \code{DNAVector} is used for storing DNA sequences,
#' \code{RNAVector} for RNA sequences and \code{AAVector} for
#' amino acid sequences. The class \code{BioVector} is derived from
#' the R base type \code{character} representing a vector of character
#' strings. It is an abstract class which can not be instantiated.
#' \code{BioVector} is the parent class for \code{DNAVector}, \code{RNAVector}
#' and \code{AAVector}. For the three derived classes identically named
#' functions exist which are constructors. It should be noted that the
#' constructors only wrap the sequence data into a class without copying or
#' recoding the data.
#' \cr\cr
#' The functions provided for \code{DNAVector}, \code{RNAVector} and
#' \code{AAVector} classes are only a very small subset compared to those
#' of \code{\link{XStringSet}} but are designed along their counterparts
#' from the Biostrings package. Assignment of \code{\link{metadata}}
#' and element metadata via \code{\link{mcols}} is supported for the
#' \code{DNAVector}, \code{RNAVector} and \code{AAVector} objects similar to
#' objects of \code{\link{XStringSet}} derived classes (for details on metadata
#' assignment see \code{\link{annotationMetadata}} and
#' \code{\link{positionMetadata}}).\cr\cr
#' In contrast to \code{\link{XStringSet}} the \code{BioVector} derived classes
#' also support the storage of lowercase characters. This can be relevant
#' for repeat regions which are often coded in lowercase characters. During
#' the creation of \code{\link{XStringSet}} derived classes the lowercase
#' characters are converted to uppercase automatically and the information
#' about repeat regions is lost. For \code{BioVector} derived classes the
#' user can specify during creation of a sequence kernel object whether
#' lowercase characters should be included as uppercase characters or
#' whether repeat regions should be ignored during sequence analysis.
#' In this way it is possible to perform both types of analysis on the same
#' set of sequences through defining one kernel object which accepts lowercase
#' characters and another one which ignores them.\cr\cr
#' @note
#' Sequence data can be processed by KeBABS in XStringSet and BioVector based
#' format. Within KeBABS except for treatment of lowercase characters both
#' formats are equivalent. It is recommended to use \code{\link{XStringSet}}
#' based formats whenever the support of lowercase characters is not of
#' interest because these classes provide in general much richer functionality
#' than the \code{BioVector} classes. String kernels provided in the
#' \code{kernlab} package (see \link[kernlab:stringdot]{stringdot}) do not
#' support \code{\link{XStringSet}} derived objects. The usage of these kernels
#' is possible in KeBABS with sequence data in \code{BioVector} based format.
#'
#' @return constructors \code{DNAVector, RNAVector, AAVector} return a 
#' sequence set of identical class name
#'
#' @examples
#' ## in general DNAStringSet should be prefered as described above
#' ## create DNAStringSet object for a set of sequences
#' x <- DNAStringSet(c("AACCGCGATTATCGatatatatatatatatTGGAAGCTAGGACTA",
#'                     "GACTTACCCgagagagagagagaCATGAGAGGGAAGCTAGTA"))
#' ## assign names to the sequences
#' names(x) <- c("Sample1", "Sample2")
#'
#' ## to show the different handling of lowercase characters
#' ## create DNAVector object for the same set of sequences and assign names
#' xv <- DNAVector(c("AACCGCGATTATCGatatatatatatatatTGGAAGCTAGGACTA",
#'                   "GACTTACCCgagagagagagagaCATGAGAGGGAAGCTAGTA"))
#' names(xv) <- c("Sample1", "Sample2")
#'
#' ## show DNAStringSet object - lowercase characters were translated
#' x
#' ## in the DNAVector object lowercase characters are unmodified
#' ## their handling can be defined at the level of the sequence kernel
#' xv
#'
#' ## show number of the sequences in the set and their number of characters
#' length(xv)
#' width(xv)
#' nchar(xv)
#'
## @aliases
## metadata,BioVector-method
## elementMetadata,BioVector-method
#' @seealso \code{\link{metadata}}, \code{\link{elementMetadata}},
#' \code{\link{XStringSet}}, \code{\link{DNAStringSet}},
#' \code{\link{RNAStringSet}}, \code{\link{AAStringSet}}
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

###################################################
##
## BioVector constructors / accessors
##
###################################################


#' @rdname BioVector
#' @name DNAVector
#' @aliases
#' DNAVector
#' @usage ## Constructors:
#' @export
#'

DNAVector <- function(x=character())
{
    if (!is.character(x))
        stop("sequence data must be a character vector\n")

    if (!(length(x) > 0))
        stop("sequence data is empty\n")

    if (length(names(x)) > 0)
        return(new("DNAVector", x, NAMES=names(x)))
    else
        return(new("DNAVector", x))
}

#' @rdname BioVector
#' @name RNAVector
#' @aliases
#' RNAVector
#' @export

RNAVector <- function(x=character())
{
    if (!is.character(x))
        stop("sequence data must be a character vector\n")

    if (!(length(x) > 0))
        stop("sequence data is empty\n")

    if (length(names(x)) > 0)
        return(new("RNAVector", x, NAMES=names(x)))
    else
        return(new("RNAVector", x))
}

#' @rdname BioVector
#' @name AAVector
#' @aliases
#' AAVector
#' @export

AAVector <- function(x=character())
{
    if (!is.character(x))
        stop("sequence data must be a character vector\n")

    if (!(length(x) > 0))
        stop("sequence data is empty\n")

    if (length(names(x)) > 0)
        return(new("AAVector", x, NAMES=names(x)))
    else
        return(new("AAVector", x))
}

#' @rdname BioVector
#' @aliases
#' length
#' length,BioVector-method
#' names
#' names,BioVector-method
#' names<-
#' names<-,BioVector-method
#' width
#' width,BioVector-method
#' @section Accessor-like methods:
#' In the code snippets below, \code{x} is a \code{BioVector}.
#'
#' \describe{
#'   \item{}{\code{length(x)}:
#'   the number of sequences in \code{x}.
#'   }
#'   \item{}{\code{width(x)}:
#'   vector of integer values with the number of bases/amino
#'   acids for each sequence in the set.
#'   }
#'   \item{}{\code{names(x)}:
#'   character vector of sample names.
#'   }
#' }
#' @usage ## Accessor-like methods: see below

setMethod("names", signature(x="BioVector"),
    function(x)
    {
        return(x@NAMES)
    }
)

setMethod("names<-", signature(x="BioVector"),
    function(x, value)
    {
        if (!is.null(value))
        {
            value <- as.character(value)

            if (length(value) > length(x))
                stop("number of names is larger than number of sequences\n")

            if (length(value) < length(x))
                value <- c(value, rep.int(NA, length(x) - length(value)))
        }

        x@NAMES <- value
        x
    }
)

## documented with the names function

setMethod("length", signature(x="BioVector"),
    function(x)
    {
        length(x@.Data)
    }
)

## documented with the names function

setMethod("width", signature(x="BioVector"),
    function(x)
    {
        nchar(x)
    }
)

#' @rdname BioVector
#' @aliases
#' [,BioVector-method
#' c,BioVector-method
## metadata
## elementMetadata
## mcols
#' @param i numeric vector with indicies or character with element names
#' @section Subsetting and concatination:
#' In the code snippets below, \code{x} is a \code{BioVector}.
#'
#' \describe{
#'   \item{}{\code{x[i]}:
#'   return a \code{BioVector} object that only contains the samples selected
#'   with the subsetting parameter \code{i}. This parameter can be a numeric
#'   vector with indices or a character vector which is matched against the
#'   names of \code{x}. Element related metadata is subsetted accordingly if
#'   available.
#'   }
#'   \item{}{\code{c(x, ...)}:
#'   return a sequence set that is a concatination of the given sequence sets.
#'   }
#' }
#'

setMethod("[", signature(x="BioVector", i="index", j="missing"),
    function(x, i)
    {
        if (is.character(i))
        {
            if (length(x@NAMES) < 1 || any(is.na(x@NAMES)))
                stop("missing names for subsetting\n")

            i1 <- which(x@NAMES %in% i)

            if (length(i) != length(i1))
                stop("invalid names specified\n")

            i <- i1
        }
        else
        {
            if (min(abs(i)) < 1 || max(abs(i)) > length(x))
                stop("row subset must be between 1 and number of sequences\n")

            ## convert negative subset
            if (all(i < 0))
                i <- (1:length(x))[i]
            else
            {
                if (min(i) < 1)
                {
                    stop("subset indices must be all positive or all",
                         " negative\n")
                }
            }
        }

        if (length(x@NAMES) > 0)
            newNames <- x@NAMES[i]
        else
            newNames <- NULL

        initialize(x, .Data=x@.Data[i], NAMES=newNames, metadata=x@metadata,
                   elementMetadata=x@elementMetadata[i,,drop=FALSE])
    }
)

setMethod("c", signature(x="BioVector"),
    function(x, ..., recursive=FALSE)
    {
        initialize(x, .Data=c(x@.Data, ...), NAMES=c(x@NAMES, names(...)),
                   metadata=x@metadata, elementMetadata=
                   rbind(x@elementMetadata, elementMetadata(...)))
    }
)

## @rdname BioVector
# @aliases
## metadata,BioVector-method
## @section Assignment of metadata and element metadata:
## For the assignment of metadata and element metadata please use the
## functions \code{\link{annotationMetadata}} and
## \code{\link{positionMetadata}}) which take care of many aspects not handled
## in the low level access methods described here.
##
## \describe{
##   \item{}{\code{metadata}:
##   return a \code{BioVector} object that only contains the samples selected
##   with the subsetting parameter \code{i}. This parameter can be a numeric
##   vector with indices or a character vector which is matched against the
##   names of \code{x}. Element related metadata is subsetted accordingly if
##   available.
##   }
##   \item{}{\code{elementMetadata}:
##   deprecated.
##   }
##   \item{}{\code{mcols}:
##   return a DataFrame object that contains the element metadata the given
##   sequence sets as columns. KeBABS uses the column names "offset" for
##   position metadata and "annotation" for annotation sequences. The use of
##   elementMetadata instead of mcols is deprecated.
##   }
## }
##

## @rdname BioVector
## @aliases metadata

setMethod("metadata", signature(x="BioVector"),
    function(x)
    {
        x@metadata
    }
)

setReplaceMethod("metadata", signature(x="BioVector"),
    function(x, value)
    {
        if (!is.list(value))
            stop("replacement 'metadata' value must be a list\n")

        if (!length(value))
            names(value) <- NULL

        if (length(names(value) > 0))
        {
            if ("annotationCharset" %in% names(value))
            {
                annCharSet <- value[["annotationCharset"]]

                if (!is.character(annCharSet) || nchar(annCharSet) >
                    maxAnnotationCharsetLength)
                {
                    stop("annotation charset must be a character string\n",
                         paste("       with up to",
                               maxAnnotationCharsetLength, "characters\n"))
                }

                if (anyDuplicated(strsplit(annCharSet, split="")[[1]]) > 0)
                {
                    stop("annotation characterset must not contain\n",
                         "       duplicate characters\n")
                }

                if (gregexpr("\\.", annCharSet)[[1]][1] != -1)
                    stop("annotation characterset must not contain '.'\n")
            }
        }

        x@metadata <- value

        x
    }
)

setMethod("elementMetadata", signature(x="BioVector"),
    function(x)
    {
        x@elementMetadata
    }
)

setReplaceMethod("elementMetadata", signature(x="BioVector"),
    function(x, ..., value)
    {
        if (!is(value, "DataTableORNULL"))
            stop("replacement value must be a DataTable object or NULL\n")

        if (!is.null(value) && length(x) != nrow(value))
        {
            stop("the number of rows in 'value' (if non-NULL) must match \n",
                 "       the length of 'x'\n")
        }

        if (!is.null(value))
            rownames(value) <- NULL

        if (!is.null(value$annotation) && !is.character(value$annotation))
        {
            ## try to convert to character vector
            value$annotation <- tryCatch(as.character(value$annotation),
                                         warning=function(w) {stop(w)},
                                         error=function(e) {stop(e)})
        }

        if (!is.null(value$annotation) &&
            (any(nchar(value$annotation) != width(x))))
            stop("annotation length is not matching sequence length\n")


        x@elementMetadata <- value

        x
    }
)

## @rdname BioVector
## @name mcols
## @aliases
## mcols
## mcols,BioVector-method
##

setMethod("mcols", signature(x="BioVector"),
    function(x)
    {
        x@elementMetadata
    }
)

## @rdname BioVector
## @name mcols<-
## @aliases
## mcols<-
## mcols<-,BioVector-method
##

setReplaceMethod("mcols", signature(x="BioVector"),
    function(x, ..., value)
    {
        `elementMetadata<-`(x, ..., value = value)
    }
)

###################################################
##
## CrossValidationResult accessors
##
###################################################

#' @rdname CrossValidationResultAccessors
#' @title CrossValidationResult Accessors
#' @name CrossValidationResultAccessors
#' @aliases
#' folds
#' folds,CrossValidationResult-method
#' performance,CrossValidationResult-method
#' @param object a cross validation result object (can be extracted from
#' KeBABS model with accessor \code{\link{cvResult}})
#' @section Accessor-like methods:
#'
#' \describe{
#'   \item{}{\code{folds}:
#'   return the CV folds.
#'   }
#'   \item{}{\code{performance}:
#'   return the collected performance parameters.
#'   }
#' }
#' @return \code{folds}: returns the folds used in CV\cr
#' \code{performance}: returns a list with the performance values
#' @examples
#' ## create kernel object for normalized spectrum kernel
#' specK5 <- spectrumKernel(k=5)
#' \dontrun{
#' ## load data
#' data(TFBS)
#'
#' ## perform training - feature weights are computed by default
#' model <- kbsvm(enhancerFB, yFB, specK5, pkg="LiblineaR",
#'                svm="C-svc", cross=10, cost=15, perfParameters="ALL")
#'
#' ## show model selection result
#' cvResult(model)
#'
#' ## extract fold AUC
#' performance(cvResult(model))$foldAUC
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setMethod("folds", "CrossValidationResult",
function(object) object@folds)

performance.CrossValidationResult <- function(object)
{
    result <- list()

    result$cvError <- object@cvError
    result$foldErrors <- object@foldErrors
    result$noSV <- object@noSV

    if ("ACC" %in% object@perfParameters)
    {
        result$ACC <- object@ACC
        result$foldACC <- object@foldACC
    }
    if ("BACC" %in% object@perfParameters)
    {
        result$BACC <- object@BACC
        result$foldBACC <- object@foldBACC
    }
    
    if ("MCC" %in% object@perfParameters)
    {
        result$MCC <- object@MCC
        result$foldMCC <- object@foldMCC
    }
    
    if ("AUC" %in% object@perfParameters)
    {
        result$AUC <- object@AUC
        result$foldAUC <- object@foldAUC
    }
    
    result
}

setMethod("performance", "CrossValidationResult",
performance.CrossValidationResult)

###################################################
##
## ModelSelectionResult accessors
##
###################################################

#' @rdname ModelSelectionResultAccessors
#' @title ModelSelectionResult Accessors
#' @name ModelSelectionResultAccessors
#' @aliases
#' gridRows
#' gridRows,ModelSelectionResult-method
#' gridColumns
#' gridColumns,ModelSelectionResult-method
#' gridErrors
#' gridErrors,ModelSelectionResult-method
#' performance
#' performance,ModelSelectionResult-method
#' selGridRow
#' selGridRow,ModelSelectionResult-method
#' selGridCol
#' selGridCol,ModelSelectionResult-method
#' fullModel
#' fullModel,ModelSelectionResult-method
#' @param object a model selection result object (can be extracted from
#' KeBABS model with accessor \code{\link{modelSelResult}})
#' @section Accessor-like methods:
#'
#' \describe{
#'   \item{}{\code{gridRows}:
#'   return the grid rows containing the kernels.
#'   }
#'   \item{}{\code{gridColumns}:
#'   return the grid columns.
#'   }
#'   \item{}{\code{gridErrors}:
#'   return the grid CV errors.
#'   }
#'   \item{}{\code{performance}:
#'   return the collected performance parameters.
#'   }
#'   \item{}{\code{selGridRow}:
#'   return the selected grid row.
#'   }
#'   \item{}{\code{selGridCol}:
#'   return the selected grid column.
#'   }
#'   \item{}{\code{fullModel}:
#'   return the full model.
#'   }
#' }
#' @return \code{gridRows}: returns a list of kernel objects\cr
#' \code{gridColumns}: returns a \code{DataFrame} object with grid column
#' parameters\cr
#' \code{gridErrors}: returns a matrix with grid errors\cr
#' \code{performance}: returns a list of matrices with performance values
#' \code{selGridRow}: returns the selected kernel
#' \code{selGridCol}: returns the selected SVM and/or hyperparameter(s)
#' \code{fullModel}: returns a kebabs model of class
#' \code{\linkS4class{KBModel}}
#' @examples
#' ## create kernel object for normalized spectrum kernel
#' specK5 <- spectrumKernel(k=5)
#' \dontrun{
#' ## load data
#' data(TFBS)
#'
#' ## perform training - feature weights are computed by default
#' model <- kbsvm(enhancerFB, yFB, specK5, pkg="LiblineaR", 
#'                svm="C-svc", cost=c(1,15,50,100), cross=10, 
#'                perfParameters="ALL", showProgress=TRUE)
#'
#' ## show model selection result
#' mres <- modelSelResult(model)
#' mres
#'
#' ## extract grid errors
#' gridErrors(mres)
#'
#' ## extract other performance parameters
#' performance(mres)
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setMethod("gridRows", "ModelSelectionResult",
          function(object) object@gridRows)
setMethod("gridColumns", "ModelSelectionResult",
          function(object) object@gridCols)
setMethod("gridErrors", "ModelSelectionResult",
          function(object) object@gridErrors)

performance.ModelSelectionResult <- function(object)
{
    result <- list()

    if ("ACC" %in% object@perfParameters)
        result$ACC <- object@gridACC

    if ("BACC" %in% object@perfParameters)
        result$BACC <- object@gridBACC

    if ("MCC" %in% object@perfParameters)
        result$MCC <- object@gridMCC

    if ("AUC" %in% object@perfParameters)
        result$AUC <- object@gridAUC

    result
}

setMethod("performance", "ModelSelectionResult",
          performance.ModelSelectionResult)

setMethod("selGridRow", "ModelSelectionResult",
          function(object) object@selGridRow)
setMethod("selGridCol", "ModelSelectionResult",
          function(object) object@selGridCol)
setMethod("fullModel", "ModelSelectionResult",
          function(object) object@fullModel)

###################################################
##
## ROCData accessors
##
###################################################

#' @rdname ROCDataAccessors
#' @name ROCDataAccessors
#' @title ROCData Accessors
#' @aliases
#' auc
#' auc,ROCData-method
#' auc<-
#' auc<-,ROCData-method
#' tpr
#' tpr,ROCData-method
#' tpr<-
#' tpr<-,ROCData-method
#' fpr
#' fpr,ROCData-method
#' fpr<-
#' fpr<-,ROCData-method
#' @param object an object of class \code{\linkS4class{ROCData}}
#' @section Accessor-like methods:
#'
#' \describe{
#'   \item{}{\code{auc}:
#'   returns the area under the ROC curve.
#'   }
#'   \item{}{\code{tpr}:
#'   returns the true positive rate values as numeric vector.
#'   }
#'   \item{}{\code{fpr}:
#'   returns the false positive rate values as numeric vector.
#'   }
#' }
#' @return \code{auc}: returns a numeric value\cr
#' \code{tpr}: returns a numeric vector\cr
#' \code{fpr}: returns a numeric vector\cr
#' @examples
#' ## create kernel object for normalized spectrum kernel
#' specK5 <- spectrumKernel(k=5)
#' \dontrun{
#' ## load data
#' data(TFBS)
#'
#' ## select 70% of the samples for training and the rest for test
#' train <- sample(1:length(enhancerFB), length(enhancerFB) * 0.7)
#' test <- c(1:length(enhancerFB))[-train]
#'
#' ## perform training - feature weights are computed by default
#' model <- kbsvm(enhancerFB[train], yFB[train], specK5, pkg="LiblineaR",
#'                svm="C-svc", cost=15)
#' preddec <- predict(model, enhancerFB[test], predictionType="decision")
#' rocdata <- computeROCandAUC(preddec, yFB[test], allLabels=unique(yFB))
#'
#' ## accessor for auc
#' auc(rocdata)
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

## accessor auc
setMethod("auc", "ROCData", function(object) object@AUC)

setReplaceMethod("auc", "ROCData",
    function(x, value)
    {
        x@auc <- value
        x
    }
)

## accessor tpr
setMethod("tpr", "ROCData", function(object) object@TPR)

setReplaceMethod("tpr", "ROCData",
    function(x, value)
    {
        x@tpr <- value
        x
    }
)

## accessor fpr
setMethod("fpr", "ROCData", function(object) object@FPR)

setReplaceMethod("fpr", "ROCData",
    function(x, value)
    {
        x@fpr <- value
        x
    }
)


###################################################
##
## KBModel accessors
##
###################################################

#' @rdname KBModelAccessors
#' @name KBModelAccessors
#' @title KBModel Accessors
#' @aliases
#' modelOffset
#' modelOffset,KBModel-method
#' modelOffset<-
#' modelOffset<-,KBModel-method
#' featureWeights
#' featureWeights,KBModel-method
#' featureWeights<-
#' featureWeights<-,KBModel-method
#' SVindex
#' SVindex,KBModel-method
#' SVindex<-
#' SVindex<-,KBModel-method
#' cvResult
#' cvResult,KBModel-method
#' cvResult<-
#' cvResult<-,KBModel-method
#' modelSelResult
#' modelSelResult,KBModel-method
#' modelSelResult<-
#' modelSelResult<-,KBModel-method
#' svmModel
#' svmModel,KBModel-method
#' svmModel<-
#' svmModel<-,KBModel-method
#' probabilityModel
#' probabilityModel,KBModel-method
#' probabilityModel<-
#' probabilityModel<-,KBModel-method
#' @param object a KeBABS model
#' @section Accessor-like methods:
#'
#' \describe{
#'   \item{}{\code{modelOffset}:
#'   returns the model offset.
#'   }
#'   \item{}{\code{featureWeights}:
#'   returns the feature weights.
#'   }
#'   \item{}{\code{SVindex}:
#'   returns the support vector indices for the training samples.
#'   }
#'   \item{}{\code{cvResult}:
#'   returns result of cross validation as object of class
#'   \code{\linkS4class{CrossValidationResult}}.
#'   }
#'   \item{}{\code{modelSelResult}:
#'   returns result of model selection as object of class
#'   \code{\linkS4class{ModelSelectionResult}}.
#'   }
#'   \item{}{\code{svmModel}:
#'   returns the native svm model stored within KeBABS model.
#'   }
#'   \item{}{\code{probabilityModel}:
#'   returns the probability model stored within KeBABS model.
#'   }
#' }
#' @examples
#' ## create kernel object for normalized spectrum kernel
#' specK5 <- spectrumKernel(k=5)
#' \dontrun{
#' ## load data
#' data(TFBS)
#'
#' ## perform training - feature weights are computed by default
#' model <- kbsvm(enhancerFB, yFB, specK5, pkg="LiblineaR", 
#'                svm="C-svc", cost=15, cross=10, showProgress=TRUE)
#'                showProgress=TRUE)
#'
#' ## show result of validation
#' cvResult(model)
#' ## show feature weights
#' featureWeights(model)[1:5]
#' ## show model offset
#' modelOffset(model)
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.


## accessor modelOffset for b
setMethod("modelOffset", "KBModel", function(object) object@b)

setReplaceMethod("modelOffset", "KBModel",
    function(x, value)
    {
        x@b <- value
        x
    }
)

setMethod("featureWeights", "KBModel", function(object) object@featureWeights)

setReplaceMethod("featureWeights", "KBModel",
    function(x, value)
    {
        x@featureWeights <- value
        x
    }
)

setMethod("SVindex", "KBModel", function(object) object@svIndex)

setReplaceMethod("SVindex", "KBModel",
    function(x, value)
    {
        x@svIndex <- value
        x
    }
)

setMethod("cvResult", "KBModel", function(object) object@cvResult)

setReplaceMethod("cvResult", "KBModel",
    function(x, value)
    {
        if (!is(value, "CrossValidationResult"))
        {
            stop("new value must be of class 'CrossValidationResult'\n",
                class. = FALSE)
        }
        x@cvResult <- value
        x
    }
)

setMethod("modelSelResult", "KBModel", function(object) object@modelSelResult)

setReplaceMethod("modelSelResult", "KBModel",
    function(x, value)
    {
        if (!is(value, "ModelSelectionResult"))
        {
            stop("new value must be of class 'ModelSelectionResult'\n",
                 class. = FALSE)
        }
        x@modelSelResult <- value
        x
    }
)

setMethod("svmModel", "KBModel", function(object) object@svmModel)

setReplaceMethod("svmModel", "KBModel",
    function(x, value)
    {
        x@svmModel <- value
        x
    }
)

setMethod("probabilityModel", "KBModel",
    function(object)
    {
        list(probA=object@probA, probB=object@probB, sigma=object@sigma)
    }
)

setReplaceMethod("probabilityModel", "KBModel",
    function(x, value)
    {
        if (!(is.list(value) &&
              (sort(names(value)) == c("probA", "probB", "sigma"))))
            stop("invalid probability model parameters")

        x@probA <- value$probA
        x@probB <- value$probB
        x@sigma <- value$sigma
        x
    }
)


###################################################
##
## KernelMatrix accessors
##
###################################################

#' @rdname KernelMatrixAccessors
#' @title KernelMatrix Accessors
#' @name KernelMatrixAccessors
#' @aliases
#' [,KernelMatrix,index,index,ANY-method
#' [,KernelMatrix,index,missing,ANY-method
#' [,KernelMatrix,missing,index,ANY-method
#' @param i numeric vector with indicies or character with element names
#' @section Accessor-like methods:
#'
#' \describe{
#'   \item{}{\code{x[i,]}:
#'   return a \code{KernelMatrix} object that only contains the rows selected
#'   with the subsetting parameter i. This parameter can be a numeric vector
#'   with indices or a character vector which is matched against the names of
#'   \code{x}.
#'   }
#'   \item{}{\code{x[,j]}:
#'   return a \code{KernelMatrix} object that only contains the columns selected
#'   with the subsetting parameter j. This parameter can be a numeric vector
#'   with indices or a character vector which is matched against the names of
#'   \code{x}.
#'   }
#'   \item{}{\code{x[i,j]}:
#'   return a \code{KernelMatrix} object that only contains the rows selected
#'   with the subsetting parameter i and columns selected by j. Both parameters
#'   can be a numeric vector with indices or a character vector which is matched
#'   against the names of \code{x}.
#'   }
#' }
#' @return see above
#' @examples
#' ## create kernel object for normalized spectrum kernel
#' specK5 <- spectrumKernel(k=5)
#' \dontrun{
#' ## load data
#' data(TFBS)
#'
#' km <- specK5(enhancerFB)
#' km1to5 <- km[1:5,1:5]
#' km1to5
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.


setMethod("[", signature(x="KernelMatrix", i="index", j="missing"),
    function(x, i)
    {
        if (is.character(i))
        {
            if (length(rownames(x)) < 1 || any(is.na(rownames(x))))
                stop("missing rownames for subsetting\n")
            else
                i1 <- which(rownames(x) %in% i)

            if (length(i) != length(i1))
                stop("invalid rownames specified\n")

            i <- i1
        }
        else
        {
            if (min(abs(i)) < 1 || max(abs(i)) > nrow(x))
                stop("row subset must be between 1 and number of rows\n")

            ## convert negative subset
            if (all(i < 0))
                i <- (1:nrow(x))[i]
            else
            {
                if (min(i) < 1)
                    stop("subset indices must be all positive or",
                         " all negative\n")
            }
        }

        initialize(x, .Data=x@.Data[i,,drop=FALSE])
    }
)

setMethod("[", signature(x="KernelMatrix", i="missing", j="index"),
    function(x, j)
    {
        if (is.character(j))
        {
            if (length(colnames(x)) < 1 || any(is.na(colnames(x))))
                stop("missing colnames for subsetting\n")
            else
                j1 <- which(colnames(x) %in% j)

            if (length(j) != length(j1))
                stop("invalid colnames specified\n")

            j <- j1
        }
        else
        {
            if (min(abs(j)) < 1 || max(abs(j)) > ncol(x))
                stop("column subset must be between 1 and",
                     " number of columns\n")

            ## convert negative subset
            if (all(j < 0))
                j <- (1:ncol(x))[j]
            else
            {
                if (min(j) < 1)
                    stop("subset indices must be all positive or",
                         " all negative\n")
            }
        }

        initialize(x, .Data=x@.Data[, j, drop=FALSE])
    }
)

setMethod("[", signature(x="KernelMatrix", i="index", j="index"),
    function(x, i, j)
    {
        if (is.character(i))
        {
            if (length(rownames(x)) < 1 || any(is.na(rownames(x))))
                stop("missing rownames for subsetting\n")
            else
                i1 <- which(rownames(x) %in% i)

            if (length(i) != length(i1))
                stop("invalid rownames specified\n")

            i <- i1
        }
        else
        {
            if (min(abs(i)) < 1 || max(abs(i)) > nrow(x))
                stop("row subset must be between 1 and number of rows\n")

            ## convert negative subset
            if (all(i < 0))
                i <- (1:nrow(x))[i]
            else
            {
                if (min(i) < 1)
                    stop("subset indices must be all positive or",
                         " all negative\n")
            }
        }

        if (is.character(j))
        {
            if (length(colnames(x)) < 1 || any(is.na(colnames(x))))
                stop("missing colnames for subsetting\n")
            else
                j1 <- which(colnames(x) %in% j)

            if (length(j) != length(j1))
                stop("invalid colnames specified\n")

            j <- j1
        }
        else
        {
            if (min(j) < 1 || max(j) > ncol(x))
                stop("column subset must be between 1 and",
                     " number of columns\n")

            ## convert negative subset
            if (all(j < 0))
                j <- (1:ncol(x))[j]
            else
            {
                if (min(j) < 1)
                    stop("subset indices must be all positive or",
                         " all negative\n")
            }
        }

        initialize(x, .Data=x@.Data[i, j, drop=FALSE])
    }
)

###################################################
##
## Explicit representation subsetting
##
###################################################

#' @rdname ExplicitRepresentationAccessors
#' @title ExplicitRepresentation Accessors
#' @name ExplicitRepresentationAccessors
#' @aliases
#' [,ExplicitRepresentation,index,index,ANY-method
#' [,ExplicitRepresentationDense,index,index,ANY-method
#' [,ExplicitRepresentationDense,index,missing,ANY-method
#' [,ExplicitRepresentationDense,missing,index,ANY-method
#' [,ExplicitRepresentationSparse,index,index,ANY-method
#' [,ExplicitRepresentationSparse,index,index,logical-method
#' [,ExplicitRepresentationSparse,index,index,missing-method
#' [,ExplicitRepresentationSparse,index,missing,ANY-method
#' [,ExplicitRepresentationSparse,index,missing,logical-method
#' [,ExplicitRepresentationSparse,index,missing,missing-method
#' [,ExplicitRepresentationSparse,missing,index,ANY-method
#' [,ExplicitRepresentationSparse,missing,index,logical-method
#' [,ExplicitRepresentationSparse,missing,index,missing-method
#'
#' @usage ## S4 methods for signature 'ExplicitRepresentation'
#' ## x[i,j]
#'
#' ## further methods see below
#' @param x an explicit representation in dense or sparse format
#' @param i integer vector or character vector with a subset of the sample
#' indices or names
#' @param j integer vector or character vector with a subset of the feature
#' indices or names
#'
#' @section Accessor-like methods:
#'
#' \describe{
#'   \item{}{\code{x[i,]}:
#'   return a \code{KernelMatrix} object that only contains the rows selected
#'   with the subsetting parameter \code{i}. This parameter can be a numeric
#'   vector with indices or a character vector which is matched against the
#'   names of \code{x}.
#'   }
#'   \item{}{\code{x[,j]}:
#'   return a \code{KernelMatrix} object that only contains the columns
#'   selected with the subsetting parameter \code{j}. This parameter can be a
#'   numeric vector with indices or a character vector which is matched against
#'   the names of \code{x}.
#'   }
#'   \item{}{\code{x[i,j]}:
#'   return a \code{KernelMatrix} object that only contains the rows selected
#'   with the subsetting parameter \code{i} and columns selected by \code{j}.
#'   Both parameters can be a numeric vector with indices or a character vector
#'   which is matched against the names of \code{x}.
#'   }
#' }
#' @return see details above 
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setMethod("[", signature(x="ExplicitRepresentationDense", i="index",
                         j="missing"),
    function(x, i)
    {
        if (is.character(i))
        {
            if (length(rownames(x)) < 1 || any(is.na(rownames(x))))
                stop("missing rownames for subsetting\n")
            else
                i1 <- which(rownames(x) %in% i)

            if (length(i) != length(i1))
                stop("invalid rownames specified\n")

            i <- i1
        }
        else
        {
            if (min(abs(i)) < 1 || max(abs(i)) > nrow(x))
                stop("row subset must be between 1 and number of rows\n")

            ## convert negative subset
            if (all(i < 0))
                i <- (1:nrow(x))[i]
            else
            {
                if (min(i) < 1)
                    stop("subset indices must be all positive or",
                         " all negative\n")
            }
        }

        initialize(x, .Data=x@.Data[i, , drop=FALSE],
                   usedKernel=x@usedKernel)
    }
)

setMethod("[", signature(x="ExplicitRepresentationDense", i="missing",
                         j="index"),
    function(x, j)
    {
        if (is.character(j))
        {
            if (length(colnames(x)) < 1 || any(is.na(colnames(x))))
                stop("missing colnames for subsetting\n")
            else
                j1 <- which(colnames(x) %in% j)

            if (length(j) != length(j1))
                stop("invalid colnames specified\n")

            j <- j1
        }
        else
        {
            if (min(abs(j)) < 1 || max(abs(j)) > ncol(x))
                stop("column subset must be between 1 and",
                     " number of columns\n")

            ## convert negative subset
            if (all(j < 0))
                j <- (1:ncol(x))[j]
            else
            {
                if (min(j) < 1)
                    stop("subset indices must be all positive or",
                         " all negative\n")
            }
        }

        ## maintain feature order
        j <- sort(j)

        initialize(x, .Data=x@.Data[, j, drop=FALSE],
                   usedKernel=x@usedKernel)
    }
)

setMethod("[", signature(x="ExplicitRepresentationDense", i="index",
                         j="index"),
    function(x, i, j)
    {
        if (is.character(i))
        {
            if (length(rownames(x)) < 1 || any(is.na(rownames(x))))
                stop("missing rownames for subsetting\n")
            else
                i1 <- which(rownames(x) %in% i)

            if (length(i) != length(i1))
                stop("invalid rownames specified\n")

            i <- i1
        }
        else
        {
            if (min(abs(i)) < 1 || max(abs(i)) > nrow(x))
                stop("row subset must be between 1 and number of rows\n")

            ## convert negative subset
            if (all(i < 0))
                i <- (1:nrow(x))[i]
            else
            {
                if (min(i) < 1)
                stop("subset indices must be all positive or all negative\n")
            }
        }

        if (is.character(j))
        {
            if (length(colnames(x)) < 1 || any(is.na(colnames(x))))
                stop("missing colnames for subsetting\n")
            else
                j1 <- which(colnames(x) %in% j)

            if (length(j) != length(j1))
                stop("invalid colnames specified\n")

            j <- j1
        }
        else
        {
            if (min(abs(j)) < 1 || max(abs(j)) > ncol(x))
                stop("column subset must be between 1 and number",
                     " of columns\n")

            ## convert negative subset
            if (all(j < 0))
                j <- (1:ncol(x))[j]
            else
            {
                if (min(j) < 1)
                    stop("subset indices must be all positive or",
                         " all negative\n")
            }
        }

        ## maintain feature order
        j <- sort(j)

        initialize(x, .Data=x@.Data[i, j, drop=FALSE],
                   usedKernel=x@usedKernel)
    }
)

subsetERSSample <- function(x, i, drop=FALSE)
{
    if (is.character(i))
    {
        if (length(i) == 0)
        {
            return(initialize(x, Dim=c(0L, x@Dim[2]), factors=list(),
                              p=0L, j=integer(0), x=numeric(0),
                              Dimnames=list(NULL, x@Dimnames[[2]]),
                              quadratic=x@quadratic, usedKernel=x@usedKernel))
        }

        if (length(rownames(x)) < 1 || any(is.na(rownames(x))))
            stop("missing rownames for subsetting\n")
        else
            i1 <- which(rownames(x) %in% i)

        if (length(i) != length(i1))
            stop("invalid rownames specified\n")

        i <- i1
    }
    else
    {
        if (!is.numeric(i))
            stop("wrong type of column indices for subsetting\n")

        if (length(i) == 0)
        {
            return(initialize(x, Dim=c(0L, x@Dim[2]), factors=list(),
                        p=0L, j=integer(0), x=numeric(0),
                        Dimnames=list(NULL, x@Dimnames[[2]]),
                        quadratic=x@quadratic, usedKernel=x@usedKernel))
        }

        if (min(abs(i)) < 1 || max(abs(i)) > nrow(x))
            stop("row subset must be between 1 and number of rows\n")
    }

    ## convert negative subset
    if (all(i < 0))
        i <- (1:nrow(x))[i]
    else
    {
        if (min(i) < 1)
            stop("subset indices must be all positive or all negative\n")
    }

    if (length(x@Dimnames[[1]]) > 0)
        dimnames1 <- x@Dimnames[[1]][i]
    else
        dimnames1 <- NULL

    initialize(x, Dim=c(length(i), x@Dim[2]), factors=list(),
        p=as.integer(c(0, cumsum(x@p[i+1] - x@p[i]))),
        j=c(unlist(sapply(i,
                function(ind){x@j[if (x@p[ind] == x@p[ind+1])
                    {numeric(0)} else {(x@p[ind]+1):x@p[ind+1]}]}))),
        x=c(unlist(sapply(i,
                function(ind){x@x[if (x@p[ind] == x@p[ind+1])
                    {numeric(0)} else {(x@p[ind]+1):x@p[ind+1]}]}))),
               Dimnames=list(dimnames1, x@Dimnames[[2]]),
        quadratic=x@quadratic, usedKernel=x@usedKernel)
}

setMethod("[", signature(x="ExplicitRepresentationSparse", i="index",
                         j="missing", drop="missing"), subsetERSSample)
setMethod("[", signature(x="ExplicitRepresentationSparse", i="index",
                         j="missing", drop="logical"), subsetERSSample)
setMethod("[", signature(x="ExplicitRepresentationSparse", i="index",
                         j="missing", drop="ANY"), subsetERSSample)

subsetERSFeature <- function(x, j, drop=FALSE)
{
    if (is.character(j))
    {
        if (length(j) == 0)
        {
            return(initialize(x, Dim=c(x@Dim[1], 0L), factors=list(),
                        p=rep(0L, x@Dim[1]+1), j=integer(0), x=numeric(0),
                        Dimnames=list(x@Dimnames[[1]], NULL),
                        quadratic=x@quadratic, usedKernel=x@usedKernel))
        }

        if (length(colnames(x)) < 1 || any(is.na(colnames(x))))
            stop("missing colnames for subsetting\n")
        else
            j1 <- which(colnames(x) %in% j)

        if (length(j) != length(j1))
            stop("invalid colnames specified\n")

        j <- j1
    }
    else
    {
        if (!is.numeric(j))
            stop("wrong type of column indices for subsetting\n")

        if (length(j) == 0)
        {
            return(initialize(x, Dim=c(x@Dim[1], 0L), factors=list(),
                        p=rep(0L, x@Dim[1]+1), j=integer(0), x=numeric(0),
                        Dimnames=list(x@Dimnames[[1]], NULL),
                        quadratic=x@quadratic, usedKernel=x@usedKernel))
        }

        if (min(abs(j)) < 1 || max(abs(j)) > ncol(x))
            stop("subsetting of explicit representaion with subscript\n",
                 "        out of bounds\n")
    }

    ## convert negative subset
    if (all(j < 0))
        j <- (1:ncol(x))[j]
    else
    {
        if (min(j) < 1)
            stop("subset indices for explicit representation must be all\n",
                 "        positive or all negative\n")

        ## maintain feature order
        j <- sort(j)
    }

    if (length(x@Dimnames[[2]]) > 0)
        dimnames2 <- x@Dimnames[[2]][j]
    else
        dimnames2 <- NULL

    if (nrow(x) == 0)
    {
        return(initialize(x, Dim=c(x@Dim[1], length(j)), factors=list(),
                          p=0L, j=integer(0), x=numeric(0),
                          Dimnames=list(NULL, dimnames2),
                          quadratic=x@quadratic, usedKernel=x@usedKernel))
    }

    colIndMap <- rep(-1, j[length(j)])

    for (k in 1:length(j))
        colIndMap[j[k]] <- k - 1

    p1 <- sapply(1:nrow(x),
                 function(ind)
                 {
                    if (x@p[ind] == x@p[ind+1])
                        return(c())
                    else
                        return(x@j[(x@p[ind]+1):x@p[ind+1]])
                 })

    if (is.matrix(p1))
    {
        p1 <- c(0,
                cumsum(apply(p1,2,
                             function(x){length(which(x %in% (j-1)))})))
    }
    else
    {
        p1 <- c(0,
                cumsum(sapply(p1,
                              function(x){length(which(x %in% (j-1)))})))
    }

    indices <- which(x@j %in% (j-1))

    initialize(x, Dim=c(x@Dim[1], length(j)), factors=list(),
               p=as.integer(p1),
               j=as.integer(colIndMap[x@j[indices]+1]),
               x=x@x[indices],
               Dimnames=list(x@Dimnames[[1]], dimnames2),
               quadratic=x@quadratic, usedKernel=x@usedKernel)
}

setMethod("[", signature(x="ExplicitRepresentationSparse", i="missing",
                         j="index", drop="missing"), subsetERSFeature)
setMethod("[", signature(x="ExplicitRepresentationSparse", i="missing",
                         j="index", drop="logical"), subsetERSFeature)
setMethod("[", signature(x="ExplicitRepresentationSparse", i="missing",
                         j="index", drop="ANY"), subsetERSFeature)

subsetERSSampleFeature <- function(x, i, j, drop=FALSE)
{
    if ((is.character(i) || is.numeric(i)) &&
        (is.character(j) || is.numeric(j)))
    {
        if (length(i) == 0 || length(j) == 0)
        {
            if (length(i) == 0 && length(j) == 0)
            {
                return(initialize(x, Dim=c(0L, 0L), factors=list(),
                       p=0L, j=integer(0), x=numeric(0),
                       Dimnames=list(NULL, NULL),
                       quadratic=x@quadratic,
                       usedKernel=x@usedKernel))

            }
        }
        else if (length(i) == 0)
        {
            return(initialize(x, Dim=c(0L, x@Dim[2]), factors=list(),
                              p=0L, j=integer(0), x=numeric(0),
                              Dimnames=list(NULL, x@Dimnames[[2]]),
                              quadratic=x@quadratic, usedKernel=x@usedKernel))
        }
        else if (length(j) == 0)
        {
            return(initialize(x, Dim=c(x@Dim[1], 0L), factors=list(),
                              p=rep(0L, x@Dim[1]+1), j=integer(0), x=numeric(0),
                              Dimnames=list(x@Dimnames[[1]], NULL),
                              quadratic=x@quadratic, usedKernel=x@usedKernel))
        }
    }
    else
        stop("wrong type of row or column indices for subsetting\n")

    if (is.character(i))
    {
        if (length(rownames(x)) < 1 || any(is.na(rownames(x))))
            stop("missing rownames for subsetting\n")
        else
            i1 <- which(rownames(x) %in% i)

        if (length(i) != length(i1))
            stop("invalid rownames specified\n")

        i <- i1
    }
    else
    {
        if (min(abs(i)) < 1 || max(abs(i)) > nrow(x))
            stop("row subset must be between 1 and number of rows\n")
    }

    ## convert negative row subset
    if (all(i < 0))
        i <- (1:nrow(x))[i]
    else
    {
        if (min(i) < 1)
            stop("subset indices must be all positive or all negative\n")
    }

    if (is.character(j))
    {
        if (length(colnames(x)) < 1 || any(is.na(colnames(x))))
            stop("missing colnames for subsetting\n")
        else
            j1 <- which(colnames(x) %in% j)

        if (length(j) != length(j1))
            stop("invalid colnames specified\n")

        j <- j1
    }
    else
    {
        if (min(abs(j)) < 1 || max(abs(j)) > ncol(x))
            stop("column subset must be between 1 and number of columns\n")
    }

    ## convert negative col subset
    if (all(j < 0))
        j <- (1:ncol(x))[j]
    else
    {
        if (min(j) < 1)
            stop("subset indices must be all positive or all negative\n")

        ## maintain feature order
        j <- sort(j)
    }

    if (length(x@Dimnames[[1]]) > 0)
        dimnames1 <- x@Dimnames[[1]][i]
    else
        dimnames1 <- NULL

    if (length(x@Dimnames[[2]]) > 0)
        dimnames2 <- x@Dimnames[[2]][j]
    else
        dimnames2 <- NULL

    pr <- as.integer(c(0, cumsum(x@p[i+1] - x@p[i])))
    jr <- c(unlist(sapply(i, function(ind){x@j[if (x@p[ind] == x@p[ind+1])
                          {numeric(0)} else {(x@p[ind]+1):x@p[ind+1]}]})))
    xr <- c(unlist(sapply(i, function(ind){x@x[if (x@p[ind] == x@p[ind+1])
                          {numeric(0)} else {(x@p[ind]+1):x@p[ind+1]}]})))

    colIndMap <- rep(-1, j[length(j)])

    for (k in 1:length(j))
        colIndMap[j[k]] <- k - 1

    p1 <- sapply(1:length(i),
                 function(ind)
                 {
                    if (pr[ind] == pr[ind+1])
                        return(c())
                    else
                        return(jr[(pr[ind]+1):pr[ind+1]])
                 })

    if (is.matrix(p1))
    {
        p1 <- c(0,
                cumsum(apply(p1,2,
                             function(x){length(which(x %in% (j-1)))})))
    }
    else
    {
        p1 <- c(0,
                cumsum(sapply(p1,
                              function(x){length(which(x %in% (j-1)))})))
    }

    indices <- which(jr %in% (j-1))

    initialize(x, Dim=c(length(i), length(j)), factors=list(),
        p=as.integer(p1),
        j=as.integer(colIndMap[jr[indices]+1]),
        x=xr[indices],
        Dimnames=list(dimnames1, dimnames2),
        quadratic=x@quadratic, usedKernel=x@usedKernel)
}

setMethod("[", signature(x="ExplicitRepresentationSparse", i="index",
                         j="index", drop="missing"), subsetERSSampleFeature)
setMethod("[", signature(x="ExplicitRepresentationSparse", i="index",
                         j="index", drop="logical"), subsetERSSampleFeature)
setMethod("[", signature(x="ExplicitRepresentationSparse", i="index",
                         j="index", drop="ANY"), subsetERSSampleFeature)


###################################################
##
## PredictionProfile accessors / subsetting
##
###################################################

#' @rdname PredictionProfileAccessors
#' @title PredictionProfile Accessors
#' @name PredictionProfileAccessors
#' @aliases
#' sequences
#' sequences,PredictionProfile-method
#' profiles
#' profiles,PredictionProfile-method
#' baselines
#' baselines,PredictionProfile-method
#' [,PredictionProfile,index,ANY,ANY-method
#' @param object a prediction profile object
#'
#' @section Accessor-like methods:
#'
#' \describe{
#'   \item{}{\code{sequences}:
#'   return the sequences.
#'   }
#'   \item{}{\code{profiles}:
#'   return the prediction profiles.
#'   }
#'   \item{}{\code{baselines}:
#'   return the baselines.
#'   }
#'   \item{}{\code{x[i]}:
#'   return a \code{PredictionProfile} object that only contains the
#'   prediction profiles selected with the subsetting parameter \code{i}. This
#'   parameter can be a numeric vector with indices or a character vector with
#'   sample names.
#'   }
#' }
#' @return \code{sequences}: sequences for which profiles were generated\cr
#' \code{profiles}: prediction profiles\cr
#' \code{baselines}: baselines for the plot, this is the model offset\cr
#' distributed to all sequence positions
#' @examples
#' ## create kernel object for gappy pair kernel
#' gappy <- gappyPairKernel(k=1,m=11, annSpec=TRUE)
#' \dontrun{
#' ## load data
#' data(CCoil)
#'
#' ## perform training - feature weights are computed by default
#' model <- kbsvm(ccseq, yCC, gappya, pkg="LiblineaR", svm="C-svc", cost=15)
#'
#' ## compute prediction profiles
#' predProf <- getPredictionProfile(ccseq, gappya, 
#'                                  featureWeights(model),
#'                                  modelOffset(model))
#' predProf15 <- predProf[c(1,5),]
#' sequences(predProf15)
#' profiles(predProf15)
#' baselines(predProf15)
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

setMethod("sequences", "PredictionProfile", function(object) object@sequences)
setMethod("profiles", "PredictionProfile", function(object) object@profiles)
setMethod("baselines", "PredictionProfile", function(object) object@baselines)

setMethod("[", signature(x="PredictionProfile", i="index"),
    function(x, i)
    {
        if (is.character(i))
        {
            if (length(names(x@sequences)) < 1 ||
                any(is.na(names(x@sequences))))
                stop("missing names for subsetting\n")
            else
                i1 <- which(names(x@sequences) %in% i)

            if (length(i) != length(i1))
                stop("invalid names specified\n")

            i <- i1
        }
        else
        {
            ## convert negative subset
            if (all(i < 0))
                i <- (1:nrow(x@profiles))[i]
            else
            {
                if (min(i) < 1)
                    stop("subset indices must be all positive or",
                         " all negative\n")
            }

            if (min(i) < 1 || max(i) > nrow(x@profiles))
                stop("column subset must be between 1 and number",
                     " of sequences\n")
        }

        maxcol <- max(width(x@sequences[i]))

        initialize(x, sequences=x@sequences[i], kernel=x@kernel,
                   baselines=x@baselines[i],
                   profiles=x@profiles[i, 1:maxcol, drop=FALSE])
    }
)
