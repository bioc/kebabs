#2345678901234567890123456789012345678901234567890123456789012345678901234567890
## like in XString
toSeqSnippet <- function(x, width)
{
    if (width < 7L)
        width <- 7L

    seqlen <- nchar(x)

    if (seqlen <= width)
        x
    else
    {
        w1 <- (width - 2) %/% 2
        w2 <- (width - 3) %/% 2
        paste(substring(x, 1, w1), "...",
              substring(x, seqlen - w2 + 1, seqlen), sep="")
    }
}

show.BioVectorSeq <- function(object, index, wIndex, wWidth, wSeq,
                              wNames, withNames)
{
    seqLength <- width(object)[index]

    if (length(names(object)) > 0)
    {
        if (length(names(object)[index]) == 0)
            currName = "<NA>"
        else
        {
            if (nchar(names(object)[index]) > 20)
            {
                currName <- paste(substr(names(object)[index], 1, 17),
                                  "...", sep="")
            }
            else
                currName <- names(object)[index]
        }
    }
    else
        wSeq <- wSeq + wNames + 1

    cat(format(paste("[", index, "]", sep=""), width=wIndex,
               justify="right"), " ",
        format(seqLength, width=wWidth, justify="right"), " ", sep="")

    if (seqLength <= wSeq)
    {
        cat(format(object[[index]], width=wSeq,
                   justify="left"))
    }
    else
    {
        seqSnippet <- toSeqSnippet(object[index], width=wSeq)
        cat(format(seqSnippet, width=wSeq))
    }

    if (withNames)
        cat(" ", format(currName, width=wNames, justify="left"), sep="")

    cat("\n")

}

## same structure as show for XStringSet

#' @rdname show-methods
#' @title Display Various KeBABS Objects
#' @aliases
#' show
#' show,BioVector-method
#'
#' @description
#' Display methods for BioVector, SpectrumKernel, MismatchKernel,
#' GappyPairKernel, MotifKernel, SymmetricPairKernel,
#' ExplicitRepresentationDense, ExplicitRepresentationSparse,
#' PredictionProfile, CrossValidationResult, ModelSelectionResult,
#' SVMInformation and KBModel objects
#'
#' @param object object of class BioVector, PredictionProfile,
#' SpectrumKernel, MismatchKernel, GappyPairKernel,
#' MotifKernel, SymmetricPairKernel, ExplicitRepresentation,
#' ExplicitRepresentationSparse, PredictionProfile, CrossValidationResult,
#' ModelSelectionResult, SVMInformation or KBModel to be displayed
#'
#' @details
#' \code{show} displays on overview of the selected object.
#' @return \code{show}: show returns an invisible \code{NULL}
#'
#' @examples
#'
#' ## load coiled coil data
#' data(CCoil)
#'
#' ## show amino acid sequences
#' ccseq
#'
#' ## define spectrum kernel object
#' specK1 <- spectrumKernel(k=1, normalized=FALSE)
#'
#' ## show kernel object
#' show(specK1)
#'
#' ## compute explicit representation for the first 5 sequences 
#' ## in dense format
#' er <- getExRep(ccseq, specK1, sel=1:5, sparse=FALSE)
#'
#' ## show dense explicit representation
#' show(er)
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.

show.BioVector <- function(object)
{
    wNames <- 20L
    numRows <- 5L

    numSeqs <- length(object)
    withNames <- length(names(object)) > 0
    wIndex <- nchar(as.character(numSeqs)) +  2
    maxSeqLength <- max(nchar(object))
    wWidth <- max(nchar(maxSeqLength), nchar("width"))

    cat("  A", class(object), "instance of length", numSeqs, "\n")

    if (numSeqs > 0)
    {
        wSeq <- getOption("width") - wIndex - wWidth - wNames - 3L

        cat(format("width", width=wIndex + wWidth + 1, justify="right"),
            format("seq", width=wSeq, justify="left"))

        if (withNames)
            cat(" ", format("names", width=wNames, justify="left"), sep="")

        cat("\n")

        if (length(object) <= (2 * numRows + 1))
        {
            for (index in 1:length(object))
            {
                show.BioVectorSeq(object, index, wIndex, wWidth, wSeq,
                                  wNames, withNames)
            }
        }
        else
        {
            for (index in 1:numRows)
            {
                show.BioVectorSeq(object, index, wIndex, wWidth, wSeq,
                                  wNames, withNames)
            }

            cat(format("...", width=wIndex, justify="right"),
                format("...", width=wWidth, justify="right"),
                format("...", width=wSeq, justify="left"), "\n")

            for (index in (length(object) - numRows + 1L):length(object))
            {
                show.BioVectorSeq(object, index, wIndex, wWidth, wSeq,
                                  wNames, withNames)
            }
        }
    }
}

setMethod("show", signature(object="BioVector"), show.BioVector)

show.PredictionProfile <- function(object)
{
    maxCol <- 5
    noOfDigits <- 9
    colWidth <- noOfDigits + 3
    noOfBlocks <- 1
    blockSize <- nrow(object@profiles)

    cat("An object of class ", dQuote(class(object)), "\n\n")

    if (!is.null(object@sequences))
    {
        if (length(object@sequences) == 1)
            cat("Sequence:\n\n")
        else
            cat("Sequences:\n\n")

        show(object@sequences)
        cat("\n")
    }

    show(object@kernel)
    if (length(object@baselines) == 1)
        cat("\nBaseline: ", object@baselines, "\n\n")
    else if (length(object@baselines) <= 5)
        cat("\nBaselines: ", object@baselines, "\n\n")
    else
    {
        cat("\nBaselines: ", object@baselines[1:2], "  ...  ",
            object@baselines[(length(object@baselines) - 1):
                             length(object@baselines)], "\n\n")

    }

    if (nrow(object@profiles) == 1)
        cat("Profile:\n")
    else
        cat("Profiles:\n")

    if (nrow(object@profiles) > 10)
    {
        noOfBlocks <- 2
        blockSize <- 5
    }

    if (length(rownames(object@profiles)) > 0)
        nwidth <- min(max(nchar(rownames(object@profiles))), 20)
    else
        nwidth <- ceiling(log10(nrow(object@profiles))) + 2

    noPos <- ncol(object@profiles)
    offset <- 0

    if (nwidth > 17 && noPos > 4)
    {
        noOfDigits <- 8
        colWidth <- noOfDigits + 3
    }


    for (i in 1:noOfBlocks)
    {
        if (i == 2)
            offset <- nrow(object@profiles) - blockSize

        if (i == 1)
        {
            if (ncol(object@profiles) > 4)
            {
                cat(format("Pos 1", width=nwidth+1+colWidth, justify="right"),
                    format("Pos 2", width=colWidth, justify="right"),
                    format(paste("Pos", noPos-1), width=5+colWidth,
                           justify="right"),
                    format(paste("Pos", noPos), width=colWidth,
                           justify="right"),
                    "\n")
            }
            else
            {
                cat(format(" ", width=nwidth))

                for (j in 1:ncol(object@profiles))
                    cat(format(paste("Pos", j), width=colWidth+1,
                               justify="right"))

                cat("\n")
            }
        }


        for (j in (1 + offset):(blockSize + offset))
        {
            if (length(rownames(object@profiles)) > 0)
            {
                sampleName <- rownames(object@profiles)[j]

                if (nchar(sampleName) > 20)
                {
                    sampleName <- paste(substring(sampleName,1,17),
                                        "...", sep="")
                }
            }
            else
            {
                sampleName <- format(paste("[", j, "]", sep=""), nwidth,
                                     justify="right")
            }

            if (ncol(object@profiles) > 4)
            {
                cat(formatC(sampleName, format="s", width=nwidth),
                    formatC(object@profiles[j,1], format="f",
                            digits=noOfDigits, width=colWidth),
                    formatC(object@profiles[j,2], format="f",
                            digits=noOfDigits, width=colWidth),
                    " ...",
                    formatC(object@profiles[j,noPos-1], format="f",
                            digits=noOfDigits, width=colWidth),
                    formatC(object@profiles[j,noPos], format="f",
                            digits=noOfDigits,
                            width=colWidth))
                cat("\n")
            }
            else
            {
                cat(formatC(sampleName, format="s", width=nwidth))

                for (k in 1:ncol(object@profiles))
                {
                    cat(formatC(object@profiles[j,k], format="f",
                                digits=noOfDigits, width=colWidth+1))
                }

                cat("\n")
            }
        }

        if (i == 1 && noOfBlocks > 1)
        {
            cat(formatC(paste(rep(".", nwidth - 2), sep="", collapse=""),
                       format="s", width=nwidth))

            if (ncol(object@profiles) > 4)
            {
                cat(formatC(paste(rep(".", 6), sep="", collapse=""),
                            format="s", width=colWidth + 1),
                    formatC(paste(rep(".", 6), sep="", collapse=""),
                            format="s", width=colWidth),
                    "    ",
                    formatC(paste(rep(".", 6), sep="", collapse=""),
                            format="s", width=colWidth),
                    formatC(paste(rep(".", 6), sep="", collapse=""),
                            format="s", width=colWidth))
            }
            else
            {
                for (j in 1:ncol(object@profiles))
                {
                    cat(formatC(paste(rep(".", 6), sep="", collapse=""),
                                format="s", width=colWidth+1))
                }
            }
            cat("\n")
        }
    }

    cat("\n")
}

#' @rdname show-methods
#' @aliases
#' show,PredictionProfile-method
#'

setMethod("show", signature(object="PredictionProfile"),
          show.PredictionProfile)

#' @rdname show-methods
#' @aliases
#' show,SpectrumKernel-method
#'

setMethod("show", signature(object="SpectrumKernel"),
    function(object)
    {
        cat("Spectrum Kernel: ")
        cat(paste("k=", object@k, sep=""))
        if (object@r > 1)
            cat(paste(", r=", object@r, sep=""))
        if (object@annSpec == TRUE)
            cat(paste(", annSpec=TRUE"))
        if (length(object@distWeight) > 0)
        {
             dwString <- distWeightKernelToString(object@distWeight)
             cat(", distWeight=", dwString, sep="")
        }
        if (object@normalized == FALSE)
            cat(", normalized=FALSE")
        if (object@exact == FALSE)
            cat(", exact=FALSE")
        if (object@ignoreLower == FALSE)
            cat(", ignoreLower=FALSE")
        if (object@presence == TRUE)
            cat(", presence=TRUE")
        if (object@revComplement == TRUE)
            cat(", revComplement=TRUE")
        cat("\n")
    }
)

#' @rdname show-methods
#' @aliases
#' show,MismatchKernel-method
#'

setMethod("show", signature(object="MismatchKernel"),
    function(object)
    {
        cat("Mismatch Kernel: ")
        cat(paste("k=", object@k, sep=""))
        cat(paste(", m=", object@m, sep=""))
        if (object@r > 1)
            cat(paste(", r=", object@r, sep=""))
        if (object@normalized == FALSE)
            cat(paste(", normalized=FALSE"))
        if (object@exact == FALSE)
            cat(", exact=FALSE")
        if (object@ignoreLower == FALSE)
            cat(", ignoreLower=FALSE")
        if (object@presence == TRUE)
            cat(", presence=TRUE")
        cat("\n")
    }
)

#' @rdname show-methods
#' @aliases
#' show,MotifKernel-method
#'

setMethod("show", signature(object="MotifKernel"),
    function(object)
    {
        cat("Motif Kernel:\n\n")
        cat("Motifs:\n")

        obj <- kernelParameters(object)

        if (length(obj$motifs) < 11)
        {
            for (i in 1:length(obj$motifs))
            {
                numChars <- nchar(obj$motifs[i])
                if (numChars <= 70)
                    cat(obj$motifs[i], "\n")
                else
                {
                    cat(substr(obj$motifs[i],1,30), "   ....   ",
                        substr(obj$motifs[i], numChars-29,
                               numChars), "\n")
                }
            }
            cat("\n")
        }
        else
        {
            for (i in 1:5)
            {
                numChars <- nchar(obj$motifs[i])
                if (numChars <= 70)
                    cat(obj$motifs[i], "\n")
                else
                {
                    cat(substr(obj$motifs[i], 1, 30), "  . . .  ",
                        substr(obj$motifs[i], numChars - 29,
                               numChars), "\n")
                }
            }
            cat(" . . .  \n")
            cat(" . . .  \n")
            cat(" . . .  \n")
            for (i in (length(obj$motifs) - 4):length(obj$motifs))
            {
                numChars <- nchar(obj$motifs[i])
                if (numChars <= 70)
                    cat(obj$motifs[i], "\n")
                else
                {
                    cat(substr(obj$motifs[i], 1, 30), "  . . .  ",
                        substr(obj$motifs[i], numChars - 29,
                               numChars), "\n")
                }
            }
            cat("\n")
        }
        cat("Kernel Parameters:\n")
        suppressColon = TRUE
        if (object@r > 1)
        {
            if (suppressColon)
                suppressColon = FALSE
            else
                cat(", ")
            cat(paste("r=", object@r, sep=""))
        }
        if (object@annSpec == TRUE)
        {
            if (suppressColon)
                suppressColon = FALSE
            else
                cat(", ")
            cat(paste("annSpec=TRUE"))
        }
        if (length(object@distWeight) > 0)
        {
            if (is.numeric(object@distWeight))
            {
                if (length(object@distWeight) == 1)
                {
                    if (suppressColon)
                        suppressColon = FALSE
                    else
                        cat(", ")
                    cat(paste("distWeight=",
                              object@distWeight, sep=""))
                }
                else
                {
                    if (suppressColon)
                        suppressColon = FALSE
                    else
                        cat(", ")
                    cat(paste("distWeight=",
                              paste("c(",paste(object@distWeight,
                                               collapse=","),")",
                                    sep=""), sep=""))
                }
            }
            else
            {
                if (suppressColon)
                    suppressColon = FALSE
                else
                    cat(", ")
                cat("\ndistWeight=")
                cat(format(object@distWeight))
            }
        }
        if (object@normalized == FALSE)
        {
            if (suppressColon)
                suppressColon = FALSE
            else
                cat(", ")
            cat("normalized=FALSE")
        }
        if (object@exact == FALSE)
        {
            if (suppressColon)
                suppressColon = FALSE
            else
                cat(", ")
            cat("exact=FALSE")
        }
        if (object@ignoreLower == FALSE)
        {
            if (suppressColon)
                suppressColon = FALSE
            else
                cat(", ")
            cat("ignoreLower=FALSE")
        }
        if (object@presence == TRUE)
        {
            if (suppressColon)
                suppressColon = FALSE
            else
                cat(", ")
            cat("presence=TRUE")
        }
        cat("\n")
    }
)

#' @rdname show-methods
#' @aliases
#' show,GappyPairKernel-method
#'

setMethod("show", signature(object="GappyPairKernel"),
    function(object)
    {
        cat("gappy pair kernel: ")
        cat(paste("k=", object@k, sep=""))
        cat(paste(", m=", object@m, sep=""))
        if (object@r > 1)
            cat(paste(", r=", object@r, sep=""))
        if (object@annSpec == TRUE)
            cat(paste(", annSpec=TRUE"))
        if (length(object@distWeight) > 0)
        {
            if (is.numeric(object@distWeight))
            {
                if (length(object@distWeight) == 1)
                    {
                    cat(paste(", distWeight=",
                              object@distWeight, sep=""))
                    }
                else
                {
                    cat(paste(", distWeight=",
                              paste("c(",paste(object@distWeight,
                                               collapse=","),")",
                            sep=""), sep=""))
                }
            }
            else
            {
                cat(", \ndistWeight=")
                cat(format(object@distWeight))
            }
        }
        if (object@normalized == FALSE)
            cat(", normalized=FALSE")
        if (object@exact == FALSE)
            cat(", exact=FALSE")
        if (object@ignoreLower == FALSE)
            cat(", ignoreLower=FALSE")
        if (object@presence == TRUE)
            cat(", presence=TRUE")
        if (object@revComplement == TRUE)
            cat(", revComplement=TRUE")
        cat("\n")
    }
)

#' @rdname show-methods
#' @aliases
#' show,SymmetricPairKernel-method
#'

setMethod("show", signature(object="SymmetricPairKernel"),
    function(object)
    {
        cat("Symmetric Pair Kernel: ")
        cat("kernelType=", object@kernelType)
        if (object@r > 1)
            cat(paste(", r=", object@r, sep=""))
        cat("\n  Single Instance Kernel:\n    ")
        show(object@siKernel)
    }
)

#' @rdname show-methods
#' @aliases
#' show,ExplicitRepresentationDense-method
#'

setMethod("show", signature(object="ExplicitRepresentationDense"),
    function(object)
    {
        cat("Dense explicit representation of class ",
            dQuote(class(object)), "\n\n")

        cat(paste("Quadratic             :", object@quadratic))
        cat("\n")
        if (length(object@usedKernel) > 0)
        {
            cat(paste("Used kernel           : \n   "))
            cat(show(object@usedKernel))
            cat("\n")
        }

        show(object@.Data)
        cat("\n")
    }
)

#' @rdname show-methods
#' @aliases
#' show,ExplicitRepresentationSparse-method
#'

setMethod("show", signature(object="ExplicitRepresentationSparse"),
    function(object)
    {
        cat("Sparse explicit representation of class ",
            dQuote(class(object)), "\n\n")

        cat(paste("Quadratic             :", object@quadratic))
        cat("\n")
        if (length(object@usedKernel) > 0)
        {
            cat(paste("Used kernel           : \n   "))
            cat(show(object@usedKernel))
            cat("\n")
        }

        show(as(object, "dgRMatrix"))
        cat("\n")
    }
)

#' @rdname show-methods
#' @aliases
#' show,CrossValidationResult-method
#'

setMethod("show", signature(object="CrossValidationResult"),
    function(object)
    {
        if (object@outerCV)
        {
            cat("\nOuter cross validation result object of class ",
                dQuote(class(object)), "\n\n")
            cat(paste("nestedCross           :", object@cross, "\n"))
            cat(paste("noNestedCross         :", object@noCross, "\n"))
        }
        else
        {
            cat("\nCross validation result object of class ",
                dQuote(class(object)), "\n\n")
            cat(paste("cross                 :", object@cross, "\n"))
            cat(paste("noCross               :", object@noCross, "\n\n"))
        }

        if (object@noCross > 1)
        {
            cat(paste("Mean CV error         :",
                      format(mean(object@cvError), digits=8), "\n"))

            if ("ACC" %in% object@perfParameters)
            {
                cat(paste("Mean accuracy         :",
                          format(mean(object@ACC), digits=8), "\n"))
            }

            if ("BACC" %in% object@perfParameters)
            {
                cat(paste("Mean balanced accuracy:",
                          format(mean(object@BACC), digits=8), "\n"))
            }

            if ("MCC" %in% object@perfParameters)
            {
                cat(paste("Mean Matthews CC      :",
                          format(mean(object@MCC), digits=8), "\n"))
            }

            if ("AUC" %in% object@perfParameters)
            {
                cat(paste("Area under the curve  :",
                          format(mean(object@AUC), digits=8), "\n"))
            }

            cat("\n")
            cat(paste("Run CV errors         :\n"))
            cat(paste(format(object@cvError, digits=8), collapse=", "))
            cat("\n\n")

            if ("ACC" %in% object@perfParameters)
            {
                cat(paste("Run accuracies        :\n"))
                cat(paste(format(object@ACC, digits=8), collapse=", "))
                cat("\n\n")
            }

            if ("BACC" %in% object@perfParameters)
            {
                cat(paste("Run bal. accuracies   :\n"))
                cat(paste(format(object@BACC, digits=8), collapse=", "))
                cat("\n\n")
            }

            if ("MCC" %in% object@perfParameters)
            {
                cat(paste("Run Matthews CC       :\n"))
                cat(paste(format(object@MCC, digits=8), collapse=", "))
                cat("\n\n")
            }

            if ("AUC" %in% object@perfParameters)
            {
                cat(paste("Run AUCs              :\n"))
                cat(paste(format(object@AUC, digits=8), collapse=", "))
                cat("\n\n")
            }
        }
        else
        {
            cat(paste("CV error:             :",
                      format(object@cvError, digits=8), collapse=", "), "\n")

            if ("ACC" %in% object@perfParameters)
            {
                cat(paste("Accuracy:             :",
                          format(object@ACC, digits=8), collapse=", "), "\n")
            }

            if ("BACC" %in% object@perfParameters)
            {
                cat(paste("Balanced accuracy:    :",
                          format(object@BACC, digits=8), collapse=", "), "\n")
            }

            if ("MCC" %in% object@perfParameters)
            {
                cat(paste("Matthews CC           :",
                          format(object@MCC, digits=8), collapse=", "), "\n")
            }

            if ("AUC" %in% object@perfParameters)
            {
                cat(paste("Area under the curve  :",
                          format(mean(object@AUC), digits=8), "\n"))
            }

            cat("\n")
        }
    }
)

#' @rdname show-methods
#' @aliases
#' show,ModelSelectionResult-method
#'

setMethod("show", signature(object="ModelSelectionResult"),
    function(object)
    {
        if (object@nestedCross > 0)
        {
            cat("\nModel selection result object of class ",
                dQuote(class(object)), "\n\n")

            cat(paste("cross                 :", object@cross, "\n"))
            cat(paste("noCross               :", object@noCross, "\n"))
            cat(paste("nestedCross           :", object@nestedCross, "\n"))
            cat(paste("noNestedCross         :", object@noNestedCross))
        }
        else
        {
            cat("\nGrid search result object of class ",
                dQuote(class(object)), "\n\n")

            cat(paste("cross                 :", object@cross, "\n"))
            cat(paste("noCross               :", object@noCross))
            ## the field smallestCVError always contains the best value
            ## as defined in the performance objective
            if (object@perfObjective == "MCC")
            {
                cat(paste("\nBest MCC value        :",
                          format(object@smallestCVError, digits=8)))
            }
            else if (object@perfObjective == "BACC")
            {
                cat(paste("\nBest bal. accuracy    :",
                          format(object@smallestCVError, digits=8)))
            }
            else if (object@perfObjective == "AUC")
            {
                cat(paste("\nBest AUC              :",
                format(object@smallestCVError, digits=8)))
            }
            else
            {
                cat(paste("\nSmallest CV error     :",
                          format(object@smallestCVError, digits=8)))
            }
        }

        cat("\n\n")

        if (length(object@groupBy) > 0)
        {
            cat(paste("groupBy:\n"))

            if (length(object@groupBy) < 100)
                cat(paste(object@groupBy, collapse=","))
            else
            {
                cat(paste(object@groupBy[1:3], collapse=","), "...",
                    paste(object@groupBy[(length(object@groupBy)-2):
                          length(object@groupBy)],
                          collapse=","))
            }

            cat("\n\n")
        }

        if (length(object@gridRows) > 0)
        {
            cat("Grid Rows:\n")
            namesGridRows <- names(object@gridRows)
            maxNameLength <- max(nchar(namesGridRows))+3

            for (i in 1:length(object@gridRows))
            {
                cat(paste(namesGridRows[i],
                          paste(rep(" ", 11-nchar(namesGridRows[i])),
                                collapse="")))
                cat(seqKernelAsChar(object@gridRows[[i]]), "\n")
            }

            cat("\n")
        }

        if (!is.null(object@gridCols))
        {
            cat("Grid Columns:\n")
            print(object@gridCols)
            cat("\n")
        }

        if (object@nestedCross == 0)
        {
            if (object@perfObjective == "MCC")
            {
                cat("Grid MCC values:\n\n")
                print(object@gridMCC)
                cat ("\n")
            }
            else if (object@perfObjective == "BACC")
            {
                cat("Grid Balanced Accuracies:\n\n")
                print(object@gridBACC)
                cat ("\n")
            }
            else if (object@perfObjective == "AUC")
            {
                cat("Grid AUCs:\n\n")
                print(object@gridAUC)
                cat ("\n")
            }
            else
            {
                cat("Grid Errors:\n\n")
                print(object@gridErrors)
                cat ("\n")
            }
        }

        if (length(object@gridRows) > 1)
        {
            cat("Selected Grid Row:\n")
            if (is.list(object@selGridRow))
            {
                for (i in 1:length(object@selGridRow))
                    cat(i, "",seqKernelAsChar(object@selGridRow[[i]]), "\n")

                cat("\n")
            }
            else
                cat(seqKernelAsChar(object@selGridRow), "\n\n")
        }

        if (length(object@gridCols) > 0 && nrow(object@gridCols) > 1)
        {
            cat("Selected Grid Column:\n")

            if (is.numeric(object@selGridCol))
                cat(names(object@selGridCol), "=", object@selGridCol)
            else
                print(object@selGridCol)
        }

        cat("\n\n")
    }
)

#' @rdname show-methods
#' @aliases
#' show,SVMInformation-method
#'

setMethod("show", signature(object="SVMInformation"),
    function(object)
    {
        cat("\nSVM info object of class ", dQuote(class(object)), "\n\n")

        if (length(object@reqKernel) > 0)
        {
            cat(paste("Kernel                : "))
            cat(show(object@reqKernel))
            cat("\n")
        }

        cat(paste("Available Packages    :", paste(object@availPackages,
                                                   collapse=", "), "\n"))
        cat(paste("Package               :", object@selPackage, "\n"))
        cat(paste("Classifier            :", object@selSVM, "\n"))

        if (length(object@selSVMPar) > 0)
        {
            svmPars <- object@selSVMPar
            cat(paste("Classifier Parameters : "))
            cat(names(svmPars)[[1]], "=", svmPars[[1]])
            if (length(svmPars) > 1)
            {
                for (i in 2:length(svmPars))
                    cat(", ", names(svmPars)[[i]], "=", svmPars[[i]])
            }
            cat("\n")
        }

        if (object@selPackage != object@reqPackage)
        {
            cat(paste("Requested Package     :", object@reqPackage, "\n"))
            cat(paste("Requested Classifier  :", object@reqSVM, "\n"))
            if (length(object@reqSVMPar) > 0)
            {
                svmPars <- object@reqSVMPar
                cat(paste("Classifier Parameters : "))
                cat(names(svmPars)[[1]], "=", svmPars[[1]])
                if (length(svmPars) > 1)
                {
                    for (i in 2:length(svmPars))
                        cat(", ", names(svmPars)[[i]], "=", svmPars[[i]])
                }
                cat("\n")
            }
        }

        if (object@selExplicit)
        {
            cat(paste("Explicit Kernel       :",
                      object@svmInfo@explicitKernel, "\n"))
            if (object@svmInfo@reqExplicit != "yes")
            {
                cat(paste("Requested explicit    :",
                          object@svmInfo@reqExplicit, "\n"))
            }
        }
        cat("\n")
    }
)

#' @rdname show-methods
#' @aliases
#' show,KBModel-method
#'

setMethod("show", signature(object="KBModel"),
    function(object)
    {
        cat("\nKEBABS result object of class \"KBModel\"","\n\n")
        cat(paste("Number of samples     :", object@numSequences, "\n"))
        if (length(object@svmInfo@reqKernel) > 0)
        {
            cat(paste("Kernel                : "))
            cat(show(object@svmInfo@reqKernel))
            cat("\n")
        }

        cat(paste("Available Packages    :",
                  paste(object@svmInfo@availPackages, collapse=", "), "\n"))

        if (length(object@svmInfo@selPackage) > 1)
        {
            cat(paste("Packages              :",
                      paste(object@svmInfo@selPackage, collapse=", "), "\n"))
        }
        else
        {
            cat(paste("Package               :", object@svmInfo@selPackage,
                      "\n"))
        }

        if (length(object@svmInfo@selSVM) > 1)
        {
            cat(paste("SVMs                  :",
                      paste(object@svmInfo@selSVM, collapse=", "), "\n"))
        }
        else
        {
            cat(paste("SVM                   :", object@svmInfo@selSVM,
                      "\n"))
        }

        if (length(object@svmInfo@selSVMPar) > 0)
        {
            svmPars <- object@svmInfo@selSVMPar
            cat(paste("Classifier Parameters : "))
            cat(names(svmPars)[[1]], "=", svmPars[[1]])
            if (length(svmPars) > 1)
            {
                for (i in 2:length(svmPars))
                cat(", ", names(svmPars)[[i]], "=", svmPars[[i]])
            }
            cat("\n")
        }

        if (!(all(object@svmInfo@selPackage == object@svmInfo@reqPackage)))
        {
            if (length(object@svmInfo@reqPackage) > 1)
            {
                cat(paste("Requested Packages    :",
                      paste(object@svmInfo@reqPackage, collapse=", "), "\n"))
            }
            else
            {
                cat(paste("Requested Package     :",
                          object@svmInfo@reqPackage, "\n"))
            }

            if (length(object@svmInfo@reqSVM) > 1)
            {
                cat(paste("Requested SVMs        :",
                      paste(object@svmInfo@reqSVM, collapse=", "), "\n"))
            }
            else
            {
                cat(paste("Requested SVM         :",
                          object@svmInfo@reqSVM, "\n"))
            }

            if (length(object@svmInfo@reqSVMPar) > 0)
            {
                svmPars <- object@svmInfo@reqSVMPar
                cat(paste("Classifier Parameters : "))
                cat(names(svmPars)[[1]], "=", svmPars[[1]])
                if (length(svmPars) > 1)
                {
                    for (i in 2:length(svmPars))
                    cat(", ", names(svmPars)[[i]], "=", svmPars[[i]])
                }
                cat("\n")
            }
        }

        if (object@svmInfo@selExplicit)
        {
            cat(paste("Explicit Kernel       :",
                      object@svmInfo@explicitKernel, "\n"))
            if (object@svmInfo@reqExplicit != "yes")
            {
                cat(paste("Requested explicit    :",
                          object@svmInfo@reqExplicit, "\n"))
            }
        }

        if (object@svmInfo@probModel)
        {
            if (!is.na(object@sigma))
                cat(paste("Laplace distr. width:", object@sigma, "\n"))
        }

        cat("\n")

        cat("Call:\n")
        cat(paste(object@call, "\n"))
        cat("\n")
        cat("Classifier specific model    :\n\n")
        cat(show(object@svmModel))
        cat("\n")
    }
)

#' @rdname show-methods
#' @aliases
#' show,ROCData-method
#'

setMethod("show", signature(object="ROCData"),
    function(object)
    {
        cat("ROC Data of class ",
        dQuote(class(object)), "\n\n")
    
        cat(paste("AUC:", round(object@AUC, 3)))
        cat("\n\n")
    }
)

