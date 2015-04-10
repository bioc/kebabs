##345678901234567890123456789012345678901234567890123456789012345678901234567890
getSlidingWindowAverage <- function(x, windowSize)
{
    if ((windowSize %% 2) == 0)
        windowSize <- windowSize + 1

    result <- rep(0, length(x))

    for (i in 1: length(x))
    {
        lower <- max(1, i - floor(windowSize/2))
        upper <- min(length(x), i + floor(windowSize/2))
        result[i] <- sum(x[lower:upper]) / (upper - lower)
    }

    result
}

plotPredictionProfile.Missing <- function(x, sel=NULL, col=c("red", "blue"),
                    standardize=TRUE, shades=NULL, legend="default",
                    legendPos="topright", ylim=NULL, xlab="", ylab="weight",
                    lwd.profile=1, lwd.axis=1, las=1, heptads=FALSE,
                    annotate=FALSE, markOffset=TRUE, windowSize=1, ...)
{
    if (length(x@sequences) < 1)
        stop("no sequence in prediction profiles\n")

    if (nrow(x@profiles) < 1)
        stop("profile data not available\n")

    if ((is.null(sel) && length(x@sequences) > 2) ||
        (!is.null(sel) && length(sel) > 2))
        stop("too many sequences selected for profile plot\n")

    if (!is.null(sel))
    {
        if (!is.numeric(sel) || min(sel) < 1 || max(sel) > nrow(x@profiles) ||
            length(sel) > 2)
            stop("parameter 'sel' must contain up to two indices to profiles\n")
    }
    else
        sel <- c(1:min(2, nrow(x@profiles)))

    if (!is.null(shades) &&
        (!is.numeric(try(col2rgb(shades[1]), silent=TRUE)[1]) ||
         !is.numeric(try(col2rgb(shades[2]), silent=TRUE)[1])))
        stop("argument `shades' must be a vector of two colors\n")

    if (!isTRUEorFALSE(heptads))
        stop("argument `heptads' must be a logical value\n")

    if (!isTRUEorFALSE(annotate))
        stop("argument `annotate' must be a logical value\n")

    if (!isTRUEorFALSE(markOffset))
        stop("argument `markOffset' must be a logical value\n")

    if (!is.numeric(windowSize) || length(windowSize) != 1)
        stop("argument `windowSize' must be a single integer value\n")

    if (!is.infinite(windowSize))
        windowSize <- as.integer(windowSize)

    if (!is.null(legend))
    {
        if (length(legend) == 1 && legend[1] == "default")
        {
            legend <- rownames(x@profiles)[sel]
            
            if (is.null(legend))
                legend <- as.character(sel)
        }
    }

    doubleProf <- length(sel) == 2

    if (inherits(x@kernel, "SequenceKernel") && !isUserDefined(x@kernel) &&
        kernelParameters(x@kernel)$annSpec ==TRUE)
        annot <- mcols(x@sequences[sel])[["annotation"]]
    else
        annot <- NULL
    
    if (inherits(x@kernel, "SequenceKernel") && !isUserDefined(x@kernel) &&
        length(kernelParameters(x@kernel)$distWeight) > 0)
        offset <- mcols(x@sequences[sel])[["offset"]]
    else
        offset <- NULL
    
    shift <- 0
    
    if (!is.null(offset))
    {
        offsetAll <- mcols(x@sequences)[["offset"]]
        startPos <- -offsetAll
        startPos[which(startPos > 0)] <- 0
        startPos <- startPos + 1
        minPos <- min(startPos)

        if (doubleProf)
        {
            shift <- offset[1] - offset[2]
            start1 <- startPos[sel[1]] - minPos + 1
            start2 <- startPos[sel[2]] - minPos + 1
        }
        else
        {
            start1 <- startPos[sel[1]] - minPos + 1
            start2 <- 0
        }
    }
    else
    {
        markOffset <- FALSE
        start1 <- 1
        start2 <- 1
    }

    if (!doubleProf)
    {
        if (!is.numeric(try(col2rgb(col), silent=TRUE)[1]))
            stop("argument `col' must be a string denoting a color\n")

        if (!isSingleString(legend))
            stop("argument `legend' must be a string\n")
    }
    else
    {
        if (!is.null(col) &&
            (!is.numeric(try(col2rgb(col[1]), silent=TRUE)[1]) ||
             !is.numeric(try(col2rgb(col[2]), silent=TRUE)[1])))
            stop("argument `col' must be a vector of two colors\n")

        if (!is.null(legend) && (!is.character(legend) || length(legend) < 2))
            stop("argument `legend' must be a vector of two strings\n")

        if (width(x@sequences[sel[1]]) != width(x@sequences[sel[2]]) &&
            standardize == FALSE)
        {
            stop("unstandardized plot not possible for sequences of\n",
                 "        unequal length\n")
        }
    }

    length1 <- width(x@sequences[sel[1]])

    if (doubleProf)
        length2 <- width(x@sequences[sel[2]])
    else
        length2 <- 0

    n <- max(width(x@sequences[sel]))

    prfl1 <- x@profiles[sel[1],start1:(start1 + length1 - 1)]

    if (doubleProf)
        prfl2 <- x@profiles[sel[2],start2:(start2 + length2 - 1)]
    else
        prfl2 <- NULL

    if (standardize)
    {
        prfl1 <- prfl1 - x@baselines[sel[1]]
        hlin <- 0

        if (doubleProf)
            prfl2 <- prfl2 - x@baselines[sel[2]]
    }
    else
        hlin <- x@baselines[sel[1]]

    if (is.infinite(windowSize))
        prfl1 <- cumsum(prfl1)
    else if (windowSize != 1)
        prfl1 <- getSlidingWindowAverage(prfl1, windowSize)

    if (doubleProf)
    {
        if (is.infinite(windowSize))
            prfl2 <- cumsum(prfl2)
        else if (windowSize != 1)
            prfl2 <- getSlidingWindowAverage(prfl2, windowSize)
    }

    if (is.null(ylim))
    {
        ylim <- rep(1.05 * max(abs(c(prfl1, prfl2, hlin))), 2)
        ylim[1] <- - ylim[1]
    }

    plot(x=NULL,y=NULL, xlim=c(1, n + 1 + abs(shift)), ylim=ylim,
         axes=FALSE, xlab=xlab, ylab=ylab, type="s", ...)

    if (!is.null(shades))
    {
        if (is.vector(shades) && length(shades) >= 2 &&
            is.character(shades[1]) && is.character(shades[2]))
        {
            polygon(c(1, n + 1 + abs(shift), n + 1 + abs(shift), 1),
                    c(hlin, hlin, ylim[2], ylim[2]), col=shades[1], border=NA)
            polygon(c(1, n + 1 + abs(shift), n + 1 + abs(shift), 1),
                    c(hlin, hlin, ylim[1], ylim[1]), col=shades[2], border=NA)
        }
        else
            stop("shades must be a vector of 2 colors\n")
    }

    lines(x=c(1, n + 1 + abs(shift)), y=c(hlin, hlin), col="lightgray")

    seq1 <- strsplit(as.character(x@sequences[sel[1]]), "")[[1]]
    seq2 <- seq1
    shift1 <- 0
    shift2 <- 0

    if (doubleProf)
    {
        seq2 <- strsplit(as.character(x@sequences[sel[2]]), "")[[1]]

        if (shift < 0)
            shift1 <- abs(shift)
        else
            shift2 <- abs(shift)

        seq1 <- c(rep(" ", shift1), seq1,
                  rep(" ", n + abs(shift) - shift1 - length(seq1)))
        seq2 <- c(rep(" ", shift2), seq2,
                  rep(" ", n + abs(shift) - shift2 - length(seq2)))
    }

    matchS <- (seq1 == seq2)

    if (length(annot) > 0 && heptads == TRUE)
    {
        sapply(gregexpr("a", annot[1])[[1]],
               function(i)
               if (i > 0)
               {
                   lines(c(i + shift1, i + shift1), ylim,
                         col="lightgray")
               })

        if (!is.null(heptads) && (substr(annot[1], n, n) == "g"))
        {
            lines(c(n + 1 + shift1, n + 1 + shift1), ylim,
                  col="lightgray")
        }

        sapply(gregexpr("(a[^b]|b[^c]|c[^d]|d[^e]|e[^f]|f[^g]|g[^a])",
                        annot[1], perl=TRUE)[[1]],
               function(i)
               if (i > 0)
               {
                   lines(c(i + 1+ shift1, i + 1 + shift1), ylim,
                         col="red")
               })
    }

    if (doubleProf)
    {
        lines(c(1:(length2 + 1)) + shift2, c(prfl2[1:length2], prfl2[length2]),
              type="s", col=col[2], lwd=lwd.profile, ...)
    }

    lines(c(1:(length1 + 1)) + shift1, c(prfl1[1:length1], prfl1[length1]),
          type="s", col=col[1], lwd=lwd.profile, ...)

    axis(side=2, lwd=lwd.axis, las=las, ...)

    if (markOffset)
    {
        if (doubleProf && (shift2 + offset[2] > 0))
        {
            mtext(side=1, at=(shift2 + offset[2] + 0.5), line=-1.3, text="v",
                  col=col[2], cex=1.2)
        }

        if (shift1 + offset[1] > 0)
        {
            mtext(side=3, at=(shift1 + offset[1] + 0.5), line=-1.3, text="^",
                  col=col[1], cex=1.3)
        }
    }

    mtext(side=1, at=(1:(n + abs(shift)) + 0.5), line=0, text=seq2,
          col=ifelse(matchS, "black", col[2]), ...)
    mtext(side=3, at=(1:(n + abs(shift)) + 0.5), line=0, text=seq1,
          col=ifelse(matchS, "black", col[1]), ...)


    if (length(annot) > 0 && annotate == TRUE)
        text(1:n + 0.5 + shift1, 0, strsplit(annot[1], "")[[1]], adj=c(0.5,0))

    if (!is.null(legend))
    {
        if (length(legendPos) > 1)
        {
            legend(x=legendPos[1], y=legendPos[2], col=col,
                   legend=legend, lwd=lwd.profile, bg="white")
        }
        else
        {
            legend(x=legendPos[1], col=col,
                   legend=legend, lwd=lwd.profile, bg="white")
        }
    }
}

#' @rdname plot-methods
#' @aliases
#' plot
#' plot,PredictionProfile-method
#' @title Plot Prediction Profiles, Cross Validation Result, Grid Search
#' Performance Parameters and Receiver Operating Characteristics
#'
#' @description Functions for visualizing prediction profiles,
#' cross validation result, grid search performance parameters and
#' receiver operating characteristics
#'
#' @param x for the first method above a prediction profile object of class
#' \code{\link{PredictionProfile}} containing the profiles to be plotted,
#' for the second method a cross validation result object usually taken
#' from the trained kebabs model object
#'
#' @param sel an integer vector with one or two entries to select samples
#' of the prediction profile matrix for plotting, if this parameter is not
#' supplied by the user the frist one or two samples are selected.
#'
#' @param col a character vector with one or two color names used for plotting
#' the samples. Default=c("red", "blue").
#'
#' @param ylim argument that allows the user to preset the y-range of the
#' profile plot.
#'
#' @param standardize logical. If \code{FALSE}, the profile values \code{s_i}
#' are displayed as they are with the value \eqn{y=-b/L} superimposed as
#' a light gray line. If \code{TRUE} (default), the whole profile is
#' shifted by \eqn{-b/L} and the light gray line is displayed at \code{y=0}.
#'
#' @param windowSize length of sliding window. When the parameter is set to
#' the default value 1 the contributions of each position are plotted
#' as step function. For kernels with multiple patterns at one position
#' (mismatch, gappy pair and motif kernel) the weight contributions of all
#' patterns at the position are summed up. Values larger than 1 define the
#' length of a sliding window. All contributions within the window are
#' averaged and the resulting value is displayed at the center position of
#' the window. For positions within half of the window size from the start
#' and end of the sequence the averaging cannot be performed over the full
#' window but just the remaining positions. This means that the variation
#' of the averaged weight contributions is higher in these border regions.
#' If an even value is specified for this parameter one is added to the
#' parameter value. When the parameter is set to \code{Inf} (infinite)
#' instead of averages cumulative values along the sequence are used, i.e.
#' at each position the sum of all contributions up to this position is
#' displayed. In this case the plot shows how the standardized or
#' unstandardized value (see parameter \code{standardize}) of the
#' discrimination function builds up along the sequence. Default=1
#'
#' @param shades vector of at least two color specifications; If not NULL,
#' the background area above and below the base line \code{y=-b/L} are shaded
#' in colors \code{shades[1]} and \code{shades[2]}, respectively. Default=NULL
#'
#' @param legend a character vector with one or two character strings
#' containing the legend/description of the profile. If empty, no legend is
#' displayed.
#'
#' @param legendPos position specification for the legend(if \code{legend} is
#' specified). Can either be a vector with coordinates or a single keyword like
#' \dQuote{topright} (see \code{\link[graphics:legend]{legend}}).
#'
#' @param xlab label of horizontal axis, empty by default.
#'
#' @param ylab label of vertical axis, defaults to "weight".
#'
#' @param heptads logical indicating whether for proteins with heptad
#' annotation (i.e. characters a to g, usually in periodic repetition) the
#' heptad structure should be indicated through vertical lightgray lines
#' each heptad. Default=FALSE
#'
#' @param annotate logical indicating whether annotation information should
#' be shown in the center of the plot; Default=FALSE
#'
#' @param markOffset logical indicating whether the start positions in the
#' sequences according to the assigned offset elmement metadata values should
#' be shown near the sequence characters; for the upper sequence the first
#' position is marked by \code{"^"} below the respective character, for the
#' lower sequence it is marked by \code{"v"} above the sequence. If no offset
#' element metadata is assigned to the sequences the marks are suppressed.
#' Default=TRUE
#'
#' @param aucDigits number of decimal places of AUC to be printed into
#' the ROC plot. If this parameter is set to 0 the AUC will not be added to
#' the plot. Default=3
#'
#' @param lwd.profile profile line width as described for parameter \code{lwd}
#' in \code{\link[graphics:par]{par}}
#'
#' @param lwd.axis axis line width as described for parameter \code{lwd}
#' in \code{\link[graphics:par]{par}}
#'
#' @param las see \code{\link[graphics:par]{par}}
#'
#' @param lwd see \code{\link[graphics:par]{par}}
#'
#' @param cex see \code{\link[graphics:mtext]{mtext}}
#'
#' @param side see \code{\link[graphics:mtext]{mtext}}
#'
#' @param line see \code{\link[graphics:mtext]{mtext}}
#'
#' @param adj see \code{\link[graphics:mtext]{mtext}}
#'
#' @param ... all other arguments are passed to the standard
#' \code{\link[graphics:plot]{plot}} command that is called internally to
#' display the graphics window.
#'
#' @details
#' Plotting of Prediction Profiles\cr\cr
#' The first variant of the \code{plot} method mentioned in the usage section
#' displays one or two prediction profiles as a step
#' function with the steps connected by vertical lines. The parameter
#' \code{sel} allows to select the sample(s) if the prediction profile object
#' contains the profiles of more than two samples. The alignment of the
#' step functions is impacted by offset metadata assigned to the sequences.
#' When offset values are assigned one sequence if shifted horizontally to
#' align the start position 1 pointed to by the offset value for each sequence.
#' (see also parameter \code{markOffset}). If no offset metadata is available
#' for the sequences both step functions start at their first position on the
#' left side of the plot. The vertical plot range can be determined by the
#' \code{rng} argument. If the plot is generated for one profile, the sequence
#' is is visualized above the plot, for two sequences the first sequence is
#' shown above, the second sequence below the plot. Matching characters at a
#' position are shown in the same color (by default in \code{"black"}, the
#' non-matching characters in the sample-specific colors (see parameter
#' \code{col}). Annotation information can also be visualized along with the
#' step function. A call with two prediction profiles should facilitate the
#' comparison of profiles (e.g. wild type versus mutated sequence).
#'
#' The baseline for the step function of a single sample represents the offset
#' b of the model distributed equally to all sequence positions according to
#' the following reformulation of the discriminant function
#'
#'
#' \deqn{f(\vec{x}) = b + \sum_{i=1}^L(s_i(\vec{x})) = \sum_{i=1}^L(s_i(\vec{x})
#' - \frac{-b}{L})}{f(x) = b + sum(si(x)) = sum(si(x) -(-b/L))
#'   for i = 1, ... L}
#'
#'
#' For standardized plots (see parameter \code{standardize} this baseline value
#' is subtracted from the weight contribution at each position. When sequences
#' of different length are plotted together only a standardized plot gives
#' compareable y ranges for both step functions. For sequences of equal length
#' the visualization can be done in non-standardized or standardized form
#' showing the lightgray horizontal baseline at positon \eqn{y=-b/L} or at
#' \eqn{y=0}. If the area between the step function and the baseline lying
#' above the baseline is larger than the area below the baseline the sample
#' is predicted as belonging to the class assciated with positive values
#' of the discrimination function, otherwise to the opposite class. (For
#' multiclass problems prediction profiles can only be generated from the
#' feature weights related to one of the classifiers in the pairwise or
#' one-against-rest approaches leaving only two classes for the profile
#' plot.)
#'
#' When plotting to a pdf it is recommended to use a height to width ratio
#' of around 1:(max sequence length/25), e.g. for a maximum sequence length
#' of 500 bases or amino acids select height=10 and width=200 when opening the
#' pdf document for plotting.\cr\cr
#' @return see details above
#' @seealso \code{\link{getPredictionProfile}},
#' \code{\link{positionDependentKernel}}, \code{\link{mcols}},
#' \code{\link{spectrumKernel}}, \code{\link{mismatchKernel}},
#' \code{\link{gappyPairKernel}}, \code{\link{motifKernel}}
#' @examples
#'
#' ## set seed for random generator, included here only to make results
#' ## reproducable for this example
#' set.seed(456)
#' ## load transcription factor binding site data
#' data(TFBS)
#' enhancerFB
#' ## select 70% of the samples for training and the rest for test
#' train <- sample(1:length(enhancerFB), length(enhancerFB) * 0.7)
#' test <- c(1:length(enhancerFB))[-train]
#' ## create the kernel object for gappy pair kernel with normalization
#' gappy <- gappyPairKernel(k=1, m=3)
#' ## show details of kernel object
#' gappy
#'
#' ## run training with explicit representation
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=gappy,
#'                pkg="LiblineaR", svm="C-svc", cost=80, explicit="yes",
#'                featureWeights="yes")
#'
#' ## compute and plot ROC for test sequences
#' preddec <- predict(model, enhancerFB[test], predictionType="decision")
#' rocdata <- computeROCandAUC(preddec, yFB[test], allLabels=unique(yFB))
#' plot(rocdata)
#'
#' ## generate prediction profile for the first three test sequences
#' predProf <- getPredictionProfile(enhancerFB, gappy, featureWeights(model),
#'                                  modelOffset(model), sel=test[1:3])
#'
#' ## show prediction profiles
#' predProf
#'
#' ## plot prediction profile to pdf
#' ## As sequences are usually very long select a ratio of height to width
#' ## for the pdf which takes care of the maximum sequence length which is
#' ## plotted. Only single or pairs of prediction profiles can be plotted.
#' ## Plot profile for window size 1 (default) and 50. Load package Biobase
#' ## for openPDF
#' \dontrun{
#' library(Biobase)
#' pdf(file="PredictionProfile1_w1.pdf", height=10, width=200)
#' plot(predProf, sel=c(1,3))
#' dev.off()
#' openPDF("PredictionProfile1_w1.pdf")
#' pdf(file="PredictionProfile1_w50.pdf", height=10, width=200)
#' plot(predProf, sel=c(1,3), windowSize=50)
#' dev.off()
#' openPDF("PredictionProfile1_w50.pdf")
#' pdf(file="PredictionProfile2_w1.pdf", height=10, width=200)
#' plot(predProf, sel=c(2,3))
#' dev.off()
#' openPDF("PredictionProfile2_w1.pdf")
#' pdf(file="PredictionProfile2_w50.pdf", height=10, width=200)
#' plot(predProf, sel=c(2,3), windowSize=50)
#' dev.off()
#' openPDF("PredictionProfile2_w50.pdf")
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics} (accepted).
#' DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btv176}{10.1093/bioinformatics/btv176}.
#' @keywords prediction profile
#' @keywords plot
#' @keywords methods
#' @export

setMethod("plot", signature(x="PredictionProfile", y="missing"),
          plotPredictionProfile.Missing)

plot.cvResult <- function(x, col="springgreen")
{
    if (nrow(x@foldErrors) > 1)
    {
        boxplot(t(x@foldErrors), horizontal=TRUE, col=col,
                main="Cross Validation Error", xlab="Cross Validation Error")
    }
    else
    {
        cvErrors <- as.numeric(x@foldErrors)
        bar <- boxplot(cvErrors, plot=FALSE)
        boxplot(cvErrors, outline=FALSE, horizontal=TRUE, col=col,
                ylim=c(min(c(bar$stats,bar$out)), max(c(bar$stats,bar$out))),
                main="Cross Validation Error", xlab="Cross Validation Error")
    }
}

#' @rdname plot-methods
#' @aliases
#' plot,CrossValidationResult-method
#' @details
#' Plotting of CrossValidation Result\cr\cr
#' The second variant of \code{plot} method shown in the usage section
#' displays the cross validation result as boxplot.
#' @export
#'

setMethod("plot", signature(x="CrossValidationResult", y="missing"),
          plot.cvResult)

plot.performance <- function(x, sel=c("ACC", "BACC", "MCC", "AUC"))
{
    if (!is(x, "ModelSelectionResult"))
        stop("'x' is not of class 'ModelSelectionResult'")

    sel <- match.arg(sel)

    if (sel == "MCC")
    {
        title <- "Grid Matthews Correlation Coefficient"
        perfData <- x@gridMCC
    }
    else if (sel == "BACC")
    {
        title <- "Grid Balanced Accuracy"
        perfData <- x@gridBACC
    }
    else if (sel == "AUC")
    {
        if (length(x@gridAUC) == 0)
            stop("AUC was not stored in model selection result\n")
        
        title <- "Grid Area Under ROC Curve"
        perfData <- x@gridAUC
    }
    else
    {
        title <- "Grid Accuracy"
        perfData <- x@gridACC
    }

    numRows <- nrow(perfData)
    numCols <- ncol(perfData)
    xat <- c(-0.5, c(0:numCols))*(numCols + 1)/numCols^2
    yat <- c(-0.5, c(0:numRows))*(numRows + 1)/numRows^2
    xtext <- c("", paste("GC", 1:numCols, sep="_"), "")
    ytext <- c("", paste("K", 1:numRows, sep="_"), "")
    image(t(perfData), main=title, axes=FALSE, xlab="Grid Column",
          ylab="Kernel")
    axis(side=1, xat, labels=xtext, tick='n', las=2, cex.axis=0.7)
    axis(side=2, yat, labels=ytext, tick='n', las=2, cex.axis=0.7)
}

#' @rdname plot-methods
#' @aliases
#' plot,ModelSelectionResult-method
#' @details
#' Plotting of Grid Performance Values\cr\cr
#' The third variant of \code{plot} method shown in the usage section
#' plots grid performance data as grid with the color of each rectange
#' corresponding to the preformance value of the grid point.
#' @export

setMethod("plot", signature(x="ModelSelectionResult", y="missing"),
          plot.performance)

plot.roc <- function(x, lwd=2, aucDigits=3, cex=0.8, side=1,
                     line=-3, adj=0.9, ...)
{
    if (length(tpr(x)) < 3 || length(fpr(x)) < 3)
        stop("missing tpr/fpr data\n")

    plot(fpr(x), tpr(x), type="l", lwd=lwd,
         xlab="False Positive Rate", ylab="True Positive Rate", ...)

    if (aucDigits != 0 && !is.na(auc(x)))
    {
        mtext(paste("AUC:", round(auc(x), aucDigits)), side=side,
              line=line, adj=adj, cex=cex)
    }
}

#' @rdname plot-methods
#' @aliases
#' plot,ROCData-method
#' @details
#' Plotting of Receiver Operating Characteristics (ROC)\cr\cr
#' The fourth variant of \code{plot} method shown in the usage section
#' plots the receiver operating characteristics for the given ROC data.
#' @export

setMethod("plot", signature(x="ROCData", y="missing"),
          plot.roc)
