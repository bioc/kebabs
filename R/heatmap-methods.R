#2345678901234567890123456789012345678901234567890123456789012345678901234567890
midDend.local <-
    function (x) if (is.null(mp <- attr(x, "midpoint"))) 0 else mp

memberDend.local <-
    function (x) if (is.null(r <- attr(x, "members"))) 1 else r

isLeaf.local <-
    function (x) (is.logical(L <- attr(x, "leaf"))) && L

midCacheDend.local <- function (x)
{
    stopifnot(inherits(x, "dendrogram"))
    
    setmid <- function(d)
    {
        if (isLeaf.local(d))
            return(d)
        
        k <- length(d)
        
        if (k < 1)
            stop("dendrogram node with non-positive #{branches}")
        
        r <- d
        midS <- 0
        
        for (j in 1:k)
        {
            r[[j]] <- unclass(setmid(d[[j]]))
            midS <- midS + midDend.local(r[[j]])
        }
        
        if (k == 2)
            attr(r, "midpoint") <- (memberDend.local(d[[1]]) + midS) / 2
        else
            attr(r, "midpoint") <- midDend.local(d)
        
        r
    }
    
    setmid(x)
}

revDend.local <- function (x)
{
    if (isLeaf.local(x))
        return(x)
    
    k <- length(x)
    
    if (k < 1)
        stop("dendrogram non-leaf node with non-positive #{branches}")
    
    r <- x
    
    for (j in 1:k)
        r[[j]] <- revDend.local(x[[k + 1 - j]])
    
    midCacheDend.local(r)
}

heatmap.predprof <- function(x, Rowv=TRUE, add.expr, margins=c(5, 5),
                             RowSideColors=NULL,
                             cexRow=max(min(35 / nrow(x@profiles), 1), 0.1),
                             cexCol=max(min(35 / ncol(x@profiles), 1), 0.1),
                             main=NULL, dendScale=1, barScale=1, startPos=1,
                             endPos=ncol(x@profiles), labels=NULL,
                             windowSize=1, ...)
{
    if (nrow(x@profiles) < 1)
        stop("'x' must contain prediction profiles\n")

    if (startPos < 1 || startPos >= ncol(x@profiles))
        stop ("'startPos' is out of range\n")

    if (endPos < 1 || endPos > ncol(x@profiles))
        stop ("'endPos' is out of range\n")

    if (endPos < startPos)
        stop("'startPos' is larger than 'endPos'\n")

    if (!is.numeric(windowSize) || length(windowSize) != 1 || windowSize < 1)
        stop("argument `windowSize' must be a single integer value\n")

    ## column subsetting of profile just for display
    
    numProfiles <- nrow(x@profiles)
    numPositions <- ncol(x@profiles)
    
    rDend <- NULL
    doRdend <- FALSE
    
    ## sliding window averaging
    if (windowSize != 1)
    {
        if (!is.infinite(windowSize))
        {
            windowSize <- as.integer(windowSize)
        
            if ((windowSize %% 2) == 0)
                windowSize <- windowSize + 1
        }

        avgProfiles <- x@profiles

        for (i in 1:numProfiles)
        {
            numPos <- width(x@sequences)[i]

            for (j in 1:numPos)
            {
                lower <- max(1, j - floor(windowSize/2))
                upper <- min(numPos, j + floor(windowSize/2))
                avgProfiles[i,j] <-
                    sum(x@profiles[i,lower:upper]) / (upper - lower)
            }
        }

        x@profiles <- avgProfiles
    }

    ## determine sample order
    if (is.null(Rowv) || Rowv == FALSE || is.na(Rowv))
        rowInd <- 1:numProfiles
    else if (is.numeric(Rowv))
        rowInd <- Rowv
    else if (Rowv == "random")
        rowInd <- sample(1:numProfiles, numProfiles)
    else if (Rowv == "decision")
        rowInd <- order(rowSums(x@profiles), decreasing=TRUE)
    else
    {
        ## perform hierarchical clustering of rows
        ## including the zeros after the sequence end
        hc <- hclust(dist(x@profiles))
        rowInd <- hc$order
        rDend <- as.dendrogram(hc)
        doRdend <- TRUE
    }
    
    if (length(names(x@sequences)) == length(x@sequences))
        labRow <- names(x@sequences)
    else
        labRow <- as.character(1:numProfiles)
    
    labCol <- colnames(x@profiles)[startPos:endPos]
    
    ## replace zeros after sequence end with NAs
    for (i in 1:numProfiles)
    {
        lastPlusOne <- width(x@sequences)[i] + 1
        
        if (lastPlusOne <= numPositions)
            x@profiles[i,lastPlusOne:numPositions] <- NA
    }
    
    ## perform layouting
    lmat <- matrix(c(NA,3,2,1), nrow=2, byrow=TRUE)
    lwid <- c(if (doRdend) dendScale else 0.05, 4)
    lhei <- c(0.05 + if (!is.null(main)) 0.3 else 0, 4)
    
    if (!is.null(labels))
    {
        if (!is.numeric(barScale) || length(barScale) != 1 ||
            barScale < 0.4 || barScale > 4)
        stop("'barScale' must be a single positive value not larger than 4\n")
        
        if (is.null(labels))
            labels <- rep(NA, numProfiles)
        else
        {
            if (length(table(labels)) != length(RowSideColors))
            {
                stop("Number of row side colors must match the\n",
                     "        number of different labels\n")
            }
            
            if (is.factor(labels))
                labels <- as.numeric(labels)
            else if (is.character(labels))
                labels <- as.numeric(factor(labels))
        }
        
        labels <- matrix(labels[rev(rowInd)], ncol=1)
        lmat <- rbind(c(lmat[1, ] + 1, lmat[1, 2] + 1), c(3, lmat[2, ]))
        lwid <- c(lwid[1L], 0.1 * barScale, lwid[2L])
    }
    
    lmat[is.na(lmat)] <- 0

    ## plot
    dev.hold()
    on.exit(dev.flush())
    op <- par(no.readonly=TRUE)
    on.exit(par(op), add=TRUE)
    layout(lmat, widths=lwid, heights=lhei, respect=TRUE)

    numPositions <- endPos - startPos + 1
    par(mar=c(margins[1], 0, if (!is.null(main)) 1 else 0, margins[2]))

    image(1:numPositions, 1:numProfiles,
          t(x@profiles[rev(rowInd), startPos:endPos]),
          xlim=(0.5 + c(0, numPositions)), ylim=(0.5 + c(0, numProfiles)),
          axes=FALSE, xlab="", ylab="", ...)

    if (cexCol > 0)
    {
        axis(1, 1:numPositions, labels=labCol, las=2, line=-0.5, tick=0,
             cex.axis=cexCol)
    }

    if (cexRow > 0)
    {
        axis(4, 1:numProfiles, labels=labRow[rev(rowInd)], las=2,
             line=-0.5, tick=0, cex.axis=cexRow)
    }

    if (!missing(add.expr))
        eval(substitute(add.expr))
        
    if (!is.null(labels))
    {
        par(mar=c(margins[1], 0, if (!is.null(main)) 1 else 0, 0.3))

        image(1, 1:numProfiles, t(labels), col=RowSideColors, axes=FALSE,
              xlab="", ylab="")
    }

    par(mar=c(margins[1], 0, if (!is.null(main)) 1 else 0, 0))

    if (doRdend)
    {
        plot(revDend.local(rDend), horiz=TRUE, axes=FALSE, yaxs="i",
             leaflab="none")
    }
    else
        frame()

    par(mar=c(0, 0, if (!is.null(main)) 1.5 else 0, margins[2]))

    if (!is.null(main))
    {
        frame()
        par(xpd=NA)
        title(main, cex.main=(1.5 * op[["cex.main"]]))
    }

    return(invisible(rDend))
}

#' @rdname heatmap-methods
#' @title Heatmap Methods
#'
#' @description Create a heat map of prediction profiles
#'
#' @param x prediction profile of class
#' \code{\linkS4class{PredictionProfile}}.
#'
#' @param Rowv determines the row order of the plot. When set to \code{TRUE}
#' the profile rows are clustered via hierarchical clustering and a row
#' dendrogram is plotted. When set to \code{FALSE, NA or NULL} the order is
#' corresponds to the order of the sequences in the profile. If this parameter
#' has a value of \code{random} rows are ordered randomly, for \code{decision}
#' the ordering is according to decreasing decision values. A user-defined
#' order can be specified through a numeric vector of indices. Default=TRUE
#'
#' @param startPos start sequence position. Together with the
#' parameter \code{endPos} a subset of sequence positions can be selected
#' for the heatmap. Default=1
#'
#' @param endPos end sequence position (see also \code{startPos}).
#' Default=maximum sequence length in the profile.
#'
#' @param labels a numeric vector, character vector or factor specifying
#' the labels for the sequences in the profile. If this parameter is
#' different from NULL the labels are plotted as side bar using the
#' colors specified in the parameter \code{RowSideColors}. Default=NULL
#'
#' @param RowSideColors a vector of color values specifying the colors for
#' the side bar. Default=NULL
#'
#' @param add.expr largely analogous to the standard
#' \code{\link[stats:heatmap]{heatmap}} function.
#'
#' @param margins largely analogous to the standard
#' \code{\link[stats:heatmap]{heatmap}} function. Default=c(5,5)
#'
#' @param main largely analogous to the standard
#' \code{\link[stats:heatmap]{heatmap}} function.
#'
#' @param cexRow largely analogous to the standard
#' \code{\link[stats:heatmap]{heatmap}} function. When set to 0 the row
#' labels are suppressed. Default=defined dependent on number of
#' profile rows
#'
#' @param cexCol largely analogous to the standard
#' \code{\link[stats:heatmap]{heatmap}} function. When set to 0 the column
#' labels are suppressed. Default=defined dependent on number of profile
#' columns
#'
#' @param dendScale factor scaling the width of the row dendrogram; values
#' have to be larger than 0 and not larger than 2. Default=1
#'
#' @param barScale factor scaling the width of the label color bar. Values
#' have to be larger than 0 and not larger than 4. Default=1
#'
#' @param windowSize numerical value specifying the window size of an optional
#' sliding window averaging of the prediction profiles. The value must be
#' larger than 0. Even values are changed internally to odd values by adding
#' 1. Default=1
#'
#' @param ... additional parameters which are passed to the \code{image}
#' method transparently.
#'
#' @details
#'
#' The \code{heatmap} function provides plotting of heatmaps from prediction
#' profiles with various possibilities for sample (=row) ordering (see
#' parameter \code{Rowv}). The heatmap is shown together with an optional
#' color sidebar showing the labels and an optional row cluster dendrogram
#' when hierarchical clustering defines the row order. For long sequences the
#' heatmap can be restricted to a subset of positions. Additionally smoothing
#' can be applied to the prediction profiles through sliding window averaging.
#' Through smoothing important regions can become better visible.
#'
#' @return
#' Invisibly, a cluster dendrogram.
#' @seealso \code{\link{getPredictionProfile}}
#' @examples
#'
#' ## load coiled coil data
#' data(CCoil)
#'
#' ## define annotation specific gappy pair kernel
#' gappya <- gappyPairKernel(k=1,m=11, annSpec=TRUE)
#'
#' ## train model
#' model <- kbsvm(x=ccseq, y=as.numeric(yCC), kernel=gappya,
#'                pkg="e1071", svm="C-svc", cost=15)
#'
#' ## generate prediction profiles
#' predProf <- getPredictionProfile(ccseq, gappya,
#'                        featureWeights(model), modelOffset(model))
#'
#' ## show prediction profiles
#' predProf
#'
#' \dontrun{
#' ## plot heatmap for the prediction profiles - random ordering of samples
#' heatmap(predProf, Rowv="random", main="Prediction Profiles", labels=yCC,
#' RowSideColors=c("blue", "red"), cexRow=0.15, cexCol=0.3)
#'
#' ## plot heatmap for the prediction profiles - ordering by decision values
#' heatmap(predProf, Rowv="decision", main="Prediction Profiles", labels=yCC,
#' RowSideColors=c("blue", "red"), cexRow=0.15, cexCol=0.3)
#'
#' ## plot heatmap for the prediction profiles - with hierarchical clustering
#' heatmap(predProf, Rowv=TRUE, main="Prediction Profiles", labels=yCC,
#' RowSideColors=c("blue", "red"), cexRow=0.15, cexCol=0.3)
#' }
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs/}\cr\cr
#' (Bodenhofer, 2009) -- U. Bodenhofer, K. Schwarzbauer, M. Ionescu and
#' S. Hochreiter. Modelling position specificity in sequence kernels by fuzzy
#' equivalence relations. \cr\cr
#' (Mahrenholz, 2011) -- C.C. Mahrenholz, I.G. Abfalter, U. Bodenhofer, R. Volkmer
#' and S. Hochreiter. Complex networks govern coiled-coil oligomerizations -
#' predicting and profiling by means of a machine learning approach.\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics}, 31(15):2574-2576, 2015.
#' DOI: \doi{10.1093/bioinformatics/btv176}.
#' @keywords prediction profiles
#' @export

#' @rdname heatmap-methods
#' @aliases
#' heatmap
#' heatmap,PredictionProfile-method
#' @export

setMethod("heatmap",
          signature=signature(x="PredictionProfile", y="missing"),
          heatmap.predprof)
