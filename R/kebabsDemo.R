##345678901234567890123456789012345678901234567890123456789012345678901234567890
#' @rdname kebabsOverview
#' @aliases
#' kebabs
#' KeBABS
#' KEBABS
#' @title kebabs\cr\cr
#'
#' @description KeBABS - An R package for kernel based analysis\cr
#' of biological sequences\cr\cr
#'
#' @details
#' Package Overview\cr\cr
#'
#' The package provides functionality for kernel based analysis of
#' DNA-, RNA- and amino acid sequences via SVM based methods. As core
#'functionality kebabs contains following sequence kernels: spectrum
#' kernel, mismatch kernel, gappy pair kernel and motif kernel. Apart
#' from an efficient implementation of position independent
#' functionality the kernels are extended in a novel way to take the
#' position of patterns into account for the similarity measure.
#' Because of the flexibility of the kernel formulation other kernels
#' like the weighted degree kernel or the shifted weighted degree
#' kernel are included as special cases. An annotation specific
#' variant of the kernels uses annotation information placed along
#' the sequence together with the patterns in the sequence. The
#' package allows generation of a kernel matrix or an explicit
#' representation for all available kernels which can be used with
#' methods implemented in other R packages. With focus on SVM based
#' methods kebabs provides a framework which simplifies the usage of
#' existing SVM implementations in kernlab, e1071 and LiblineaR.
#' Binary and multiclass classification as well as regression tasks
#' can be used in a unified way without having to deal with the
#' different functions, parameters and formats of the selected SVM.
#' As support for choosing hyperparameters the package provides cross
#' validation, grid search and model selection functions.For easier
#' biological interpretation of the results the package computes
#' feature weights for all SVMs and prediction profiles, which show
#' the contribution of individual sequence positions to the
#' prediction result and give an indication about the relevance of
#' sequence sections for the learning result and the underlying
#' biological functions.
#' @return
#' see above
#' @examples
#' ## load package provided sequence dataset
#' data(TFBS)
#'
#' ## display sequences
#' enhancerFB
#'
#' ## display part of label vector
#' head(yFB, 20)
#'
#' ## display no of samples of positive and negative class
#' table(yFB)
#'
#' ## split dataset into training and test samples
#' train <- sample(1:length(enhancerFB), 0.7*length(enhancerFB))
#' test <- c(1:length(enhancerFB))[-train]
#'
#' ## create the kernel object for the normalized spectrum kernel
#' spec <- spectrumKernel(k=5)
#'
#' ## train model
#' ## pass sequence subset, label subset, kernel object, the package and
#' ## svm which should be used for training together with the SVM parameters
#' model <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=spec,
#'                pkg="LiblineaR", svm="C-svc", cost=10)
#'
#' ## predict the test samples
#' pred <- predict(model, enhancerFB, sel=test)
#'
#' ## evaluate the prediction result
#' evaluatePrediction(pred, yFB[test], allLabels=unique(yFB))
#'
#' @author Johannes Palme <kebabs@@bioinf.jku.at>
#' @references
#' \url{http://www.bioinf.jku.at/software/kebabs/}\cr\cr
#' J. Palme, S. Hochreiter, and U. Bodenhofer (2015) KeBABS: an R package
#' for kernel-based analysis of biological sequences.
#' \emph{Bioinformatics}, 31(15):2574-2576, 2015.
#' DOI: \doi{10.1093/bioinformatics/btv176}.
#' @keywords kebabs


kebabsDemo <- function()
{
    ## nullify objects loaded from TFBS data set in order to avoid warnings about
    ## global variables without visible binding
    enhancerFB <- NULL
    yFB <- NULL
    ccseq <- NULL
    yCC <- NULL

    ## load sequence data and labels
    data(TFBS, envir=environment())

    ## show sequences
    cat("\nExample1: Classification of DNA sequence as enhancer region\n",
        "with binding of EP300/CREBBP\n\n")
    print(enhancerFB)

    ## show first 20 labels
    cat("\n## first 20 label values (1=Enhancer Region, -1=No Enhancer)\n\n")
    print(head(yFB, 20))

    cat("\n## number of samples for positive and negative class\n\n")
    print(table(yFB))

    readline("\nType <Return> \t to continue\n")

    ## split dataset randomly into 80% training data and 20% test data
    train <- sample(1:length(enhancerFB), 0.8 * length(enhancerFB))
    test <- c(1:length(enhancerFB))[-train]

    ## define spectrum kernel object
    specK5 <- spectrumKernel(k=5)

    cat("\n## training is performed with 80% of the sequences\n",
        "## selected randomly from the dataset with following kernel:\n",
        sep="")

    print(specK5)

    cat("\ntraining on LiblineaR C-SVM\n",
          "---------------------------\n", sep="")
    ## train model with C-SVM in package LiblineaR
    ttimeL <- system.time(
               modelL <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=specK5,
                              pkg="LiblineaR", svm="C-svc", cost=10))[3]

    cat("\ntraining time:", ttimeL, "sec\n")
    ## predict test sequences
    predL <- predict(modelL, enhancerFB[test])

    cat("\nprediction result:\n\n")
    ## evaluate prediction results
    print(evaluatePrediction(predL, yFB[test], allLabels=unique(yFB)))

    readline("\nType <Return> \t to continue\n")

    cat("\ntraining on libsvm C-SVM via e1071\n",
          "----------------------------------\n", sep="")
    ## train model with C-SVM
    ttimeE <- system.time(
               modelE <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=specK5,
                              pkg="e1071", svm="C-svc", cost=10))[3]

    cat("\ntraining time:", ttimeE, "sec\n")

    ## predict test sequences
    predE <- predict(modelE, enhancerFB[test])

    cat("\nprediction result:\n\n")
    ## evaluate prediction results
    print(evaluatePrediction(predE, yFB[test], allLabels=unique(yFB)))

    readline("\nType <Return> \t to continue\n")

    cat("\ntraining on kernlab C-SVM\n",
          "-------------------------\n", sep="")
    ## train model with C-SVM
    ttimeK <- system.time(
               modelK <- kbsvm(x=enhancerFB[train], y=yFB[train], kernel=specK5,
                         pkg="kernlab", svm="C-svc", cost=10, explicit="yes"))[3]

    cat("\ntraining time:", ttimeK, "sec\n")

    ## predict test sequences
    predK <- predict(modelK, enhancerFB[test])

    cat("\n## prediction result:\n\n")
    ## evaluate prediction results
    print(evaluatePrediction(predK, yFB[test], allLabels=unique(yFB)))

    readline("\nType <Return> \t to continue\n")

    ## load sequence data and labels
    data(CCoil, envir=environment())

    cat("\n---------------------------------------------------------------\n\n")

    ## show seqeuences
    cat("\nExample2: Computation of prediction profile for dimeric\n",
        "coil protein GCN4 wild type and trimeric mutation\n\n")
    print(ccseq)

    ## show first 20 labels
    cat("\n## first 20 label values (1=Trimer, -1=Dimer)\n\n")
    print(head(yCC, 20))

    cat("\n## number of samples for positive and negative class\n\n")
    print(table(yCC))

    readline("\nType <Return> \t to continue\n")

    ## define annotation specific gappy pair kernel object
    gappya <- gappyPairKernel(k=1, m=11, annSpec=TRUE)

    cat("\n## training is performed with all sequences of the\n",
        "## dataset with following annotation specific kernel:\n",
        sep="")

    print(gappya)

    cat("\ntraining on libsvm via e1071 C-SVM\n",
          "----------------------------------\n", sep="")
    ## train model with C-SVM in package e1071 on full data
    model <- kbsvm(x=ccseq, y=yCC, kernel=gappya,
                   pkg="e1071", svm="C-svc", cost=15)

    cat("\n## create histogram of feature weights\n\n")
    hist(featureWeights(model), breaks=30, main="Histogram of Feature Weights",
         xlab="Feature Weights")

    readline("\nType <Return> \t to continue\n")

    GCN4 <- AAStringSet(c("MKQLEDKVEELLSKNYHLENEVARLKKLV",
                          "MKQLEDKVEELLSKYYHTENEVARLKKLV"))
    names(GCN4) <- c("GCN4wt", "GCN_N16Y,L19T")
    annCharset <- annotationCharset(ccseq)
    annot <- c("abcdefgabcdefgabcdefgabcdefga",
               "abcdefgabcdefgabcdefgabcdefga")

    annotationMetadata(GCN4, annCharset=annCharset) <- annot

    print(GCN4)

    predProf <- getPredictionProfile(GCN4, gappya,
                                     featureWeights(model),
                                     modelOffset(model))
    cat("\n## show prediction profile\n\n")
    print(predProf)

    readline("\nType <Return> \t to continue\n")

    cat("\n## plot prediction profile\n\n")
    plot(x=predProf, sel=c(1,2), ylim=c(-0.4, 0.2), heptads=TRUE,
         annotate=TRUE)
}

