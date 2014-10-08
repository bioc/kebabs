##345678901234567890123456789012345678901234567890123456789012345678901234567890


#initKebabs <- function(lib, pkg)
#{
#    library.dynam(chname="kebabs",package=pkg,lib.loc=lib)

getSupportedPackages <- function()
{
    c("kernlab", "e1071", "LiblineaR")
}

initKebabs <- function()
{

    ## $$$ TODO extend classifier mapping table with psvm

    ## multiclass pairwise via C-svc, l2rl2l-svc, l2rl2lp-svc or l1rl2l-svc
    ## with more than two classes in label vector
    ## multiclass in LiblineaR is one against the rest
    ## default multiclass in kernlab and libsvm is pairwise
    classifierMap <- matrix(c(
        "C-svc",        "C-classification",        "3",
        "",             "",                        "1",
        "",             "",                        "2",
        "",             "",                        "5",
        "nu-svc",       "nu-classification",       "",
        "C-bsvc",       "",                        "",
        # mc pairwise - no special type needed
        # mc native - Crammer/Singer
        "spoc-svc",     "",                        "4",
        # mc native - Weston/Watkins
        "kbb-svc",      "",                        "",
        "one-svc",      "one-classification",      "",
        "eps-svr",      "eps-regression",          "",
        "nu-svr",       "nu-regression",           "",
        "eps-bsvr",     "",                        ""
        ), ncol=length(getSupportedPackages()), byrow=TRUE)

    colnames(classifierMap) <- getSupportedPackages()
    rownames(classifierMap) <- c(
        "C-svc",      # C-SVC, L2 regularization, L1 loss
        "l2rl2l-svc", # L2RL2L-SVC (dual), L2 regularization, L2 loss
        "l2rl2lp-svc",# L2RL2L-SVC (primal), L2 regularization, L2 loss
        "l1rl2l-svc", # L1RL2L-SVC, L1 regularization, L2 loss
        "nu-svc",     # nu-SVC
        "C-bsvc",     # bound constraint-SVC
        "mc-natC",    # Crammer, Singer native multiclass
        "mc-natW",    # Weston, Watkins native multiclass
        "one-svc",    # one-class SVC
        "eps-svr",    # epsilon SVR
        "nu-svr",     # nu SVR
        "eps-bsvr"    # bound-constraint SVR
                                 )

    svmParMap <- matrix(c(
        "C",              "cost",                    "cost",
        "nu",             "nu",                      "",
        "epsilon",        "epsilon",                 "",
        "class.weights",  "class.weights",           "wi",
        "tol",            "tolerance",               "epsilon",
        "cross",          "cross",                   "cross",
        # next 3 are for polynomial kernel
        "degree",         "degree",                  "1",
        "scale",          "gamma",                   "1",
        "offset",         "coef0",                   "1"
        ), ncol=length(getSupportedPackages()), byrow=TRUE)

    colnames(svmParMap) <- getSupportedPackages()
    rownames(svmParMap) <- c(
                    "cost",          # cost parameter of C-SVM
                    "nu",            # nu parameter of nu-SVM
                    "eps",           # epsilon parameter of eps-SVR, nu-SVR
                    "classWeights",  # weights for asymmetrical class size
                    "tolerance",     # tolerance of termination criterion
                    "cross",         # k-fold cross validation
                    "degree",        # kernel param of polynomial kernel
                    "scale",         # kernel param of polynomial kernel
                    "offset")        # kernel param of polynomial kernel

    modelParMap <- matrix(c(
        "SVindex",      "index",        "",
        "alphaindex",   "index",        "",
        "vanilladot",   "linear",       "",
        "coef",         "coefs",        "",
        "b",            "rho",          "",
        "nclass",       "nclasses",     "NbClass",
        "lev",          "levels",       "ClassNames",
        "",             "",             "W",
        "",             "",             "Bias",
        "prob.model",   "probA",        "",
        "prob.model",   "probB",        "",
        "prob.model",   "sigma",        "",
        "next1",        "next2",        "next3"
        ), ncol=length(getSupportedPackages()), byrow=TRUE)

    colnames(modelParMap) <- getSupportedPackages()
    rownames(modelParMap) <- c("svIndex", "alphaIndex", "linearKernel", 
                               "coef", "b", "numClasses", "classNames",
                               "weights", "bias", "probA", "probB", "sigma",
                               "next12")

    allowedSeqClasses <- c("DNAVector", "RNAVector", "AAVector",
                           "DNAStringSet", "RNAStringSet", "AAStringSet",
                           "DNAString", "RNAString", "AAString")

    allowedSeqSetClasses <- allowedSeqClasses[1:6]

    allowedSeqTypes <- c(1:9)
    names(allowedSeqTypes) <- c("DNAexact", "DNAiupac",
                                "RNAexact", "RNAiupac",
                                "AAexact", "AAiupac",
                                "unknown")


    kInfo <- new("KebabsInfo", supportedPkgs=getSupportedPackages(),
                 classifierMap=classifierMap, svmParMap=svmParMap,
                 modelParMap=modelParMap, kebabsDebug=FALSE,
                 allowedSeqClasses=allowedSeqClasses,
                 allowedSeqSetClasses=allowedSeqSetClasses,
                 allowedSeqTypes=allowedSeqTypes)

    return(kInfo)
}

kebabsInfo <- initKebabs()

cleanupKebabs <- function(lib)
{
    ## currently not used
}

svmSupportsOnlyExplicitRep <- function(pkg, svm)
{
    if (pkg %in% c("e1071","LiblineaR"))
        return(TRUE)
    else
        return(FALSE)
}

mapModelParamName <- function(package, param)
{
    if (exists("kebabsInfo"))
    {
        if (package %in% colnames(kebabsInfo@modelParMap))
        {
            if (param %in% rownames(kebabsInfo@modelParMap))
                return(kebabsInfo@modelParMap[param, package])

            stop(paste("Parameter", param, "not found in\n"),
                 "parameter map\n")
        }
        else
        {
            stop(paste("Package", package, "not found in\n"),
                 "parameter map\n")
        }
    }
    else
        stop("kebabsInfo is missing\n")
}

mapSVMParamName <- function(package, param, oldPackage=NULL)
{
    if (exists("kebabsInfo"))
    {
        if (package %in% colnames(kebabsInfo@svmParMap))
        {
            if (param == "C")
                param <- "cost"

            if (param %in% rownames(kebabsInfo@svmParMap))
                return(kebabsInfo@svmParMap[param, package])

            if (!is.null(oldPackage))
            {
                oldRow <- which(kebabsInfo@svmParMap[, oldPackage] == param)

                if (length(oldRow) == 1)
                    return(kebabsInfo@svmParMap[oldRow, package])
            }

            stop(paste("Parameter", param, "not found in\n"),
                 "parameter map\n")
        }
        else
        {
            stop(paste("Package", package, "not found in\n"),
                 "parameter map\n")
        }
    }
    else
        stop("kebabsInfo is missing\n")
}

getBioCharset <- function(x, exact)
{
    charsets <- c("ACGT",                      ## DNA_EXACT
                  "ACGTMRWSYKVHDBN-+",         ## DNA_IUPAC
                  "ACGU",                      ## RNA_EXACT
                  "ACGUMRWSYKVHDBN-+",         ## RNA_IUPAC
                  "ACDEFGHIKLMNPQRSTUVWY",     ## AA_EXACT
                  "ABCDEFGHIJKLMNPQRSTUVWXYZ") ## AA_IUPAC

    names(charsets) <- c("DNAexact", "DNAiupac",
                         "RNAexact", "RNAiupac",
                         "AAexact", "AAiupac")

    bioCharset <- switch(class(x),
        "DNAStringSet" = ifelse(exact,
                                kebabsInfo@allowedSeqTypes["DNAexact"],
                                kebabsInfo@allowedSeqTypes["DNAiupac"]),
        "RNAStringSet" = ifelse(exact,
                                kebabsInfo@allowedSeqTypes["RNAexact"],
                                kebabsInfo@allowedSeqTypes["RNAiupac"]),
        "AAStringSet"  = ifelse(exact,
                                kebabsInfo@allowedSeqTypes["AAexact"],
                                kebabsInfo@allowedSeqTypes["AAiupac"]),
        "DNAVector"       = ifelse(exact,
                                kebabsInfo@allowedSeqTypes["DNAexact"],
                                kebabsInfo@allowedSeqTypes["DNAiupac"]),
        "RNAVector"    = ifelse(exact,
                                kebabsInfo@allowedSeqTypes["RNAexact"],
                                kebabsInfo@allowedSeqTypes["RNAiupac"]),
        "AAVector"     = ifelse(exact,
                                kebabsInfo@allowedSeqTypes["AAexact"],
                                kebabsInfo@allowedSeqTypes["AAiupac"]),
        stop("please use XStringSet or BioVector classes for\n",
             "       sequence data\n")
    )

    list(charsets[bioCharset], bioCharset)
}

