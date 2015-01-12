##345678901234567890123456789012345678901234567890123456789012345678901234567890

checkKBSVMParams <- function(x, y, sel, kernel, svm, pkg, explicit,
                       explicitType, featureType, featureWeights, weightLimit,
                       classWeights, probModel, features=NULL, cross, noCross,
                       groupBy, nestedCross, noNestedCross, perfParameters,
                       perfObjective, runtimeWarning, addArgs)
{
    ## allocate new model
    model <- new("KBModel")
    model@svmInfo <- new("SVMInformation")
    model@ctlInfo <- new("ControlInformation")

    ## detect installed svms
    model@svmInfo@availPackages <-
                       detectInstalledPackages(getSupportedPackages())

    ## basic parameter check
    if (!is.character(pkg))
        stop("'pkg' parameter must be a character vector\n", call.=FALSE)

    ## load SVM packages
    loadSVMPackages(pkg)

    if (!is.character(svm))
        stop("'svm' parameter must be a character vector\n", call.=FALSE)

    if (inherits(x, "ExplicitRepresentation"))
    {
        if (length(sel) > 0)
            stop("parameter 'sel' only allowed for sequences\n", call.=FALSE)

        if ((length(explicit) > 0) && (!isSingleString(explicit) ||
            !(explicit %in% c("auto", "yes", "no"))))
        {
            stop("'explicit' must be auto, yes or no\n",
                 call.=FALSE)
        }

        model@numSequences <- nrow(x)
    }
    else if (is(x, "KernelMatrix"))
    {
        if (length(sel) > 0)
            stop("parameter 'sel' only allowed for sequences\n", call.=FALSE)

        if (length(explicit) > 0 && !(explicit == "no"))
        {
            stop("explicit representation not available for call with\n",
                 "       kernel matrix\n", call.=FALSE)
        }

        if (any(pkg == "LiblineaR"))
            stop("kernel matrix via LiblineaR is not supported\n", call.=FALSE)

        if (any(pkg == "e1071"))
        {
            stop("kernel matrix via e1071 is currently not supported\n",
                 call.=FALSE)
        }
        model@numSequences <- nrow(x)
    }
    else
    {
        if(length(sel) > 0)
        {
            if (!is.numeric(sel) || length(sel) > length(x) ||
                max(abs(sel)) > length(x) || min(abs(sel)) < 1)
            {
                stop("'sel' must contain subset indices for 'x'\n",
                     call.=FALSE)
            }

            if (all(sel < 0))
            {
                sel1 <- c(1:length(x))[-sel]
                sel <- sel1
            }
            if (length(sel) != length(y))
            {
                stop("length of 'sel' does not match length of 'y'\n",
                     call.=FALSE)
            }

            model@numSequences <- length(sel)
        }
        else
            model@numSequences <- length(x)

        model@sel <- sel

        if (is.null(kernel))
            stop("'kernel' parameter is missing", call.=FALSE)

        if ((length(kernel) == 1) &&
            !is(kernel, "SequenceKernel") && !is(kernel, "stringkernel"))
        {
            stop("'kernel' must be a SequenceKernel or stringkernel\n",
                 call.=FALSE)
        }

        if (!isSingleString(explicit) ||
            !(explicit %in% c("auto", "yes", "no")))
        {
            stop("'explicit' must be auto, yes or no\n",
                 call.=FALSE)
        }

        if (explicit == "no")
        {
            if (any(pkg == "LiblineaR"))
            {
                stop("kernel matrix via LiblineaR is not supported\n",
                     call.=FALSE)
            }

            if (any(pkg == "e1071"))
            {
                stop("kernel matrix via e1071 is currently not supported\n",
                     call.=FALSE)
            }
        }
    }


    if (!isSingleString(featureType) ||
        !(featureType %in% c("linear", "quadratic")))
    {
        stop("'featureType' must be linear or quadratic\n",
             call.=FALSE)
    }

    if (!is.numeric(y) && !is.character(y) && !is.factor(y))
        stop("'y' must be a numeric vector, a character vector\n",
             "        or a factor\n", call.=FALSE)

    if (model@numSequences != length(y))
    {
        stop("number of samples does not match with length of 'y'\n",
             call.=FALSE)
    }

    model@y <- y
    model@classWeights <- classWeights

    if (isClassification(pkg, svm))
    {
        if (is.factor(y))
            model@levels <- levels(y)
        else
            model@levels <- sort(unique(y))
            
            ## $$$ TODO AUC for multiclass is currently not supported
            if (("AUC" %in% perfParameters) && length(model@levels) > 2)
            {
                if (length(setdiff(c("ACC","BACC","MCC","AUC"), perfParameters)) == 0)
                    perfParameters <- setdiff(perfParameters, "AUC")
                else
                    stop("AUC is currently not supported for multiclass")
            }
    }

    if (!isSingleString(explicitType) ||
        !(explicitType %in% c("auto", "sparse", "dense")))
        stop("'explicitType' must be auto, sparse or dense\n",call.=FALSE)

    if ((length(features) > 0) && !is.character(features))
        stop("'features' must be a character vector\n", call.=FALSE)

    if ((!(is.null(kernel) && is.list(kernel))) ||
        length(svm) > 1 || length(pkg) > 1 ||
        (!is.null(addArgs) && length(addArgs) > 0))
    {
        model <- checkModelSelParams(model=model, x=x, y=y, kernel=kernel,
                        svm=svm, pkg=pkg, explicit=explicit,
                        featureWeights=featureWeights, sel=sel,
                        features=features, weightLimit=weightLimit,
                        cross=cross, noCross=noCross, groupBy=groupBy,
                        nestedCross=nestedCross, noNestedCross=noNestedCross,
                        addArgs=addArgs, perfParameters=perfParameters,
                        perfObjective=perfObjective)
    }

    if (!isSingleString(featureWeights) ||
        !(featureWeights %in% c("auto", "yes", "no")))
        stop("'featureWeights' must be auto, yes or no\n", call.=FALSE)

    ## feature weights currently not supported for mc-natW (kbb-svc in kernlab)
    if ("mc-natW" %in% svm)
    {
        if (featureWeights == "yes")
        {
            stop("feature weights currently not supported for mc-natW\n",
                 call.=FALSE)
        }

        if (featureWeights == "auto")
            featureWeights <- "no"
    }

    if (!is.null(weightLimit) && !isSingleNumber(weightLimit))
        stop("'weightLimit' must be a number\n", call.=FALSE)

    if (!is.logical(probModel) || length(probModel) != 1)
        stop("'probModel' must be TRUE or FALSE\n", call.=FALSE)

    ## platt scaling could be implemented in kebabs but is problematic anyway
    if    (probModel && !(pkg == "kernlab" || pkg == "e1071"))
    {
        stop("probability model requires package kernlab or e1071\n",
             call.=FALSE)
    }

    if (!is.numeric(cross) || !is.numeric(noCross))
        stop("'cross' and 'noCross' must be numeric values\n", call.=FALSE)

    if (cross < -1 || cross == 1)
        stop("'cross' must be -1, 0 or larger than 1\n", call.=FALSE)

    if (noCross < 1)
        stop("'noCross' must be larger than 0\n", call.=FALSE)

    if (nestedCross < -1 || nestedCross == 1)
        stop("'nestedCross' must be -1, 0 or larger than 1\n", call.=FALSE)

    if (noNestedCross < 1)
        stop("'noNestedCross' must be larger than 0\n", call.=FALSE)

    if (!is.null(groupBy))
    {
        if (!(is.numeric(groupBy) || is.factor(groupBy)))
            stop("'groupBy' must be a numeric vector or a factor\n",
                 call.=FALSE)

        if ((length(sel) > 0 && length(groupBy) != length(sel)) ||
            (length(sel) == 0 && length(groupBy) != getNoOfElementsOfSeqRep(x)))
        {
            stop("length of 'groupBy' is not matching number of training\n",
                 "       samples\n", call.=FALSE)
        }

        if (cross == -1)
            stop("leave one out CV is not possible for 'groupBy'\n",
                 call.=FALSE)

        numGroups <- length(table(as.numeric(groupBy)))

        if (numGroups < cross)
        {
            stop("the number of groups must be larger than or equal to\n",
                 "        the number of folds\n", call.=FALSE)
        }

        if (nestedCross == -1 || nestedCross > 1)
        {
            if (numGroups < nestedCross)
            {
                stop("the number of groups must be larger than or equal to\n",
                     "        the number of folds\n", call.=FALSE)
            }
        }
    }

    model@ctlInfo@runtimeWarning <- runtimeWarning

    ## store user requests to model
    model@svmInfo@reqSVM <- svm
    model@svmInfo@reqPackage <- pkg
    model@svmInfo@reqKernel <- kernel
    model@svmInfo@reqExplicit <- explicit
    model@svmInfo@reqExplicitType <- explicitType
    model@svmInfo@reqFeatureType <- featureType
    model@svmInfo@featureWeights <- featureWeights
    model@svmInfo@weightLimit <- weightLimit
    model@svmInfo@probModel <- probModel

    if(!is.null(addArgs))
        model@svmInfo@reqSVMPar <- as.list(addArgs)

    if (length(model@classWeights > 0))
        model@svmInfo@reqSVMPar$classWeights <- model@classWeights

    if (model@svmInfo@reqExplicitType=="sparse" ||
        is(x, "ExplicitRepresentationSparse"))
        model@ctlInfo@sparse=TRUE

    return(model)
}

checkModelSelParams <- function(model, x, y, sel, features, kernel, svm, pkg,
                                explicit, featureWeights, weightLimit, cross,
                                noCross, groupBy, nestedCross, noNestedCross,
                                addArgs, perfParameters, perfObjective)
{
    modelSel <- FALSE

    if (is.null(model@modelSelResult))
    {
        if (!(perfObjective %in% perfParameters))
            perfParameters <- c(perfParameters, perfObjective)

        model@modelSelResult <- new("ModelSelectionResult")
        model@modelSelResult@cross <- cross
        model@modelSelResult@noCross <- noCross
        model@modelSelResult@groupBy <- groupBy
        model@modelSelResult@nestedCross <- nestedCross
        model@modelSelResult@noNestedCross <- noNestedCross
        model@modelSelResult@perfParameters <- perfParameters
        model@modelSelResult@perfObjective <- perfObjective

        if (!is.null(kernel))
        {
            if (is.list(kernel))
            {
                if (length(features) > 0)
                {
                    stop("feature subsetting is not possible with multiple\n",
                         "        kernels\n", call.=FALSE)
                }

                for (i in 1:length(kernel))
                {
                    if (!(is(kernel[[i]], "stringkernel") ||
                          is(kernel[[i]], "SequenceKernel")))
                    {
                        stop("'kernel' must be a SequenceKernel or ",
                             "stringkernel\n", call.=FALSE)
                    }
                }

                model@modelSelResult@gridRows <- kernel
                modelSel <- TRUE
            }
            else
                model@modelSelResult@gridRows <- list(kernel)

            names(model@modelSelResult@gridRows) <-
                paste("Kernel", 1:length(model@modelSelResult@gridRows),
                      sep="_")
        }

        if (length(pkg) > 1 || length(svm) > 1)
        {
            if (length(pkg) != length(svm))
            {
                stop("'pkg' and 'svm' must be character vectors of\n",
                     "        identical length\n", call.=FALSE)
            }

            modelSel <- TRUE
            count <- length(which(tolower(pkg) == "kernlab"))

            if (count == length(pkg))
            {
                if (length(model@svmInfo@reqExplicitType) > 0 &&
                    model@svmInfo@reqExplicitType == "sparse")
                {
                    stop("Package kernlab does not support sparse data\n",
                         call.=FALSE)
                }

                model@ctlInfo@onlyDense <- TRUE
            }
            else
                model@ctlInfo@onlyDense <- FALSE

            for (i in 1:length(pkg))
            {
                if (!(pkg[i] %in% model@svmInfo@availPackages))
                    stop("not installed or unsupported package\n", call.=FALSE)

                if (tolower(pkg[i]) == "liblinear")
                    pkg[i] = "LiblineaR"
                else
                    pkg[i] = tolower(pkg[i])

                classifierType <- kebabsInfo@classifierMap[svm[i], pkg[i]]

                if (nchar(classifierType) < 1)
                    stop("invalid svm selected\n", call.=FALSE)

                if (!require(pkg[i], character.only=TRUE))
                {
                    stop(paste("Could not load package", pkg[i]), "\n",
                         call.=FALSE)
                }
            }

            model@modelSelResult@gridCols <- data.frame(pkg=pkg, svm=svm)
        }
        else
        {
            if (!require(pkg[1], character.only=TRUE))
            {
                stop(paste("Could not load package", pkg[1]), "\n",
                     call.=FALSE)
            }
        }

        if(!is.null(addArgs) && length(addArgs) > 0)
        {
            paramNames <- c()
            paramList <- list()

            ## check parameter vectors
            for (i in 1:length(addArgs))
            {
                currName <- names(addArgs[i])
                currElem <- addArgs[[i]]

                if (currName %in% c("C", "cost", "nu", "eps"))
                {
                    if (is.numeric(currElem) && length(currElem) > 1)
                    {
                        paramNames <- c(paramNames, currName)
                        paramList[[currName]] <- currElem
                        tryCatch(paramList[[currName]] <- currElem,
                                 error=function(e) {stop(e)})
                    }
                }
            }

            if (length(paramList) == 1)
            {
                modelSel <- TRUE

                if (!(is.null(model@modelSelResult@gridCols)))
                {
                    if (nrow(model@modelSelResult@gridCols) !=
                        length(paramList[[1]]))
                    {
                        stop("inconsistent length of svm parameter\n",
                             call.=FALSE)
                    }

                    model@modelSelResult@gridCols[paramNames] <-
                                                           paramList[[1]]
                }
                else
                {
                    model@modelSelResult@gridCols <-
                                    data.frame(expand.grid(paramList))
                    colnames(model@modelSelResult@gridCols) <- paramNames
                }
            }
            else if (length(paramList) > 1)
            {
                modelSel <- TRUE

                if (!(is.null(model@modelSelResult@gridCols)))
                {
                    stop("only single svm for parameter expansion possible\n",
                         call.=FALSE)
                }

                model@modelSelResult@gridCols <-
                        data.frame(expand.grid(paramList))
                colnames(model@modelSelResult@gridCols) <- paramNames
            }
        }

        if (!(is.null(model@modelSelResult@gridCols)))
        {
            if ((colnames(model@modelSelResult@gridCols)[1] == "pkg") &&
                (ncol(model@modelSelResult@gridCols) < 3))
            {
                stop("missing or incomplete SVM parameters",
                     " for model selection\n",call.=FALSE)
            }

            rownames(model@modelSelResult@gridCols) <-
                paste("GridCol", 1: nrow(model@modelSelResult@gridCols),
                      sep="_")
        }


        if (!modelSel)
            model@modelSelResult <- NULL
        else
        {
            if (cross == 0)
            {
                warning("grid search without cross validation (cross=0)\n",
                        call.=FALSE)
            }
        }
    }

    return(model)
}

## $$$ TODO Remove ???
getSVMParams <- function(model, ...)
{
    ## $$$ TODO
    return(model@svmInfo)
}

detectInstalledPackages <- function(lookFor)
{
    ## check which of the packages are available
    detected <- lookFor[which(lookFor %in% installed.packages())]
    return(detected)
}

loadSVMPackages <- function(pkgs)
{
    ## neutralization of case sensitivity
    neededPkgs = tolower(unique(pkgs))
    neededPkgs[which(neededPkgs == "liblinear")] <- "LiblineaR"
    unsuppPkgs <- setdiff(neededPkgs, getSupportedPackages())

    if (length(unsuppPkgs) > 0)
    {
        stop("\npackage list contains following unsupported packages:",
             unsuppPkgs, call.=FALSE)
    }

    presentPkgs <- detectInstalledPackages(neededPkgs)

    if (length(presentPkgs) < length(neededPkgs))
    {
        stop("\nplease install following packages:", 
             setdiff(neededPkgs, presentPkgs))
    }

    loadSuccess <- TRUE
    for (i in 1:length(neededPkgs))
    {
        tryCatch({library(neededPkgs[i], character.only=TRUE,
                          warn.conflicts=FALSE)},
                warning =
                  function(w)
                  {
                    message(w)
                  },
                error =
                  function(e)
                  {
                    message(paste("package", neededPkgs[i], 
                                  "is not installed"))
                   loadSuccess <- FALSE
                  })
    }

    if (!loadSuccess)
        stop("\n", call.=FALSE)
}

selectSVMMethod <- function(model, supportedPackages, x)
{
    ## determine control information

    ## $$$ TODO

    kmSampleLimit = 99
    exRepLimit = 100

    ## $$$ TODO decide whether dense or sparse ER
    ## $$$ TODO for probModel only kernlab or e1071 is possible as
    ##          automatically selected package

    ## select svm and method dependent on data size, compatibility of
    ## requested paramters and knowledge of performance characteristics
    ## warning if high performance demand expected

    ## string kernel does not support explicit representation
    ## => always go via kernel matrix in this case

    ## decide whether explicit representation or kernel matrix
    if (model@svmInfo@reqExplicit == "yes" ||
        (model@svmInfo@reqExplicit=="auto" &&
         (!is(model@svmInfo@reqKernel, "stringkernel") &&
          !(is(model@svmInfo@reqKernel, "SequenceKernel") &&
            length(kernelParameters(model@svmInfo@reqKernel)$distWeight) > 0))
            && (model@numSequences > kmSampleLimit ||
                model@svmInfo@reqPackage %in% c("e1071", "LiblineaR"))))
    {
        model@svmInfo@selExplicit <- TRUE
        model@svmInfo@explicitKernel <- model@svmInfo@reqFeatureType
    }
    else
    {
        model@svmInfo@selExplicit <- FALSE
        model@svmInfo@explicitKernel <- "no"
    }

    if (model@svmInfo@selExplicit == TRUE)
        model@ctlInfo@selMethod <- "explicitRep"
    else
        model@ctlInfo@selMethod <- "KernelMatrix"

    if (kebabsInfo@kebabsDebug == "TRUE")
    {
        print(model@ctlInfo@selMethod)
        print(model@svmInfo@selExplicit)
        print(model@svmInfo@explicitKernel)
    }

    ## decide which package should be used and load it
    if (model@svmInfo@reqPackage != "auto")
    {
        if (!require(model@svmInfo@reqPackage, character.only=TRUE))
        {
            stop(paste("Could not load package", model@svmInfo@reqPackage,
                       "\n", call.=FALSE))
        }

        model@svmInfo@selPackage <- model@svmInfo@reqPackage
        model@svmInfo@selSVM <- model@svmInfo@reqSVM
        model@svmInfo@selKernel <- model@svmInfo@reqKernel
        model@svmInfo@selSVMPar <- model@svmInfo@reqSVMPar
    }
    else
    {
        if (model@svmInfo@selExplicit == TRUE &&
            model@numSequences > exRepLimit &&
            "LiblineaR" %in% model@svmInfo@availPackages)
            model@svmInfo@selPackage <- "LiblineaR"
        else
            model@svmInfo@selPackage <- "e1071"

        if (model@svmInfo@selPackage != model@svmInfo@reqPackage)
        {
            # model@svmInfo <- transSVMAndKernelPars(model)
            model@svmInfo <- convertSVMParameters(model)
        }
        else
        {
            model@svmInfo@selSVM <- model@svmInfo@reqSVM
            model@svmInfo@selKernel <- model@svmInfo@reqKernel
            model@svmInfo@selSVMPar <- model@svmInfo@reqSVMPar
        }
    }

    ## $$$ TODO extend sparsity control

    if (model@svmInfo@selExplicit == TRUE)
    {
        if (model@svmInfo@reqExplicitType == "auto")
        {
            if (model@svmInfo@selPackage == "kernlab")
                model@ctlInfo@sparse <- FALSE
            else
                model@ctlInfo@sparse <- TRUE
        }
        else if (model@svmInfo@reqExplicitType == "sparse")
            model@ctlInfo@sparse <- TRUE
        else
            model@ctlInfo@sparse <- FALSE
    }

#    if (model@svmInfo@selExplicit == TRUE &&
#        !model@svmInfo@reqExplicitType == "dense" &&
#        (is.null(model@modelSelResult) &&
#         model@svmInfo@selPackage != "kernlab") &&
#        !is.null(model@svmInfo@selKernel))
#    {
#        if ((class(model@svmInfo@selKernel) == "MotifKernel") ||
#            (kernelParameters(model@svmInfo@selKernel)$k >= 5))
#            model@ctlInfo@sparse=TRUE
#    }

    model@ctlInfo@classification <-
        isClassification(model@svmInfo@selPackage,
                         model@svmInfo@selSVM)

    if (kebabsInfo@kebabsDebug == "TRUE")
    {
        print(model@ctlInfo@selMethod)
        print(model@svmInfo@selExplicit)
        print(model@svmInfo@explicitKernel)
        print(model@svmInfo@selSVM)
        print(model@svmInfo@selPackage)
        print(model@svmInfo@selSVMPar)
    }

    return(list(model@svmInfo, model@ctlInfo))
}

addOrReplaceSVMParameters <- function(newParName, newPar, model)
{
    ## add parameter if not existing or replace
    ## only pure SVM parameters are treated
    ## kernel parameters are not relevant here

    if (length(model@svmInfo@reqSVMPar) == 0)
    {
        model@svmInfo@reqSVMPar <- list(newPar)
        names(model@svmInfo@reqSVMPar) <- newParName
        return(model@svmInfo)
    }

    namesPar <- names(model@svmInfo@reqSVMPar)

    if (length(namesPar) < 1)
        stop("missing name in grid parameter\n", call.=FALSE)

    parameterFound <- FALSE

    for (i in 1:length(model@svmInfo@reqSVMPar))
    {
        if (namesPar[i] == newParName)
        {
            model@svmInfo@reqSVMPar[[i]] <- newPar
            names(model@svmInfo@reqSVMPar)[i] <- newParName
            parameterFound <- TRUE
            break
        }
    }

    if (!parameterFound)
    {
        ## append to existing parameters
        model@svmInfo@reqSVMPar[length(model@svmInfo@reqSVMPar) + 1] <-
                                                                list(newPar)
        names(model@svmInfo@reqSVMPar)[length(model@svmInfo@reqSVMPar) + 1] <-
                                                                newParName
    }

    return(model@svmInfo)
}

convertSVMParameters <- function(model)
{
    ## TODO $$$ remove parameters that are not supported in target SVM

    ## adapt SVM parameters to target SVM
    ## only pure SVM parameters are treated
    ## kernel parameters are not relevant here

    oldPars <- model@svmInfo@reqSVMPar
    model@svmInfo@selSVMPar <- list()

    if (length(oldPars) == 0)
        return(model@svmInfo)

    for (i in 1:length(oldPars))
    {
        newParName <- mapSVMParamName(model@svmInfo@selPackage,
                                      names(oldPars)[i],
                                      model@svmInfo@reqPackage)
        if (length(newParName) > 0)
            model@svmInfo@selSVMPar[[newParName]] <- oldPars[[i]]
        else
            ## keep unknown parameters
            model@svmInfo@selSVMPar[[names(oldPars)[i]]] <- oldPars[[i]]
    }

    return(model@svmInfo)
}
