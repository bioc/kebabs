##345678901234567890123456789012345678901234567890123456789012345678901234567890

checkKBSVMParams <- function(x, y, sel, kernel, svm, pkg, explicit,
                       explicitType, featureType, featureWeights, weightLimit,
                       classWeights, probModel, features=NULL, cross, noCross,
                       groupBy, nestedCross, noNestedCross, perfParameters,
                       perfObjective, runtimeWarning, addArgs)
{
    ## pkg can still be "auto" for basic training and CV in this
    ## function and is set after this function in selectSVMMethod
    
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


    ## check and load SVM packages
    if (!(length(pkg) == 1 && pkg[1] == "auto"))
    {
        pkg <- checkSVMPackagePresence(pkg)
        loadSVMPackages(pkg)
        
        if (all(pkg == "kernlab"))
            model@ctlInfo@onlyDense <- TRUE
    }
    
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

        if (length(explicit) > 0 && (explicit == "yes"))
        {
            stop("explicit representation not available for call with\n",
                 "       kernel matrix\n", call.=FALSE)
        }

        if (any(pkg == "LiblineaR"))
            stop("kernel matrix via LiblineaR is not supported\n", call.=FALSE)

        model@numSequences <- nrow(x)
        
        if (featureWeights == "auto")
            featureWeights <- "no"
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
            if (length(setdiff(c("ACC","BACC","MCC","AUC"),
                               perfParameters)) == 0)
                perfParameters <- setdiff(perfParameters, "AUC")
            else
                stop("AUC is currently not supported for multiclass")
        }
    }

    if (!isSingleString(explicitType) ||
        !(explicitType %in% c("auto", "sparse", "dense")))
        stop("'explicitType' must be auto, sparse or dense\n",call.=FALSE)

    if (explicitType == "sparse" &&
        model@ctlInfo@onlyDense == TRUE)
        stop("package kernlab does not support sparse data\n",call.=FALSE)

    if (!isSingleString(featureWeights) ||
        !(featureWeights %in% c("auto", "yes", "no")))
        stop("'featureWeights' must be auto, yes or no\n", call.=FALSE)

    ## feature weights not supported for user-defined kernels
    if (anyUserDefinedKernel(kernel))
    {
        if (featureWeights == "yes")
        {
            stop("feature weights not supported for user-defined kernels\n",
                 call.=FALSE)
        }

        if (featureWeights == "auto")
            featureWeights <- "no"

        if (explicit == "auto")
            explicit <- "no"
    }

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

    if ((length(features) > 0) && !is.character(features))
        stop("'features' must be a character vector\n", call.=FALSE)

    if ((!is.null(kernel) && is.list(kernel)) ||
        length(svm) > 1 || length(pkg) > 1 ||
        (!is.null(addArgs) && any(lapply(addArgs, length) > 1)))
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

    ## no automatic package selection for grid search and model selection
    if (length(pkg) == 1 && pkg == "auto")
        stop("missing SVM package name(s)\n")

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

            for (i in 1:length(pkg))
            {
                if (tolower(pkg[i]) == "liblinear")
                    pkg[i] = "LiblineaR"
                else
                    pkg[i] = tolower(pkg[i])

                classifierType <- kebabsInfo@classifierMap[svm[i], pkg[i]]

                if (nchar(classifierType) < 1)
                    stop("invalid svm selected\n", call.=FALSE)
            }

            model@modelSelResult@gridCols <- data.frame(pkg=pkg, svm=svm)
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

detectInstalledPackages <- function(lookFor)
{
    ## check which of the packages are available
    detected <- lookFor[which(lookFor %in% installed.packages())]
    return(detected)
}

checkSVMPackagePresence <- function(pkgs)
{
    if (length(pkgs) == 1 && pkgs[1] == "auto")
        return(pkgs)
    
    ## neutralization of case sensitivity
    newPkgs <- tolower(pkgs)
    newPkgs[which(newPkgs == "liblinear")] <- "LiblineaR"
    neededPkgs <- unique(newPkgs)
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
    
    return(newPkgs)
}

loadSVMPackages <- function(pkgs)
{
    loadSuccess <- TRUE
    neededPkgs <- unique(pkgs)

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
    sampleLimitForKM <- 5000
    fastSVMLimitForER <- 1000
    denseERFeatureLimit <- 1000
    featureSpaceLimit <- 100000

    ## $$$ TODO Remove
    kebabsInfo@kebabsDebug <- FALSE
    
    ## support automatic package selection only for C-svc
    if (model@svmInfo@reqPackage == "auto" &&
        model@svmInfo@reqSVM != "C-svc")
        stop("please specific an SVM package\n")

    model@svmInfo@selKernel <- model@svmInfo@reqKernel
    model@svmInfo@selSVM <- model@svmInfo@reqSVM

    ## select method and package dependent on data size, capabilities of
    ## SVM, compatibility with requested parameters and knowledge of
    ## performance characteristics
    
    ## 1. Determine whether explicit representation or kernel matrix

    ## user requests ER or package needs ER
    if (model@svmInfo@reqExplicit == "yes" ||
        model@svmInfo@reqPackage %in% c("LiblineaR"))
    {
        model@svmInfo@selExplicit <- TRUE
        model@svmInfo@explicitKernel <- model@svmInfo@reqFeatureType
    }
    ## user requests KM or passes kernel matrix or kernel allows only KM
    else if (model@svmInfo@reqExplicit == "no" ||
             is(x, "KernelMatrix") ||
             ((is(x, "XstringSet") || is(x, "BioVector")) &&
              !supportsExplicitRep(model@svmInfo@selKernel)))
    {
        model@svmInfo@selExplicit <- FALSE
        model@svmInfo@explicitKernel <- "no"
        
        if (model@svmInfo@reqExplicit == "yes")
            stop("processing via explicit representation is not possible\n")
    }
    else if (model@svmInfo@reqExplicit == "auto")
    {
        ## all packages support ER but kernlab only dense
        if (!is(x, "KernelMatrix") &&
            (is(x, "ExplicitRepresentation") ||
             supportsExplicitRep(model@svmInfo@selKernel) &&
             (model@numSequences > sampleLimitForKM ||
              (!isUserDefined(model@svmInfo@selKernel) &&
               (is(x, "XStringSet") || is(x, "BioVector")) &&
               getFeatureSpaceDimension(model@svmInfo@selKernel, x) <
                  featureSpaceLimit))))
        {
            model@svmInfo@selExplicit <- TRUE
            model@svmInfo@explicitKernel <- model@svmInfo@reqFeatureType
        }
        else
        {
            model@svmInfo@selExplicit <- FALSE
            model@svmInfo@explicitKernel <- "no"
        }
    }

    if (model@svmInfo@selExplicit == TRUE)
        model@ctlInfo@selMethod <- "explicitRep"
    else
        model@ctlInfo@selMethod <- "KernelMatrix"

    ## 2. If ER determine dense or sparse

    if (model@svmInfo@selExplicit == TRUE)
    {
        model@ctlInfo@sparse <- TRUE

        if (model@svmInfo@reqExplicitType == "dense")
            model@ctlInfo@sparse <- FALSE
        else if (model@svmInfo@reqExplicitType == "auto")
        {
            if (model@svmInfo@reqPackage == "kernlab")
                model@ctlInfo@sparse <- FALSE
        }
    }

    ## 3. Select package if not predefined

    if (model@svmInfo@reqPackage != "auto")
    {
        model@svmInfo@selPackage <- model@svmInfo@reqPackage
        model@svmInfo@selSVM <- model@svmInfo@reqSVM
        model@svmInfo@selSVMPar <- model@svmInfo@reqSVMPar

        if (model@svmInfo@selExplicit == TRUE &&
            !supportsExplicitRep(model@svmInfo@selKernel))
            stop("kernel does not support an explicit representation\n")

        if (model@svmInfo@selExplicit == FALSE &&
            model@svmInfo@selPackage %in% c("LiblineaR"))
            stop("SVM does not support kernel matrices\n")
    }
    else
    {
        if (length(model@svmInfo@reqSVM) > 1 ||
            !(model@svmInfo@reqSVM %in%
              c("C-svc", "nu-svc", "eps-svr", "nu-svr")))
            stop("please specific an R package")

        ## probability model only available in kernlab and e1071
        if (model@svmInfo@selExplicit == FALSE)
        {
            model@svmInfo@selPackage <- "e1071"
            
            if (model@svmInfo@reqExplicitType == "auto")
                model@ctlInfo@sparse <- FALSE
        }
        else
        {
            if (model@numSequences > fastSVMLimitForER &&
                "LiblineaR" %in% model@svmInfo@availPackages &&
                !(model@svmInfo@probModel == TRUE))
            {
                model@svmInfo@selPackage <- "LiblineaR"

                if (model@svmInfo@reqExplicitType == "auto")
                    model@ctlInfo@sparse <- FALSE
            }
            else
            {
                model@svmInfo@selPackage <- "e1071"
            
                if (model@svmInfo@reqExplicitType == "auto")
                {
                    if (inherits(x, "ExplicitRepresentationDense") ||
                        (length(model@svmInfo@selKernel) == 1 &&
                        inherits(model@svmInfo@selKernel, "SequenceKernel") &&
                        !isUserDefined(model@svmInfo@selKernel) &&
                        (is(x, "XstringSet") || is(x, "BioVector")) &&
                        getFeatureSpaceDimension(model@svmInfo@selKernel, x) <
                        denseERFeatureLimit))
                        model@ctlInfo@sparse <- FALSE
                    else
                        model@ctlInfo@sparse <- TRUE
                }
            }
        }

        if (model@svmInfo@selPackage != model@svmInfo@reqPackage)
            model@svmInfo <- convertSVMParameters(model)
        else
            model@svmInfo@selSVMPar <- model@svmInfo@reqSVMPar
    }

    ## 4. Load package if a new one was selected
    
    if (model@svmInfo@reqPackage == "auto")
        loadSVMPackages(model@svmInfo@selPackage)

    model@ctlInfo@classification <-
        isClassification(model@svmInfo@selPackage,
                         model@svmInfo@selSVM)

    if (kebabsInfo@kebabsDebug == "TRUE")
    {
        cat("Classification: ", model@ctlInfo@classification, "\n")
        cat("Method:         ", model@ctlInfo@selMethod, "\n")
        cat("Sparse:         ", model@ctlInfo@sparse, "\n")
        cat("Explicit Type:  ", model@svmInfo@explicitKernel, "\n")
        cat("Package:        ", model@svmInfo@selPackage, "\n")
        cat("SVM:            ", model@svmInfo@selSVM, "\n")
        cat("SVM Parameters:\n")
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
