##345678901234567890123456789012345678901234567890123456789012345678901234567890

sortWith_LC_Collate_C <- function(x, decreasing=FALSE, ...)
{
    orig_LC_COLLATE <- Sys.getlocale(category="LC_COLLATE")
    Sys.setlocale(category="LC_COLLATE", locale="C")
    res <- base::sort(x=x, decreasing=decreasing, ...)
    Sys.setlocale(category="LC_COLLATE", locale=orig_LC_COLLATE)
    res
}

getSlotValue <- function(slotname, model)
{
    if (isS4(model))
        return(slot(model, slotname))

    else if (length(attr(model, which="class", exact=TRUE)) > 0 &&
             length(names(model)) > 0 && is.list(model))
    {
        parVal <- NULL
        parIndex <- which(names(model) == slotname)

        if (length(parIndex > 0))
            parVal <- model[[parIndex]]

        return(parVal)
    }
    else
        stop("wrong svm model object\n")
}

#' @rdname KBModelAccessors
#'
#' @param paramName unified name of an SVM model data element
#' @param model a KeBABS model
#' @param raw when set to \code{TRUE} the parameter value is delivered in
#'        exactly the way as it is stored in the SVM specific model, when
#'        set to \code{FALSE} it is delivered in unified format
#' @return \code{getSVMSlotValue}:
#' value of requested parameter in unified or native format dependent on
#' parameter \code{raw}.
#' @export
#'

getSVMSlotValue <- function(paramName, model, raw=FALSE)
{
    slotName <- mapModelParamName(model@svmInfo@selPackage, paramName)

    if (slotName == "")
    {
        if (paramName == "b" && model@svmInfo@selPackage == "LiblineaR")
        {
            if (model@svmModel$Bias)
            {
                parValue <-
                     as.numeric(model@svmModel$W[, ncol(model@svmModel$W)])

                if (!raw)
                    parValue <- formatSVMSlotValue(paramName, model, parValue)
            }
            else
                parValue <- rep(0, nrow(model@svmModel$W))

            return(parValue)
        }

        ## $$$ TODO Remove warning
        if (paramName != "svIndex")
        {
            warning(paste("SVM model parameter", paramName,
                          "not found in mapping\n"),
                    paste("for package", model@svmInfo@selPackage), "\n")
        }
        return(NULL)
    }

    if (isS4(model@svmModel))
    {
        parValue <- slot(model@svmModel, slotName)

        if (!raw)
            parValue <- formatSVMSlotValue(paramName, model, parValue)

        return(parValue)
    }
    else if (length(attr(model@svmModel, which="class", exact=TRUE)) > 0 &&
             length(names(model@svmModel)) > 0 && is.list(model@svmModel))
    {
        parValue <- NULL

        if (class(model@svmModel) == "LiblineaR")
        {
            if (paramName == "svIndex")
            {
                ## TODO $$$ Remove warning
                warning("Attempt to read SV indices from LiblineaR model\n")
                return(NULL)
            }

            if (paramName == "b")
            {
                parValue <- 0

                if (model@svmModel$Bias == TRUE)
                    parValue <- model@svmMode$W[length(model@svmModel$W)]
            }
            else
            {
                parIndex <- which(names(model@svmModel) == slotName)

                if (length(parIndex > 0))
                    parValue <- model@svmModel[[parIndex]]
            }
        }
        else
        {
            parIndex <- which(names(model@svmModel) == slotName)

            if (length(parIndex > 0))
                parValue <- model@svmModel[[parIndex]]
        }

        if (!raw)
            parValue <- formatSVMSlotValue(paramName, model, parValue)

        return(parValue)
    }
    else
        stop("wrong svm model object\n")
}

formatSVMSlotValue <- function(paramName, model, parValue)
{
    if (paramName == "svIndex" && class(model@svmModel) == "svm")
        parValue <- sortWith_LC_Collate_C(parValue)

    if (paramName == "classNames")
    {
        if (length(parValue) == 2)
            parValue <- sortWith_LC_Collate_C(parValue, decreasing=TRUE)
    }

    if (paramName == "weights" && class(model@svmModel) == "LiblineaR" &&
        model@ctlInfo@classification)
    {
        classNames <- getSVMSlotValue("classNames", model, raw=TRUE)

        if (length(classNames) == 2 &&
            (any(sortWith_LC_Collate_C(classNames, decreasing=TRUE) !=
                 classNames)))
            ## as b is last element of W also b is inverted
            parValue <- - parValue
    }

    if (paramName == "coef")
    {
        if (class(model@svmModel) == "svm")
        {
            if (model@ctlInfo@classification && model@svmModel$nclasses > 2)
            {
                numClasses <- model@svmModel$nclasses
                coefs <- matrix(0,nrow=nrow(model@svmModel$SV),
                                ncol=numClasses*(numClasses - 1)  / 2)
                colnames(coefs) <- rep("", ncol(coefs))
                start <- c(0, cumsum(model@svmModel$nSV[1:
                                     length(model@svmModel$nSV)-1])) + 1
                nSV_1 <- model@svmModel$nSV - 1

                ## retrieve SV indices sorted per pairwise SVM
                svIndexRaw <- getSVMSlotValue("svIndex", model, raw=TRUE)
                ## sort em
                svIndexSorted <- sort(svIndexRaw)
                ## generate index mapping from svm-wise to sorted
                indexMap <- sapply(svIndexRaw, function(x){which(svIndexSorted
                                   %in% x)})
                pair <- 1

                for (i in 1:(numClasses-1))
                {
                    for (j in (i+1):numClasses)
                    {
                        coefs[indexMap[start[i]:(start[i]+nSV_1[i])], pair] <-
                            parValue[start[i]:(start[i] + nSV_1[i]), j-1]
                        coefs[indexMap[start[j]:(start[j]+nSV_1[j])], pair] <-
                            parValue[start[j]:(start[j] + nSV_1[j]), i]

                        if (model@svmModel$label[i] < model@svmModel$label[j])
                        {
                            ## reverse colnames and change sign
                            colnames(coefs)[pair] <- paste(
                                model@svmModel$level[model@svmModel$label[i]],
                                model@svmModel$level[model@svmModel$label[j]],
                                sep="/")
                            coefs[,pair] <- - coefs[,pair]
                        }
                        else
                        {
                            colnames(coefs)[pair] <- paste(
                                model@svmModel$level[model@svmModel$label[j]],
                                model@svmModel$level[model@svmModel$label[i]],
                                sep="/")
                        }

                        pair <- pair + 1
                    }
                }

                ## sort columns for consistency with kernlab
                coefs <- coefs[, order(colnames(coefs))]

                if (length(rownames(model@svmModel$SV)) > 0)
                {
                    rownames(coefs)[indexMap[1:nrow(model@svmModel$SV)]] <-
                                        rownames(model@svmModel$SV)
                }

                return(coefs)
            }

            parValue <- parValue[sort(model@svmModel$index,
                                      index.return=TRUE)$ix,]

            if (model@svmModel$nclasses == 2)
            {
                if (model@svmModel$labels[1] == 1 &&
                    model@ctlInfo@classification)
                    parValue <- -parValue
            }

            if (model@svmModel$nclasses <= 2 || !model@ctlInfo@classification)
                parValue <- matrix(parValue, ncol=1)
        }
        else if (class(model@svmModel) == "ksvm")
        {
            if (model@ctlInfo@classification && model@svmModel@nclass > 2)
            {
                coefs <- matrix(0, nrow=length(model@svIndex),
                                ncol=length(parValue))
                colnames(coefs) <- rep("", ncol(coefs))

                if (model@ctlInfo@multiclassType == "pairwise")
                {
                    if (is.list(model@svmModel@xmatrix) &&
                        (length(rownames(model@svmModel@xmatrix[[1]])) > 0))
                        rownames(coefs) <- rep("", nrow(coefs))

                    for (i in 1:ncol(coefs))
                    {
                        coefs[which(model@svmModel@SVindex %in%
                                    model@svmModel@alphaindex[[i]]),i] <-
                        model@svmModel@coef[[i]]

                        if (is.list(model@svmModel@xmatrix) &&
                            (length(rownames(model@svmModel@xmatrix[[1]])) > 0))
                        {
                            rownames(coefs)[which(model@svmModel@SVindex %in%
                                model@svmModel@alphaindex[[i]])] <-
                                    rownames(model@svmModel@xmatrix[[i]])
                        }
                    }

                    colIndex <- 1

                    for (i in 1:(length(model@svmModel@lev) - 1))
                    {
                        for (j in (i+1):length(model@svmModel@lev))
                        {
                            colnames(coefs)[colIndex] <-
                            paste(model@svmModel@lev[i],
                                  model@svmModel@lev[j], sep="/")
                            colIndex <- colIndex + 1
                        }
                    }
                }
                else if (model@ctlInfo@multiclassType == "CrammerSinger")
                {
                    if (length(rownames(model@svmModel@xmatrix)) > 0)
                    {
                        rownames(coefs) <- rownames(
                            model@svmModel@xmatrix)[model@svmModel@SVindex]
                    }

                    for (i in 1:ncol(coefs))
                    {
                        coefs[which(model@svmModel@SVindex %in%
                            model@svmModel@alphaindex[[i]]),i] <-
                                model@svmModel@coef[[i]]
                    }

                    colnames(coefs) <- names(parValue)
                }
                else if (model@ctlInfo@multiclassType == "WestonWatkins")
                {
                    ## $$$ TODO implement coef for WestonWatkins
                    stop("coef currently not implemented for Weston/Watkins\n",
                         "        multiclass")
                }

                return(coefs)
            }

            if (is.list(parValue) &&
                model@ctlInfo@classification &&
                model@svmModel@nclass == 2)
                parValue <- parValue[[1]]

            if (model@svmModel@nclass <= 2 || !model@ctlInfo@classification)
                parValue <- matrix(parValue, ncol=1)
        }

        return(parValue)
    }

    if (paramName == "b")
    {
        if (class(model@svmModel) == "svm")
        {
            ## labels is the index into the level vector
            if (!(model@ctlInfo@classification &&
                  model@svmModel$nclasses == 2 &&
                  (model@svmModel$labels[1] == 1)))
            parValue <- - parValue
        }
        else if (class(model@svmModel) == "ksvm")
                parValue <- - parValue
        else if (class(model@svmModel) == "LiblineaR")
        {
            classNames <- getSVMSlotValue("classNames", model, raw=TRUE)

            if (length(classNames) == 2 &&
                model@ctlInfo@classification &&
                (any(sortWith_LC_Collate_C(classNames, decreasing=TRUE) !=
                     classNames)))
                parValue <- - parValue
        }

        return(parValue)
    }

    if (paramName == "probA" && class(model@svmModel) == "ksvm")
        parValue <- parValue$A

    if (paramName == "probB" && class(model@svmModel) == "ksvm")
        parValue <- parValue$B

    if (paramName == "sigma" && class(model@svmModel) == "ksvm")
        parValue <- parValue[[1]]

    return(parValue)
}

getMulticlassType <- function(model)
{
    if (model@numClasses == 2)
        return(character(0))

    svmType <- model@svmInfo@selSVM

    if (model@svmInfo@selPackage == "LiblineaR")
    {
        if (svmType %in% c("C-svc", "l1rl2l-svc", "l2rl2l-svc", "l2rl2lp-svc"))
            return("oneAgainstRest")
    }

    mcType <- c("pairwise", "pairwise", "CrammerSinger", "WestonWatkins")
    mcSVMs <- c("C-svc", "nu-svc", "mc-natC", "mc-natW")
    ind <- which(mcSVMs %in% svmType)

    if (length(ind) > 0)
        return(mcType[ind])
    else
        return(character(0))
}

isClassification <- function(package, svm)
{
    ## $$$ TODO include PSVM
    if (svm[1] %in%
        c("C-svc", "l2rl2l-svc", "l2rl2lp-svc", "l1rl2l-svc", "nu-svc",
          "C-bsvc", "mc-natC", "mc-natW", "one-svc"))
        return(TRUE)
    else
        return(FALSE)
}

