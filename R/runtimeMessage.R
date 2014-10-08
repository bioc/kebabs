#2345678901234567890123456789012345678901234567890123456789012345678901234567890

generateRuntimeMessage <- function(model, x)
{
    ## $$$ TODO - soecific runtime warnings are not included in the
    ## initial release
    model@ctlInfo@runtimeWarning <- FALSE
    return(model@ctlInfo)

    if (!model@ctlInfo@runtimeWarning)
        return(model@ctlInfo)

    limitFeatSpaceDim <- 100000
    limitNumSamples <- 3000
    runtimeMsg <- ""

    if (getFeatureSpaceDimension(model@svmInfo@selKernel, x) >
        limitFeatSpaceDim || model@numSequences > limitNumSamples)
    {
        runtimeMsg <-
        paste("Long runtime could occur because of large feature space\n",
              "        and / or large number of samples\n\n")
    }

    if (nchar(runtimeMsg) > 0)
    {
        runtimeMsg <- paste(runtimeMsg,
        "        Processing in KeBABS can be interrupted with Ctrl-C at any\n",
        "        point. Please note that SVM processing running in other\n",
        "        packages sometimes does not take care of user interrupt\n",
        "        requests.\n")

        ## no warning but message because warning is output after the execution
        ##warning(runtimeMsg, call.=FALSE)
        message(runtimeMsg)
        model@ctlInfo@runtimeWarning <- FALSE
    }

    return(model@ctlInfo)
}
