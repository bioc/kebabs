if (!isGeneric("as.KernelMatrix"))
{
    setGeneric("as.KernelMatrix",function(x, center = FALSE)
               standardGeneric("as.KernelMatrix"))
}

if (!isGeneric("as.KernelMatrix"))
{
    setGeneric("as.KernelMatrix",function(x, center = FALSE)
                                          standardGeneric("as.KernelMatrix"))
}

if (!isGeneric("kernelParameters"))
{
    setGeneric("kernelParameters", function(object, ...)
                                         standardGeneric("kernelParameters"))
}

if (!isGeneric("getFeatureSpaceDimension"))
{
    setGeneric("getFeatureSpaceDimension", function(kernel, ...)
                          standardGeneric("getFeatureSpaceDimension"))
}

if (!isGeneric("kbsvm"))
    setGeneric("kbsvm", function(x, ...) standardGeneric("kbsvm"))

if (!isGeneric("predict"))
    setGeneric("predict", function(object, ...) standardGeneric("predict"))

if (!isGeneric("trainSVM"))
    setGeneric("trainSVM", function(x, ...) standardGeneric("trainSVM"))

if (!isGeneric("predictSVM"))
    setGeneric("predictSVM", function(x, ...) standardGeneric("predictSVM"))

if (!isGeneric("performCrossValidation"))
{
    setGeneric("performCrossValidation",
           function(object, ...) standardGeneric("performCrossValidation"))
}

if (!isGeneric("profiles"))
    setGeneric("profiles", function(object, ...) standardGeneric("profiles"))

if (!isGeneric("baselines"))
    setGeneric("baselines", function(object, ...) standardGeneric("baselines"))

if (!isGeneric("sequences"))
    setGeneric("sequences", function(object, ...) standardGeneric("sequences"))

if (!isGeneric("getPredictionProfile"))
{
    setGeneric("getPredictionProfile",
           function(object, ...) standardGeneric("getPredictionProfile"))
}

if (!isGeneric("annotationMetadata"))
{
    setGeneric("annotationMetadata", function(x, ..., value)
                                 standardGeneric("annotationMetadata"))
}

if (!isGeneric("annotationMetadata<-"))
{
    setGeneric("annotationMetadata<-", function(x, ..., value)
                                 standardGeneric("annotationMetadata<-"))
}

if (!isGeneric("annotationCharset"))
{
    setGeneric("annotationCharset", function(x, ..., value)
                                 standardGeneric("annotationCharset"))
}

if (!isGeneric("positionMetadata"))
{
    setGeneric("positionMetadata", function(x, ..., value)
                                 standardGeneric("positionMetadata"))
}

if (!isGeneric("positionMetadata<-"))
{
    setGeneric("positionMetadata<-", function(x, ..., value)
                                 standardGeneric("positionMetadata<-"))
}

## performance method for model selection result
if (!isGeneric("gridRows"))
{
    if (is.function("gridRows"))
        fun <- gridRows
    else
        fun <- function(object, ...) standardGeneric("gridRows")

    setGeneric("gridRows", fun)
}

if (!isGeneric("gridColumns"))
{
    if (is.function("gridColumns"))
        fun <- gridColumns
    else
        fun <- function(object, ...) standardGeneric("gridColumns")

    setGeneric("gridColumns", fun)
}

if (!isGeneric("gridErrors"))
{
    if (is.function("gridErrors"))
        fun <- gridErrors
    else
        fun <- function(object, ...) standardGeneric("gridErrors")

    setGeneric("gridErrors", fun)
}

if (!isGeneric("performance"))
{
    if (is.function("performance"))
        fun <- performance
    else
        fun <- function(object, ...) standardGeneric("performance")

    setGeneric("performance", fun)
}

## accessors for KBModel
if (!isGeneric("numseq"))
{
    if (is.function("numseq"))
        fun <- numseq
    else
        fun <- function(object, ...) standardGeneric("numseq")

    setGeneric("numseq", fun)
}

if (!isGeneric("numseq<-"))
    setGeneric("numseq<-", function(x, value) standardGeneric("numseq<-"))

if (!isGeneric("selected"))
{
    if (is.function("selected"))
        fun <- selected
    else
        fun <- function(object, ...) standardGeneric("selected")

    setGeneric("selected", fun)
}

if (!isGeneric("selected<-"))
    setGeneric("selected<-", function(x, value) standardGeneric("selected<-"))

if (!isGeneric("modelOffset"))
{
    if (is.function("modelOffset"))
        fun <- modelOffset
    else
        fun <- function(object, ...) standardGeneric("modelOffset")

    setGeneric("modelOffset", fun)
}

if (!isGeneric("modelOffset<-"))
{
    setGeneric("modelOffset<-", function(x, value)
               standardGeneric("modelOffset<-"))
}

if (!isGeneric("featureWeights"))
{
    if (is.function("featureWeights"))
        fun <- featureWeights
    else
        fun <- function(object, ...) standardGeneric("featureWeights")

    setGeneric("featureWeights", fun)
}

if (!isGeneric("featureWeights<-"))
{
    setGeneric("featureWeights<-", function(x, value)
                               standardGeneric("featureWeights<-"))
}

if (!isGeneric("cvResult"))
{
    if (is.function("cvResult"))
        fun <- cvResult
    else
        fun <- function(object, ...) standardGeneric("cvResult")

    setGeneric("cvResult", fun)
}

if (!isGeneric("cvResult<-"))
{
    setGeneric("cvResult<-", function(x, value)
               standardGeneric("cvResult<-"))
}

if (!isGeneric("modelSelResult"))
{
    if (is.function("modelSelResult"))
        fun <- modelSelResult
    else
        fun <- function(object, ...) standardGeneric("modelSelResult")

    setGeneric("modelSelResult", fun)
}

if (!isGeneric("modelSelResult<-"))
{
    setGeneric("modelSelResult<-", function(x, value)
               standardGeneric("modelSelResult<-"))
}

if (!isGeneric("svmModel"))
{
    if (is.function("svmModel"))
        fun <- svmModel
    else
        fun <- function(object, ...) standardGeneric("svmModel")

    setGeneric("svmModel", fun)
}

if (!isGeneric("svmModel<-"))
{
    setGeneric("svmModel<-", function(x, value)
               standardGeneric("svmModel<-"))
}

if (!isGeneric("probabilityModel"))
{
    if (is.function("probabilityModel"))
        fun <- svmModel
    else
        fun <- function(object, ...) standardGeneric("probabilityModel")

    setGeneric("probabilityModel", fun)
}

if (!isGeneric("probabilityModel<-"))
{
    setGeneric("probabilityModel<-", function(x, value)
               standardGeneric("probabilityModel<-"))
}

## accessors for PredictionProfile
if (!isGeneric("profiles"))
{
    if (is.function("profiles"))
        fun <- profiles
    else
        fun <- function(object, ...) standardGeneric("profiles")

    setGeneric("profiles", fun)
}

if (!isGeneric("profiles<-"))
    setGeneric("profiles<-", function(x, value) standardGeneric("profiles<-"))
